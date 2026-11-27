#pragma once

#include "psrio/detail/exceptions.hpp"

#include <algorithm>
#include <cerrno>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <deque>
#include <exception>
#include <fcntl.h>
#include <filesystem>
#include <fmt/format.h>
#include <mutex>
#include <span>
#include <sys/types.h>
#include <thread>
#include <unistd.h>
#include <utility>
#include <vector>

namespace psrio::detail {

/// Sequential regular-file writer backed by POSIX `write()` with large
/// buffering, optional background write thread, and optional Linux page-cache
/// eviction (set_drop_cache()) for cgroup-limited HPC jobs.
class WriteFile {
public:
    static constexpr std::size_t kBufferBytes = std::size_t{16} << 20U;

    WriteFile() noexcept = default;

    /// Create or truncate @p path for sequential writing.
    /// @throws psrio::IoError if the file cannot be created.
    explicit WriteFile(std::filesystem::path path,
                       std::size_t buffer_bytes = kBufferBytes)
        : m_buffer_limit(std::max<std::size_t>(buffer_bytes, 4096U)),
          m_path(std::move(path)) {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg)
        m_fd = ::open(m_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC,
                      0666);
        if (m_fd < 0) {
            error_check::throw_errno_io("open", m_path);
        }
    }

    WriteFile(const WriteFile&)            = delete;
    WriteFile& operator=(const WriteFile&) = delete;

    WriteFile(WriteFile&& other) noexcept { move_from(other); }

    WriteFile& operator=(WriteFile&& other) noexcept {
        if (this != &other) {
            try {
                close();
            } catch (...) {
                // Move assignment must not throw; close() reports errors.
            }
            move_from(other);
        }
        return *this;
    }

    ~WriteFile() noexcept {
        try {
            close();
        } catch (...) {
            // Destructors must not throw; call close() to observe failures.
        }
    }

    [[nodiscard]] bool is_open() const noexcept { return m_fd >= 0; }
    [[nodiscard]] const std::filesystem::path& path() const noexcept {
        return m_path;
    }

    void set_drop_cache(bool enable) noexcept { m_drop_cache = enable; }
    [[nodiscard]] bool drop_cache() const noexcept { return m_drop_cache; }

    /// Overlap file writes with caller computation on a dedicated thread.
    /// Safe to call right after construction before writing payload blocks.
    void set_async(bool enable, std::size_t queue_depth = 3U) {
        if (enable == m_async) {
            return;
        }
        if (enable) {
            m_depth  = std::max<std::size_t>(queue_depth, 1U);
            m_async  = true;
            m_stop   = false;
            m_error  = nullptr;
            m_thread = std::thread([this] { worker(); });
        } else {
            stop_worker();
        }
    }

    /// Append @p bytes sequentially.
    /// @throws psrio::IoError if closed or if a write fails.
    void write(std::span<const std::byte> bytes) {
        if (m_fd < 0) {
            throw IoError(
                fmt::format("psrio: writer for {} is closed", m_path.string()));
        }
        if (bytes.empty()) {
            return;
        }
        if (m_async) {
            write_async(bytes);
            return;
        }
        if (m_pending.empty() && bytes.size() >= m_buffer_limit) {
            write_all(bytes);
            return;
        }
        if (m_pending.size() + bytes.size() > m_buffer_limit &&
            !m_pending.empty()) {
            write_all(m_pending);
            m_pending.clear();
        }
        if (bytes.size() >= m_buffer_limit) {
            write_all(bytes);
            return;
        }
        m_pending.insert(m_pending.end(), bytes.begin(), bytes.end());
    }

    /// Flush all buffered and queued bytes to the OS file descriptor.
    /// @throws psrio::IoError if closed or if a write fails.
    void flush() {
        if (m_fd < 0) {
            throw IoError(
                fmt::format("psrio: writer for {} is closed", m_path.string()));
        }
        if (m_async) {
            submit_pending();
            std::unique_lock lock(m_mutex);
            m_space_cv.wait(lock, [this] {
                return (m_queue.empty() && !m_writing) || m_error;
            });
            rethrow_locked();
            return;
        }
        if (!m_pending.empty()) {
            write_all(m_pending);
            m_pending.clear();
        }
    }

    /// Flush, stop any background thread, and close the file descriptor.
    /// @throws psrio::IoError if flushing or closing fails.
    void close() {
        if (m_fd < 0) {
            return;
        }
        std::exception_ptr flush_err;
        try {
            if (m_async) {
                submit_pending();
                stop_worker();
            } else if (!m_pending.empty()) {
                write_all(m_pending);
                m_pending.clear();
            }
            evict_cache_final();
        } catch (...) {
            flush_err = std::current_exception();
        }
        const int fd    = std::exchange(m_fd, -1);
        const int rc    = ::close(fd);
        const int saved = errno;
        m_pending.clear();
        if (flush_err) {
            std::rethrow_exception(flush_err);
        }
        if (rc != 0) {
            error_check::throw_errno_io(saved, "close", m_path);
        }
    }

private:
    void move_from(WriteFile& other) noexcept {
        if (other.m_async) {
            try {
                other.submit_pending();
                other.stop_worker();
            } catch (...) {
                // Best-effort drain on move.
            }
        }
        m_fd             = std::exchange(other.m_fd, -1);
        m_drop_cache     = other.m_drop_cache;
        m_buffer_limit   = other.m_buffer_limit;
        m_written        = std::exchange(other.m_written, 0U);
        m_prev_evict_end = std::exchange(other.m_prev_evict_end, 0U);
        m_path           = std::move(other.m_path);
        m_pending        = std::move(other.m_pending);
    }

    void write_async(std::span<const std::byte> bytes) {
        if (!m_pending.empty() &&
            m_pending.size() + bytes.size() > m_buffer_limit) {
            submit_pending();
        }
        if (m_pending.empty() && bytes.size() >= m_buffer_limit) {
            std::vector<std::byte> chunk = take_spare(bytes.size());
            chunk.assign(bytes.begin(), bytes.end());
            enqueue(std::move(chunk));
            return;
        }
        m_pending.insert(m_pending.end(), bytes.begin(), bytes.end());
        if (m_pending.size() >= m_buffer_limit) {
            submit_pending();
        }
    }

    std::vector<std::byte> take_spare(std::size_t reserve_bytes) {
        std::vector<std::byte> buf;
        {
            const std::scoped_lock lock(m_mutex);
            if (!m_spare.empty()) {
                buf = std::move(m_spare.back());
                m_spare.pop_back();
            }
        }
        if (buf.capacity() < reserve_bytes) {
            buf.reserve(reserve_bytes);
        }
        buf.clear();
        return buf;
    }

    void submit_pending() {
        if (m_pending.empty()) {
            return;
        }
        std::vector<std::byte> chunk = std::move(m_pending);
        m_pending                    = take_spare(m_buffer_limit);
        enqueue(std::move(chunk));
    }

    void enqueue(std::vector<std::byte> chunk) {
        std::unique_lock lock(m_mutex);
        m_space_cv.wait(lock,
                        [this] { return m_queue.size() < m_depth || m_error; });
        rethrow_locked();
        m_queue.push_back(std::move(chunk));
        lock.unlock();
        m_work_cv.notify_one();
    }

    void stop_worker() {
        if (!m_async) {
            return;
        }
        {
            const std::scoped_lock lock(m_mutex);
            m_stop = true;
        }
        m_work_cv.notify_all();
        if (m_thread.joinable()) {
            m_thread.join();
        }
        m_async = false;
        m_spare.clear();
        if (m_error) {
            const auto err = std::exchange(m_error, nullptr);
            std::rethrow_exception(err);
        }
    }

    void rethrow_locked() {
        if (m_error) {
            const auto err = std::exchange(m_error, nullptr);
            std::rethrow_exception(err);
        }
    }

    void worker() {
        while (true) {
            std::vector<std::byte> chunk;
            {
                std::unique_lock lock(m_mutex);
                m_work_cv.wait(lock,
                               [this] { return !m_queue.empty() || m_stop; });
                if (m_queue.empty()) {
                    return;
                }
                chunk = std::move(m_queue.front());
                m_queue.pop_front();
                m_writing = true;
            }
            m_space_cv.notify_one();
            try {
                write_all(chunk);
            } catch (...) {
                const std::scoped_lock lock(m_mutex);
                m_error   = std::current_exception();
                m_writing = false;
                m_queue.clear();
                m_space_cv.notify_all();
                return;
            }
            {
                const std::scoped_lock lock(m_mutex);
                m_writing = false;
                if (m_spare.size() < m_depth + 1U) {
                    chunk.clear();
                    m_spare.push_back(std::move(chunk));
                }
            }
            m_space_cv.notify_all();
        }
    }

    void write_all(std::span<const std::byte> bytes) {
        const auto start_offset = m_written;
        std::size_t done        = 0;
        while (done < bytes.size()) {
            const auto want =
                std::min<std::size_t>(bytes.size() - done, kMaxRequestBytes);
            const ssize_t got = ::write(m_fd, bytes.data() + done, want);
            if (got < 0) {
                if (errno == EINTR) {
                    continue;
                }
                error_check::throw_errno_io("write", m_path);
            }
            if (got == 0) {
                throw IoError(
                    fmt::format("psrio: {}: short zero-byte write at byte {}",
                                m_path.string(), m_written));
            }
            done += static_cast<std::size_t>(got);
            m_written += static_cast<std::uint64_t>(got);
        }
        evict_cache_window(start_offset, bytes.size());
    }

    // NOLINTNEXTLINE(readability-make-member-function-const): Linux mutates
    void evict_cache_window([[maybe_unused]] std::uint64_t start,
                            [[maybe_unused]] std::uint64_t nbytes) noexcept {
        if (!m_drop_cache || m_fd < 0 || nbytes == 0U) {
            return;
        }
#if defined(__linux__) && defined(SYNC_FILE_RANGE_WRITE)
        // Start async writeback on the just-written window, then evict the
        // previous window once its writeback has had time to complete.
        (void)::sync_file_range(m_fd, static_cast<off_t>(start),
                                static_cast<off_t>(nbytes),
                                SYNC_FILE_RANGE_WRITE);
        if (start > m_prev_evict_end) {
            const auto prev_len = start - m_prev_evict_end;
            (void)::sync_file_range(m_fd, static_cast<off_t>(m_prev_evict_end),
                                    static_cast<off_t>(prev_len),
                                    SYNC_FILE_RANGE_WAIT_BEFORE |
                                        SYNC_FILE_RANGE_WRITE |
                                        SYNC_FILE_RANGE_WAIT_AFTER);
#ifdef POSIX_FADV_DONTNEED
            (void)::posix_fadvise(m_fd, static_cast<off_t>(m_prev_evict_end),
                                  static_cast<off_t>(prev_len),
                                  POSIX_FADV_DONTNEED);
#endif
            m_prev_evict_end = start;
        }
#elif defined(POSIX_FADV_DONTNEED)
        if (start > m_prev_evict_end) {
            (void)::posix_fadvise(m_fd, static_cast<off_t>(m_prev_evict_end),
                                  static_cast<off_t>(start - m_prev_evict_end),
                                  POSIX_FADV_DONTNEED);
            m_prev_evict_end = start;
        }
#endif
    }

    void evict_cache_final() noexcept {
        if (!m_drop_cache || m_fd < 0 || m_written <= m_prev_evict_end) {
            return;
        }
#ifdef POSIX_FADV_DONTNEED
        (void)::posix_fadvise(m_fd, static_cast<off_t>(m_prev_evict_end),
                              static_cast<off_t>(m_written - m_prev_evict_end),
                              POSIX_FADV_DONTNEED);
#endif
        m_prev_evict_end = m_written;
    }

    static constexpr std::size_t kMaxRequestBytes = std::size_t{1} << 30U;

    int m_fd{-1};
    bool m_drop_cache{false};
    std::size_t m_buffer_limit{kBufferBytes};
    std::uint64_t m_written{0};
    std::uint64_t m_prev_evict_end{0};
    std::filesystem::path m_path;
    std::vector<std::byte> m_pending;

    bool m_async{false};
    std::size_t m_depth{3};
    std::mutex m_mutex;
    std::condition_variable m_work_cv;
    std::condition_variable m_space_cv;
    std::deque<std::vector<std::byte>> m_queue;
    std::vector<std::vector<std::byte>> m_spare;
    std::exception_ptr m_error;
    bool m_writing{false};
    bool m_stop{false};
    std::thread m_thread;
};

} // namespace psrio::detail
