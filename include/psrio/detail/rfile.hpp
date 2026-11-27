#pragma once

#include "psrio/detail/exceptions.hpp"

#include <algorithm>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <fmt/format.h>
#include <limits>
#include <span>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <utility>

namespace psrio::detail {

/// Access pattern hint given to the kernel when a file is opened.
enum class Access : std::uint8_t {
    /// Mostly forward reads: the kernel reads further ahead.
    kSequential,
    /// Default read-ahead, for scattered small reads.
    kNormal,
};

/// Read-only regular file read with positioned reads (`pread`).
///
/// Streaming a large file through `pread` with large requests reads networked
/// and parallel file systems (BeeGFS, Lustre, GPFS) several times faster than
/// page faults on a memory mapping, and keeps the resident set to the
/// caller's buffers. The file is opened with a sequential access hint
/// unless Access::kNormal is asked for.
///
/// With drop-cache enabled via set_drop_cache(true), release() asks the kernel
/// to drop the cached pages of a range that has been consumed, for jobs whose
/// memory limit counts the page cache (Linux).
class ReadFile {
public:
    ReadFile() noexcept = default;

    /// @throws psrio::IoError if @p path cannot be opened or is not a regular
    ///         file.
    explicit ReadFile(std::filesystem::path path,
                      Access access = Access::kSequential)
        : m_access(access),
          m_path(std::move(path)) {
        open_fd(true);
    }

    ~ReadFile() noexcept { reset(); }

    ReadFile(const ReadFile&)            = delete;
    ReadFile& operator=(const ReadFile&) = delete;

    ReadFile(ReadFile&& other) noexcept
        : m_fd{std::exchange(other.m_fd, -1)},
          m_size{std::exchange(other.m_size, 0)},
          m_drop_cache{other.m_drop_cache},
          m_access{other.m_access},
          m_path{std::move(other.m_path)} {}

    ReadFile& operator=(ReadFile&& other) noexcept {
        if (this != &other) {
            reset();
            m_fd         = std::exchange(other.m_fd, -1);
            m_size       = std::exchange(other.m_size, 0);
            m_drop_cache = other.m_drop_cache;
            m_access     = other.m_access;
            m_path       = std::move(other.m_path);
        }
        return *this;
    }

    [[nodiscard]] bool is_open() const noexcept { return m_fd >= 0; }

    /// Close the OS file descriptor while retaining path, size, and options
    /// so ensure_open() can reopen it on demand.
    void close_fd() noexcept {
        if (m_fd >= 0) {
            (void)::close(m_fd);
            m_fd = -1;
        }
    }

    /// Reopen the file descriptor if it was closed via close_fd().
    /// @throws psrio::IoError if the file cannot be reopened or changed size.
    void ensure_open() {
        if (m_fd >= 0) {
            return;
        }
        if (m_path.empty()) {
            throw IoError("psrio: cannot reopen an uninitialized ReadFile");
        }
        open_fd(false);
    }

    /// A second descriptor for the same file, so another thread can read it
    /// while this one is moved or closed.
    /// @throws psrio::IoError if the descriptor cannot be duplicated.
    [[nodiscard]] ReadFile duplicate() const {
        ReadFile copy;
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg)
        copy.m_fd = ::fcntl(m_fd, F_DUPFD_CLOEXEC, 0);
        if (copy.m_fd < 0) {
            error_check::throw_errno_io("dup", m_path);
        }
        copy.m_size       = m_size;
        copy.m_drop_cache = m_drop_cache;
        copy.m_access     = m_access;
        copy.m_path       = m_path;
        return copy;
    }

    /// File size in bytes at open.
    [[nodiscard]] std::uint64_t size() const noexcept { return m_size; }

    [[nodiscard]] const std::filesystem::path& path() const noexcept {
        return m_path;
    }

    /// Fill @p dest with the bytes at @p offset.
    /// @throws psrio::FormatError if the range passes the end of the file.
    /// @throws psrio::IoError if a read fails or the file shrinks.
    void read_at(std::uint64_t offset, std::span<std::byte> dest) const {
        if (offset > m_size || dest.size() > m_size - offset) {
            throw FormatError(fmt::format(
                "psrio: {}: read of {} bytes at {} passes the end of the "
                "file ({} bytes)",
                m_path.string(), dest.size(), offset, m_size));
        }
        std::size_t done = 0;
        while (done < dest.size()) {
            const auto want =
                std::min<std::size_t>(dest.size() - done, kMaxRequestBytes);
            const auto at = offset + done;
            if (at >
                static_cast<std::uint64_t>(std::numeric_limits<off_t>::max())) {
                throw IoError(fmt::format("psrio: {}: offset {} is too large",
                                          m_path.string(), at));
            }
            const ssize_t got =
                ::pread(m_fd, dest.data() + done, want, static_cast<off_t>(at));
            if (got < 0) {
                if (errno == EINTR) {
                    continue;
                }
                error_check::throw_errno_io("pread", m_path);
            }
            if (got == 0) {
                throw IoError(fmt::format(
                    "psrio: {}: file ended at byte {} while reading",
                    m_path.string(), at));
            }
            done += static_cast<std::size_t>(got);
        }
    }

    /// Drop or keep cached pages once a range has been consumed.
    void set_drop_cache(bool enable) noexcept { m_drop_cache = enable; }
    [[nodiscard]] bool drop_cache() const noexcept { return m_drop_cache; }

    /// Tell the kernel that [@p offset, @p offset + @p nbytes) will not be
    /// read again. A no-op unless drop-cache is enabled.
    void release([[maybe_unused]] std::uint64_t offset,
                 std::uint64_t nbytes) const noexcept {
        if (!m_drop_cache || m_fd < 0 || nbytes == 0U) {
            return;
        }
#ifdef POSIX_FADV_DONTNEED
        (void)::posix_fadvise(m_fd, static_cast<off_t>(offset),
                              static_cast<off_t>(nbytes), POSIX_FADV_DONTNEED);
#endif
    }

    void reset() noexcept {
        if (m_fd >= 0) {
            (void)::close(m_fd);
        }
        m_fd   = -1;
        m_size = 0;
    }

private:
    void open_fd(bool initial) {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg)
        m_fd = ::open(m_path.c_str(), O_RDONLY | O_CLOEXEC);
        if (m_fd < 0) {
            error_check::throw_errno_io("open", m_path);
        }
        struct stat file_stat{};
        if (::fstat(m_fd, &file_stat) != 0) {
            const int saved = errno;
            close_fd();
            error_check::throw_errno_io(saved, "fstat", m_path);
        }
        if (!S_ISREG(file_stat.st_mode) || file_stat.st_size < 0) {
            close_fd();
            throw IoError(fmt::format("psrio: '{}' is not a regular file",
                                      m_path.string()));
        }
        const auto current_size = static_cast<std::uint64_t>(file_stat.st_size);
        if (initial) {
            m_size = current_size;
        } else if (current_size < m_size) {
            close_fd();
            throw IoError(fmt::format(
                "psrio: '{}' shrank from {} to {} bytes before reopening",
                m_path.string(), m_size, current_size));
        }
#ifdef POSIX_FADV_SEQUENTIAL
        if (m_access == Access::kSequential) {
            (void)::posix_fadvise(m_fd, 0, 0, POSIX_FADV_SEQUENTIAL);
        }
#endif
    }

    /// Linux returns at most about 2 GiB from one read.
    static constexpr std::size_t kMaxRequestBytes = std::size_t{1} << 30U;

    int m_fd{-1};
    std::uint64_t m_size{0};
    bool m_drop_cache{false};
    Access m_access{Access::kSequential};
    std::filesystem::path m_path;
};

} // namespace psrio::detail
