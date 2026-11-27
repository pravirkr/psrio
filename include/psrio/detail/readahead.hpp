#pragma once

#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/rfile.hpp"

#include <algorithm>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <deque>
#include <exception>
#include <mutex>
#include <span>
#include <thread>
#include <utility>
#include <vector>

namespace psrio::detail {

/// Streams a byte range of a file on a background thread, ahead of the
/// caller.
///
/// One thread reads [start, end) back to back in blocks of kBlockBytes that
/// sit on kBlockBytes boundaries of the file, keeping up to `depth` blocks
/// ahead of the caller. Parallel file systems serve such aligned requests
/// fastest, and the file sees one continuous stream while the caller works.
/// read() copies any forward-moving sequence of ranges out of the blocks;
/// bytes the caller skips are read but not used. Destroying the object stops
/// the thread after the block in progress.
class ReadAhead {
public:
    static constexpr std::uint64_t kBlockBytes = std::uint64_t{16} << 20U;

    /// @param file A descriptor this object owns (see ReadFile::duplicate()).
    /// @param start First byte the caller will read.
    /// @param end End of the bytes to stream (clamped to the file size).
    /// @param depth Blocks kept read ahead of the caller.
    /// @param block_bytes Block size and alignment.
    ReadAhead(ReadFile file,
              std::uint64_t start,
              std::uint64_t end,
              std::size_t depth,
              std::uint64_t block_bytes = kBlockBytes)
        : m_file(std::move(file)),
          m_end(std::min(end, m_file.size())),
          m_depth(std::max<std::size_t>(depth, 1U)),
          m_block(std::max<std::uint64_t>(block_bytes, 1U)),
          m_next((start / m_block) * m_block),
          m_position(start) {
        m_file.advise_sequential();
        m_thread = std::thread([this] { work(); });
    }

    ReadAhead(const ReadAhead&)            = delete;
    ReadAhead& operator=(const ReadAhead&) = delete;
    ReadAhead(ReadAhead&&)                 = delete;
    ReadAhead& operator=(ReadAhead&&)      = delete;

    ~ReadAhead() {
        {
            const std::scoped_lock lock(m_mutex);
            m_stop = true;
        }
        m_wake.notify_all();
        if (m_thread.joinable()) {
            m_thread.join();
        }
    }

    /// True when read(@p offset, ...) can be served: @p offset is not before
    /// the bytes already handed out and lies inside the streamed range.
    [[nodiscard]] bool covers(std::uint64_t offset,
                              std::uint64_t nbytes) const noexcept {
        return offset >= m_position && offset <= m_end &&
               nbytes <= m_end - offset;
    }

    /// Copy the bytes at @p offset into @p dest, waiting for them to be read.
    /// Offsets must not move backward between calls.
    /// @throws ValidationError if the range is not covered.
    /// @throws The error of a failed read.
    void read(std::uint64_t offset, std::span<std::byte> dest) {
        if (!covers(offset, dest.size())) {
            throw ValidationError("psrio: read-ahead cannot serve this range");
        }
        std::size_t done = 0;
        std::unique_lock lock(m_mutex);
        while (done < dest.size()) {
            const auto at = offset + done;
            m_ready_cv.wait(lock,
                            [this] { return !m_ready.empty() || m_error; });
            if (m_ready.empty()) {
                std::rethrow_exception(m_error);
            }
            auto& block     = m_ready.front();
            const auto last = block.offset + block.bytes.size();
            if (last <= at) {
                // Every byte of this block lies before the caller.
                m_spare.push_back(std::move(block.bytes));
                m_ready.pop_front();
                m_wake.notify_all();
                continue;
            }
            const auto within = static_cast<std::size_t>(at - block.offset);
            const auto take =
                std::min(dest.size() - done, block.bytes.size() - within);
            std::memcpy(dest.data() + done, block.bytes.data() + within, take);
            done += take;
        }
        m_position = offset + dest.size();
    }

private:
    struct Block {
        std::uint64_t offset{0};
        std::vector<std::byte> bytes;
    };

    void work() {
        while (true) {
            std::vector<std::byte> buffer;
            std::uint64_t at = 0;
            {
                std::unique_lock lock(m_mutex);
                m_wake.wait(lock, [this] {
                    return m_stop || m_ready.size() < m_depth;
                });
                if (m_stop || m_next >= m_end) {
                    return;
                }
                if (!m_spare.empty()) {
                    buffer = std::move(m_spare.back());
                    m_spare.pop_back();
                }
                at = m_next;
            }
            const auto stop = std::min(m_end, ((at / m_block) + 1U) * m_block);
            try {
                buffer.resize(static_cast<std::size_t>(stop - at));
                m_file.read_at(at, buffer);
                m_file.release(at, stop - at);
            } catch (...) {
                const std::scoped_lock lock(m_mutex);
                m_error = std::current_exception();
                m_ready_cv.notify_all();
                return;
            }
            {
                const std::scoped_lock lock(m_mutex);
                m_ready.push_back(
                    Block{.offset = at, .bytes = std::move(buffer)});
                m_next = stop;
            }
            m_ready_cv.notify_all();
        }
    }

    ReadFile m_file;
    std::uint64_t m_end;
    std::size_t m_depth;
    std::uint64_t m_block;
    std::uint64_t m_next;
    std::uint64_t m_position;
    std::mutex m_mutex;
    std::condition_variable m_wake;
    std::condition_variable m_ready_cv;
    std::deque<Block> m_ready;
    std::vector<std::vector<std::byte>> m_spare;
    std::exception_ptr m_error;
    bool m_stop{false};
    // Last: started after, and joined before, everything it uses.
    std::thread m_thread;
};

} // namespace psrio::detail
