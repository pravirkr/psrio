#pragma once

#include <cstddef>
#include <filesystem>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace psrio::detail {

namespace mmap_detail {

inline int open_read_only(const std::filesystem::path& path) {
    return ::open(path.c_str(), O_RDONLY); // NOLINT(cppcoreguidelines-pro-type-vararg)
}

} // namespace mmap_detail

/// Read-only memory mapping of a regular file (POSIX `mmap`).
class MappedFile {
public:
    MappedFile() = default;

    explicit MappedFile(const std::filesystem::path& path) {
        const int file_descriptor = mmap_detail::open_read_only(path);
        if (file_descriptor < 0) {
            throw std::runtime_error(
                "psrio: failed to open file for mapping: " + path.string());
        }

        struct stat stat_buffer{};
        if (::fstat(file_descriptor, &stat_buffer) != 0) {
            ::close(file_descriptor);
            throw std::runtime_error("psrio: failed to stat file: " +
                                     path.string());
        }

        if (stat_buffer.st_size < 0) {
            ::close(file_descriptor);
            throw std::runtime_error("psrio: negative file size: " +
                                     path.string());
        }

        const auto mapped_size = static_cast<std::size_t>(stat_buffer.st_size);
        void* mapped_addr      = MAP_FAILED;
        if (mapped_size > 0) {
            mapped_addr = ::mmap(nullptr, mapped_size, PROT_READ, MAP_PRIVATE,
                                 file_descriptor, 0);
            if (mapped_addr == MAP_FAILED) {
                ::close(file_descriptor);
                throw std::runtime_error("psrio: mmap failed: " +
                                         path.string());
            }
        }

        m_fd   = file_descriptor;
        m_addr = mapped_addr;
        m_size = mapped_size;
    }

    ~MappedFile() { reset(); }

    MappedFile(const MappedFile&)            = delete;
    MappedFile& operator=(const MappedFile&) = delete;

    MappedFile(MappedFile&& other) noexcept { swap(other); }

    MappedFile& operator=(MappedFile&& other) noexcept {
        if (this != &other) {
            reset();
            swap(other);
        }
        return *this;
    }

    [[nodiscard]] auto bytes() const noexcept -> std::span<const std::byte> {
        if (m_size == 0 || m_addr == MAP_FAILED) {
            return {};
        }
        return {static_cast<const std::byte*>(m_addr), m_size};
    }

    [[nodiscard]] auto size() const noexcept -> std::size_t { return m_size; }

    [[nodiscard]] auto empty() const noexcept -> bool { return m_size == 0; }

    void reset() noexcept {
        if (m_addr != MAP_FAILED && m_size > 0) {
            ::munmap(m_addr, m_size);
        }
        if (m_fd >= 0) {
            ::close(m_fd);
        }
        m_addr = MAP_FAILED;
        m_size = 0;
        m_fd   = -1;
    }

private:
    void swap(MappedFile& other) noexcept {
        std::swap(m_fd, other.m_fd);
        std::swap(m_addr, other.m_addr);
        std::swap(m_size, other.m_size);
    }

    int m_fd           = -1;
    void* m_addr       = MAP_FAILED;
    std::size_t m_size = 0;
};

} // namespace psrio::detail
