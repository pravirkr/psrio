#pragma once

#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/file.hpp"

#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <span>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <utility>

namespace psrio::detail {

/// Read-only memory mapping of a regular file (POSIX `mmap`).
///
/// The descriptor used to establish the mapping is closed after a successful
/// `mmap()`. An existing mapping remains valid when its descriptor is closed.
///
/// Empty files produce an empty byte span and do not call `mmap()`.
///
/// The backing file must not be truncated while mapped; access beyond the new
/// end-of-file may raise `SIGBUS`.
class MappedFile {
public:
    MappedFile() noexcept = default;

    /// @throws psrio::IoError if the path cannot be opened, inspected, or
    /// mapped.
    explicit MappedFile(const std::filesystem::path& path) {
        FileHandle file{std::fopen(path.c_str(), "rb")};
        if (!file) {
            error_check::throw_errno_io("fopen", path);
        }

        const int fd = ::fileno(file.get());
        if (fd < 0) {
            error_check::throw_errno_io("fileno", path);
        }

        struct stat file_stat{};
        if (::fstat(fd, &file_stat) != 0) {
            error_check::throw_errno_io("fstat", path);
        }

        if (!S_ISREG(file_stat.st_mode)) {
            throw IoError(
                std::format("psrio: cannot memory-map non-regular file '{}'",
                            path.string()));
        }

        if (file_stat.st_size < 0) {
            throw IoError(std::format(
                "psrio: file reports a negative size: '{}'", path.string()));
        }

        const auto file_size = static_cast<std::uintmax_t>(file_stat.st_size);
        constexpr auto kMaxMappingSize = static_cast<std::uintmax_t>(
            std::numeric_limits<std::size_t>::max());
        if (file_size > kMaxMappingSize) {
            throw IoError(std::format(
                "psrio: file is too large to map in this process: '{}'",
                path.string()));
        }

        const auto mapping_size = static_cast<std::size_t>(file_size);
        if (mapping_size == 0) {
            return;
        }

        void* const address =
            ::mmap(nullptr, mapping_size, PROT_READ, MAP_PRIVATE, fd, 0);
        if (address == MAP_FAILED) {
            error_check::throw_errno_io("mmap", path);
        }

        m_addr = address;
        m_size = mapping_size;
    }

    ~MappedFile() noexcept { reset(); }

    MappedFile(const MappedFile&)            = delete;
    MappedFile& operator=(const MappedFile&) = delete;

    MappedFile(MappedFile&& other) noexcept
        : m_addr{std::exchange(other.m_addr, MAP_FAILED)},
          m_size{std::exchange(other.m_size, 0)} {}

    MappedFile& operator=(MappedFile&& other) noexcept {
        if (this != &other) {
            reset();
            m_addr = std::exchange(other.m_addr, MAP_FAILED);
            m_size = std::exchange(other.m_size, 0);
        }
        return *this;
    }

    [[nodiscard]] std::span<const std::byte> bytes() const noexcept {
        if (m_addr == MAP_FAILED || m_size == 0) {
            return {};
        }
        return {static_cast<const std::byte*>(m_addr), m_size};
    }

    [[nodiscard]] std::size_t size() const noexcept { return m_size; }

    [[nodiscard]] bool empty() const noexcept { return m_size == 0; }

    [[nodiscard]] bool is_mapped() const noexcept {
        return m_addr != MAP_FAILED;
    }

    void advise_sequential() const noexcept {
        if (m_addr == MAP_FAILED || m_size == 0) {
            return;
        }
#ifdef POSIX_MADV_SEQUENTIAL
        (void)::posix_madvise(m_addr, m_size, POSIX_MADV_SEQUENTIAL);
#elif defined(MADV_SEQUENTIAL)
        (void)::madvise(m_addr, m_size, MADV_SEQUENTIAL);
#endif
    }

    void reset() noexcept {
        if (m_addr != MAP_FAILED) {
            (void)::munmap(m_addr, m_size);
        }
        m_addr = MAP_FAILED;
        m_size = 0;
    }

private:
    void* m_addr       = MAP_FAILED;
    std::size_t m_size = 0;
};

} // namespace psrio::detail
