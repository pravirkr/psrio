#pragma once

#include <cstdio>

namespace psrio::detail {

/// RAII owner for a `std::FILE*` opened during POSIX file setup.
class FileHandle {
public:
    explicit FileHandle(std::FILE* file) noexcept : m_file{file} {}

    ~FileHandle() noexcept {
        if (m_file != nullptr) {
            (void)std::fclose(m_file);
        }
    }

    FileHandle(const FileHandle&)            = delete;
    FileHandle& operator=(const FileHandle&) = delete;
    FileHandle(FileHandle&&)                 = delete;
    FileHandle& operator=(FileHandle&&)      = delete;

    [[nodiscard]] explicit operator bool() const noexcept {
        return m_file != nullptr;
    }

    [[nodiscard]] std::FILE* get() const noexcept { return m_file; }

private:
    std::FILE* m_file = nullptr;
};

} // namespace psrio::detail
