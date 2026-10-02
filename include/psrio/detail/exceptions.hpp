#pragma once

#include <cerrno>
#include <filesystem>
#include <format>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>

namespace psrio::error_check {

/// Base error for psrio. Format, validation, and I/O failures derive from this.
class Error : public std::runtime_error {
public:
    explicit Error(const std::string& message) : std::runtime_error(message) {}
};

/// A file could not be opened, stat-ed, or mapped.
class IoError : public Error {
public:
    using Error::Error;
};

/// The bytes are not a readable SIGPROC header or payload.
class FormatError : public Error {
public:
    using Error::Error;
};

/// The header parsed, but its geometry or the caller's request cannot be used.
class ValidationError : public Error {
public:
    using Error::Error;
};

/// Append `file:line` and function name to @p user_msg (developer diagnostics).
[[nodiscard]] inline std::string
located_message(std::string_view user_msg, const std::source_location& loc) {
    return std::format("{}:{}:{}: error: {}\n  in function '{}'",
                       loc.file_name(), loc.line(), loc.column(), user_msg,
                       loc.function_name());
}

/// Internal invariant failure with capture-site metadata for debugging.
///
/// Prefer @ref ValidationError or @ref FormatError for user-facing data errors.
class InvariantError : public Error {
public:
    explicit InvariantError(
        std::string_view user_msg,
        const std::source_location& loc = std::source_location::current())
        : Error(located_message(user_msg, loc)) {}
};

template <typename... Args>
[[noreturn]] inline void throw_io(std::format_string<Args...> fmt,
                                  Args&&... args) {
    throw IoError(std::format(fmt, std::forward<Args>(args)...));
}

template <typename... Args>
[[noreturn]] inline void throw_format(std::format_string<Args...> fmt,
                                      Args&&... args) {
    throw FormatError(std::format(fmt, std::forward<Args>(args)...));
}

template <typename... Args>
[[noreturn]] inline void throw_validation(std::format_string<Args...> fmt,
                                          Args&&... args) {
    throw ValidationError(std::format(fmt, std::forward<Args>(args)...));
}

/// @throws psrio::error_check::IoError with the platform message for @p err.
[[noreturn]] inline void throw_errno_io(int err,
                                        std::string_view operation,
                                        const std::filesystem::path& path) {
    throw IoError(std::format("psrio: {} failed for '{}': {}", operation,
                              path.string(),
                              std::generic_category().message(err)));
}

/// @throws psrio::error_check::IoError when a syscall returns -1 and sets
/// `errno`.
[[noreturn]] inline void throw_errno_io(std::string_view operation,
                                        const std::filesystem::path& path) {
    throw_errno_io(errno, operation, path);
}

} // namespace psrio::error_check

namespace psrio {

using error_check::Error;
using error_check::FormatError;
using error_check::InvariantError;
using error_check::IoError;
using error_check::located_message;
using error_check::throw_errno_io;
using error_check::throw_format;
using error_check::throw_io;
using error_check::throw_validation;
using error_check::ValidationError;

} // namespace psrio
