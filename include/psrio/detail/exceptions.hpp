#pragma once

#include <stdexcept>
#include <string>

namespace psrio {

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

} // namespace psrio
