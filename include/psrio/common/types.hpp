#pragma once

#include "psrio/detail/exceptions.hpp"

#include <cstdint>
#include <string_view>

namespace psrio {

/// Sub-byte bit ordering convention.
enum class BitOrder : std::uint8_t {
    kLsbFirst, ///< Least significant bit first (DSPSR convention)
    kMsbFirst, ///< Most significant bit first (PRESTO / sigpyproc convention)
};

/// On-disk sample type for one intensity value.
///
/// `sample_type_from_nbits` follows SIGPROC: `nbits == 32` is IEEE float.
/// `kUInt32` is only for an integer dataset that is actually unsigned 32-bit,
/// such as a Breakthrough Listen HDF5 file. It is not produced from a SIGPROC
/// `nbits` field.
enum class SampleType : std::uint8_t {
    kUInt1,
    kUInt2,
    kUInt4,
    kUInt8,
    kInt8,
    kInt16,
    kUInt16,
    kUInt32,
    kFloat32,
};

/// Bits occupied by one sample of @p type.
[[nodiscard]] constexpr int bits_of(SampleType type) noexcept {
    switch (type) {
    case SampleType::kUInt1:
        return 1;
    case SampleType::kUInt2:
        return 2;
    case SampleType::kUInt4:
        return 4;
    case SampleType::kUInt8:
    case SampleType::kInt8:
        return 8;
    case SampleType::kInt16:
    case SampleType::kUInt16:
        return 16;
    case SampleType::kUInt32:
    case SampleType::kFloat32:
        return 32;
    }
    return 0;
}

/// True if @p type represents IEEE floating point data.
[[nodiscard]] constexpr bool is_floating(SampleType type) noexcept {
    return type == SampleType::kFloat32;
}

/// True if @p type represents an integer quantity (signed or unsigned).
[[nodiscard]] constexpr bool is_integral(SampleType type) noexcept {
    return !is_floating(type);
}

/// True if @p type represents signed numerical values.
[[nodiscard]] constexpr bool is_signed(SampleType type) noexcept {
    return type == SampleType::kInt8 || type == SampleType::kInt16 ||
           type == SampleType::kFloat32;
}

/// Bit depth of @p type as an integer count (negative for floating-point, e.g.
/// -32 for float32).
[[nodiscard]] constexpr int to_nbits(SampleType type) noexcept {
    if (type == SampleType::kFloat32) {
        return -32;
    }
    return bits_of(type);
}

/// Infer the sample type from bit depth and signedness.
[[nodiscard]] constexpr SampleType
sample_type_from_nbits(int nbits, bool is_signed = false) {
    if (is_signed && nbits != 8) {
        throw ValidationError(
            "psrio: signed samples are only supported for 8-bit depth");
    }
    switch (nbits) {
    case 1:
        return SampleType::kUInt1;
    case 2:
        return SampleType::kUInt2;
    case 4:
        return SampleType::kUInt4;
    case 8:
        return is_signed ? SampleType::kInt8 : SampleType::kUInt8;
    case 16:
        return SampleType::kUInt16;
    case 32:
        return SampleType::kFloat32;
    default:
        throw ValidationError("psrio: unsupported sample bit depth");
    }
}

/// Printable string representation of @p type.
[[nodiscard]] constexpr std::string_view
sample_type_name(SampleType type) noexcept {
    switch (type) {
    case SampleType::kUInt1:
        return "uint1";
    case SampleType::kUInt2:
        return "uint2";
    case SampleType::kUInt4:
        return "uint4";
    case SampleType::kUInt8:
        return "uint8";
    case SampleType::kInt8:
        return "int8";
    case SampleType::kInt16:
        return "int16";
    case SampleType::kUInt16:
        return "uint16";
    case SampleType::kUInt32:
        return "uint32";
    case SampleType::kFloat32:
        return "float32";
    }
    return "unknown";
}

} // namespace psrio
