#pragma once

/// Umbrella header for the psrio public API.
#include <psrio/detail/concepts.hpp>
#include <psrio/detail/endian.hpp>
#include <psrio/detail/mmap.hpp>
#include <psrio/detail/unpack.hpp>

#include <string_view>

namespace psrio {

/// @return Semantic version string matching the CMake project version.
constexpr auto version() noexcept -> std::string_view { return "0.1.0"; }

} // namespace psrio
