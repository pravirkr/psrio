#pragma once

/// Umbrella header for the psrio public API.
#include "psrio/astro.hpp"
#include "psrio/detail/concepts.hpp" // IWYU pragma: export
#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/mmap.hpp"
#include "psrio/detail/unpack.hpp"
#include "psrio/formats/sigproc.hpp"

#include <string_view>

namespace psrio {

/// @return Semantic version string matching the CMake project version.
constexpr std::string_view version() noexcept { return "0.1.0"; }

} // namespace psrio
