#pragma once

/// Umbrella header for the psrio public API.
/// FBH5, PSRFITS, and PSRDADA stay out of this header.
#include "psrio/astro.hpp"             // IWYU pragma: export
#include "psrio/block_source.hpp"      // IWYU pragma: export
#include "psrio/common/concepts.hpp"   // IWYU pragma: export
#include "psrio/common/types.hpp"      // IWYU pragma: export
#include "psrio/detail/exceptions.hpp" // IWYU pragma: export
#include "psrio/formats/guppi.hpp"     // IWYU pragma: export
#include "psrio/formats/memory.hpp"    // IWYU pragma: export
#include "psrio/formats/presto.hpp"    // IWYU pragma: export
#include "psrio/formats/sigproc.hpp"   // IWYU pragma: export
#include "psrio/header.hpp"            // IWYU pragma: export
#include "psrio/packed.hpp"            // IWYU pragma: export
#include "psrio/timeseries.hpp"        // IWYU pragma: export

#include <string_view>

namespace psrio {

/// @return Semantic version string matching the CMake project version.
constexpr std::string_view version() noexcept { return "0.2.0"; }

// Main convenient public types
using formats::sigproc::FilterbankReader;

} // namespace psrio
