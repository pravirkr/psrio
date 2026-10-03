#pragma once

/// Umbrella header for the psrio public API.
/// FBH5, PSRFITS, and the PSRDADA ring stay out of this header.
/// DADA files are included.
#include "psrio/astro.hpp"             // IWYU pragma: export
#include "psrio/baseband.hpp"          // IWYU pragma: export
#include "psrio/block_source.hpp"      // IWYU pragma: export
#include "psrio/common/concepts.hpp"   // IWYU pragma: export
#include "psrio/common/types.hpp"      // IWYU pragma: export
#include "psrio/detail/exceptions.hpp" // IWYU pragma: export
#include "psrio/formats/dada.hpp"      // IWYU pragma: export
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
constexpr std::string_view version() noexcept { return "0.3.0"; }

// Main convenient public types
using formats::dada::DadaReader;
using formats::guppi::GuppiReader;
using formats::guppi::GuppiSet;
using formats::sigproc::FilterbankReader;

} // namespace psrio
