#pragma once

/// GUPPI RAW baseband. `RawReader` walks header and payload bytes.
/// `GuppiReader` presents canonical time samples, including a time sequence
/// of one band. `GuppiReader::open` stitches files that differ in `OBSFREQ`.
#include "psrio/formats/guppi/header.hpp" // IWYU pragma: export
#include "psrio/formats/guppi/reader.hpp" // IWYU pragma: export
