#pragma once

/**
 * @file presto.hpp
 * @brief Main entry point for the PRESTO time series (`.dat`, `.inf`) reader
 * and writer.
 *
 * @see psrio::formats::presto::TimeSeriesReader
 * @see psrio::formats::presto::PrestoTimeSeriesWriter
 * @see psrio::formats::presto::parse_inf
 */
#include "psrio/formats/presto/inf.hpp"    // IWYU pragma: export
#include "psrio/formats/presto/reader.hpp" // IWYU pragma: export
#include "psrio/formats/presto/writer.hpp" // IWYU pragma: export
