#pragma once

/**
 * @file sigproc.hpp
 * @brief Main entry point for the SIGPROC filterbank (`.fil`) reader and
 * writer.
 *
 * This file includes the header, reader, and writer for the SIGPROC
 * filterbank format.
 *
 * @see psrio::formats::sigproc::FilterbankHeader
 * @see psrio::formats::sigproc::FilterbankReader
 * @see psrio::formats::sigproc::FilterbankWriter
 */
#include "psrio/formats/sigproc/header.hpp" // IWYU pragma: export
#include "psrio/formats/sigproc/reader.hpp" // IWYU pragma: export
#include "psrio/formats/sigproc/writer.hpp" // IWYU pragma: export
