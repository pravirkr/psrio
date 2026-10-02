# Changelog

All notable changes to this project will be documented in this file.

## [0.2.0] - 2026-10-02

- **Intensity block front**:
  - `psrio::concepts::BlockReader` and move-only `psrio::BlockSource` for filterbank-like streams.
  - `read_block` copies packed time samples and may return a short count. `read_bytes` is unchanged and still requires an exact byte count.
  - `skip` moves the sample cursor forward or backward. `read_samples` unpacks into caller storage.
  - Metadata getters: `nchans`, `nifs`, `nbits`, `sample_type`, `bytes_per_sample`, `tsamp`, `tstart`, `fch1`, `foff`, `beam`, `spectra_rate`, `utc_start`, and `has_nsamples`.
  - `MemoryBlock` for an in-memory packed stream. `reverse_channels`, and `astro::mjd_to_time`.

- **GUPPI RAW baseband front** (`psrio::formats::guppi::RawReader`):
  - Parses 80-byte headers, DIRECTIO padding, and payload bytes. It does not model `BlockReader`.

- **Optional readers** (off unless the matching CMake option is on; neither is included by `psrio.hpp`):
  - FBH5 via `PSRIO_WITH_HDF5` and `psrio::fbh5`.
  - PSRDADA via `PSRIO_WITH_PSRDADA` and `psrio::psrdada`. `RingReader` streams one intensity observation from a ring.

## [0.1.0] - 2026-10-01

### Added

- **Core Architecture**:
  - Pure header-only C++20 library architecture with zero compiled binary artifacts.
  - Public umbrella header `<psrio/psrio.hpp>` exposing core types and readers.
  - Generic single-channel streaming reader concept `psrio::concepts::TimeSeriesReader`.
  - Memory-mapped I/O (`mmap`) support with zero-copy views and caller-managed buffer unpacking.

- **SIGPROC Filterbank (`.fil`)**:
  - Binary header parser (`formats::sigproc::FilterbankHeader`) supporting all standard keywords, frequency tables, and vendor extra keys.
  - Chunked streaming reader (`formats::sigproc::FilterbankReader`) supporting 1, 2, 4, 8 (signed/unsigned), 16, and 32-bit samples.
  - Sub-byte bit ordering options: `BitOrder::kLsbFirst` (DSPSR default) and `BitOrder::kMsbFirst` (PRESTO convention).

- **Unified Observational Metadata (`psrio::Header`)**:
  - Flat, idiomatic C++20 struct normalizing physical, timing, celestial, and telescope parameters.
  - Astronomical conversion helpers (MJD to Gregorian date, RA/Dec conversions, duration formatting, frequency center/bandwidth).
  - Factory methods: `Header::from_sigproc`, `Header::from_inffile`, and auto-detecting `Header::from_file`.
  - PRESTO `.inf` exporter: `Header::make_inf`.

- **Single-Channel Time Series (`psrio::timeseries`)**:
  - `psrio::TimeSeriesReader`: Lite streaming reader that parses metadata on construction and unpacks directly into user-provided `std::span<float>` with zero heap allocations.
  - `psrio::TimeSeries`: In-memory 1D container with data span access, metadata inspection, and disk exporters (`to_tim`, `to_dat`).
  - Robust PRESTO `.inf` metadata parser supporting multiple PRESTO tool output dialects.

- **Testing & Quality Assurance**:
  - Comprehensive Catch2 v3 test suite.
  - Real observational test verification using `sigpyproc3` datasets.
  - Unit tests for endian swapping, bit packing/unpacking, sexagesimal coordinate conversion.

- **CMake Packaging & Integration**:
  - Support for `CPM.cmake`, `FetchContent`, and `find_package(psrio CONFIG)`.
  - Automatic `cxx_std_20` requirement propagation to consumer targets.
  - Isolation of Catch2 and tests when included as a subproject (`PSRIO_BUILD_TESTS=OFF`).
