# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

- **Behaviour changes**:
  - `psrio digifil` on PSRFITS applies `DAT_SCL`/`DAT_OFFS`/`ZERO_OFF` only with `--scloffs` and never `DAT_WTS`, as dspsr and `spp_digifil`. Library readers still decode calibrated, weighted Stokes I.
  - `ftop()`/`fbottom()` (SIGPROC and `psrio::Header`) return the band top and bottom whichever way `foff` points; they were swapped for `foff > 0`.
- **I/O**:
  - `FilterbankReader` and `SearchFile` read with `pread` instead of `mmap`, with a 16 MiB-aligned read-ahead thread once reads run forward. A full 3-file Parkes UWL `digifil` goes from 258 s to 31 s and 15.5 GB to 0.6 GB peak RSS, byte-identical to 0.5.0.
  - Internal buffers are bounded (16 MiB staging, 64 MiB read-ahead) and freed when a stream moves past a file. `set_prefetch()`, `release_buffers()`, and `PSRIO_DROP_CACHE=1` / `set_drop_cache()`.
  - `psrio::psrio` links `Threads::Threads`.
- **Library**:
  - `psrfits::DecodeOptions` (`scloffs`, `weights`, `product`) for `SearchFile` and `SearchStream`.
  - `digifil::Source` (digifil as a `BlockReader`, no intermediate file) and opt-in `Config::rescale_clip`.
  - `stats::ChannelStats`: one-pass per-channel mean, variance, skewness, kurtosis, min and max, independent of the thread count.
  - `psrio/clean.hpp`: double-MAD and IQRM channel masks, zero-DM, sample zapping, and `clean::Source`.
  - `BlockSource::refdm()`, `BlockSource::sigproc_header()`, and `astro::mjd_to_iso8601`.
  - Hot loops use `PSRIO_OMP(...)`, active only when the consumer builds with OpenMP, and vectorise; results do not depend on the thread count.
- **CLI**:
  - `stats`, `bandpass`, and `clean rfi|zap-channels|zap-times|zerodm`, matching `spp_stats`, `spp_bandpass`, and `spp_clean`. Every `clean` subcommand writes a new file.
  - `header --set key=value [--dry-run]` edits a `.fil` header in place, as `spp_header update`: all or nothing, same header size, values checked as sigpyproc.
  - `digifil`: `--scloffs`, `-P` stored product, `-U`, `-B`, `-d 1`, `--rescale-clip`, optional `-o`. Matches `spp_digifil` byte for byte in a native build.
  - `dedisperse --device cpu|cuda[:N]|hip[:N]` through dmt (`PSRIO_DMT_GPU`; untested on a GPU). `PSRIO_CLI_WITH_DMT=OFF` builds without dmt.
  - `-j,--threads` follows `OMP_NUM_THREADS`; `psrio` runs OpenMP threads with a passive wait policy.
  - `inspect` shows the ISO UTC start and PSRFITS frontend, backend, mode, `NSBLK`, subints per file, and `NSUBOFFS`.
  - `PSRIO_CLI_CPU_ARCH` (default `native`, CI `x86-64-v3`) sets the CLI's `-march`.
- **Fixes**:
  - `dedisperse --dm 0,411 in.fil` no longer reads `in.fil` as a DM.
  - Installed `PSRIO_WITH_HDF5=ON` packages built with a fetched HighFive now pass `find_package(psrio)`.

## [0.5.0] - 2026-11-11

- **Search-mode PSRFITS reader (`psrio::formats::psrfits`)**:
  - Header-only `SearchFile` reads uncompressed `SUBINT` binary tables without CFITSIO and streams calibrated float32 Stokes I (`BlockReader`), matching `sigpyproc` / `dspsr` scale, offset, zero-offset, weight, and 2-bit level rules.
  - Supports 1, 2, 4, 8 (signed/unsigned), and 16-bit samples, arbitrary `TDIM` axis permutations, and `1PB`/`1QB`/`1PE`/`1PD` FITS heap columns (`THEAP`).
- **Contiguous multi-file streams**:
  - `FilterbankStream` (SIGPROC `.fil`) and `SearchStream` (search-mode PSRFITS `.sf`/`.fits`/`.sfits`) join time-ordered split observations into a single `BlockReader` after verifying header identity, `sigpyproc` sample-gap tolerance (`1e-8 + 0.1 * tsamp`), and `NSUBOFFS` progression.
  - `psrio::open_block_source` and `psrio::formats::fbh5::open_block_source` open single files or contiguous multi-file spans as a type-erased `BlockSource`.
- **Filterbank and time-series writers**:
  - `FilterbankHeader::to_bytes()` and `write_to()` serialize SIGPROC headers while preserving unknown vendor keys in order and omitting legacy `nsamples`.
  - `FilterbankWriter` streams packed bytes (`write_block`, `write_bytes`) or unpacked samples (`write_samples`) across 1, 2, 4, 8 (signed/unsigned), 16, and 32-bit depths.
  - `SigprocTimeSeriesWriter` (`.tim`), `PrestoTimeSeriesWriter` (`.dat`/`.inf`), and format-dispatching `TimeSeriesWriter` stream single-channel 32-bit float time series and write the final sample count on `close()`.
- **Processing & redigitisation (`scrunch.hpp`, `digifil.hpp`)**:
  - `fscrunch_inplace` and gulp-invariant `TimeScruncher`.
  - `psrio::digifil::Pipeline` and `digifil::run`: port of `dspsr`'s `digifil` / `sigpyproc`'s `spp_digifil` (IF selection, `fscrunch`, `tscrunch`, running or constant rescaling, and 1/2/4/8/16/-32-bit digitisation).
- **CLI application (`PSRIO_BUILD_APPLICATIONS`)**:
  - Unified `psrio` executable (built only when `PSRIO_BUILD_APPLICATIONS=ON`; `psrio::psrio` remains strictly header-only) with subcommands `inspect`, `slice`, `splice`, `decimate`, `convert`, `digifil`, and `dedisperse` (powered by `dmt` `ddmt`, `sdmt`, and `fdmt`).

## [0.4.0] - 2026-11-03

- **Lower compiler floor**: GCC 11+ and Clang 15+ (was GCC 13+ / Clang 18+).
  GCC 11 is the oldest libstdc++ with `std::source_location`, `std::bit_cast`,
  `std::numbers`, and a usable `<ranges>`. Clang 15 is the floor from before
  `std::format` required GCC 13 and Clang 18.
  - Diagnostic and metadata strings use `fmt::format`. `psrio::psrio` links
    `fmt::fmt-header-only`. A build fetches fmt 12.1.0 (hash-pinned) unless that
    target already exists. psrio does not install fmt.
  - An installed package calls `find_dependency(fmt 12.1 CONFIG)`. fmt 12.1 or
    newer must already be on the prefix.
  - A numeric HDF5 `VERSION` is formatted with `{:.6g}`.
  - New option `PSRIO_INSTALL` (default: top-level build only) guards the install rules.
  - CI builds GCC 11 and Clang 15 on Ubuntu 22.04, GCC 13 and Clang 18 on
    Ubuntu 24.04, and a consumer of the installed package.

## [0.3.0] - 2026-10-03

- **Baseband block front**:
  - `psrio::concepts::BasebandReader`, move-only `psrio::BasebandSource`, and `psrio::FrequencyStitch`.
  - Canonical reads are time-major: `[time][antenna][channel][polarization][component]`.
  - `GuppiReader` indexes GUPPI RAW blocks, drops `OVERLAP`, transposes `1SFA` payloads, and unpacks 2, 4,
  8, and 16-bit complex samples. `GuppiReader::open` joins a time sequence and stitches files that differ
  in `OBSFREQ`.
  - `DadaReader` reads header-only DADA files (`NDIM` 1 or 2, `NBIT` 1/2/4/8/16/32/-32) and joins a time
  sequence by `OBS_OFFSET`.
  - `SampleType::kInt16` for signed 16-bit voltages. SIGPROC 16-bit intensity stays `kUInt16`.
  - `read_voltages` exports a dual-polarisation block as `TFPRI` or `FTPRI`, with channel 0 at the lowest
  frequency when asked. `read_block` stays time-major in native channel order.
  - `psrio::BitsInfo` (and `psrio::DigitizationInfo`) `constexpr` class providing standard SIGPROC/PRESTO bit packing factors, byte item sizes, and digitization statistics (min, max, mean, sigma, scale) across 1, 2, 4, 8, 16, and 32 bits.
  - Public `psrio::pack_sub_byte` with both caller-managed span and convenience allocating overloads.
  - In-place packing (`psrio::pack_inplace`) and backwards in-place unpacking (`psrio::unpack_inplace`) with zero memory allocations.

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
