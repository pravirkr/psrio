# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### sigpyproc `spp_digifil` parity

- **Behaviour change**: `psrio digifil` on search-mode PSRFITS no longer
  applies `DAT_SCL`/`DAT_OFFS`/`ZERO_OFF` or `DAT_WTS` unless `--scloffs` is
  given, as dspsr and `spp_digifil`. Pass `--scloffs` for the previous
  output. Library readers are unchanged: `SearchFile`/`SearchStream` still
  decode calibrated, weighted Stokes I by default.
- `psrfits::DecodeOptions` (`scloffs`, `weights`, `product`) for
  `SearchFile` and `SearchStream`. `product` reads one stored polarisation
  (AA, BB, CR, CI) instead of Stokes I. `digifil::decode_options(config)`
  gives dspsr's choice.
- `psrio digifil`:
  - `-P` on PSRFITS selects a stored product.
  - New `--scloffs`, `-U,--max-ram` (`Config::max_ram_mb`), `-B` (ignored)
    and `-d 1`.
  - `-o` is optional and defaults to the input stem with `.fil`.
- On a Parkes UWL file psrio matches `spp_digifil` to within 1e-5 of output
  bytes, all at ±1 level, for the default, `--scloffs`, `-P`, `-U` and
  `-I 0` cases. On SIGPROC input 1-, 2- and float32-bit output is identical;
  4-, 8- and 16-bit output differs only at .5 rounding boundaries, because
  sigpyproc's `fastmath` kernel fuses the multiply-add.

### Earlier unreleased work

Fixes and performance work from validation on Parkes UWL search-mode data
(3 × 12 GB `.sf`, 6656 channels, `AABBCRCI`, on BeeGFS).

- **I/O: positioned reads instead of whole-file mappings**:
  - `SearchFile` and `FilterbankReader` read with `pread` (`detail::ReadFile`).
    On a parallel file system page faults on a mapping ran at 150–180 MB/s;
    large aligned reads run at 1.1–1.4 GB/s. The resident set no longer grows
    with the file: a full 3-file `digifil` peaked at 15.5 GB and now stays
    near 0.6 GB.
  - Once reads run forward, `detail::ReadAhead` streams the rest of the file
    on a background thread in 16 MiB blocks on 16 MiB file boundaries
    (alignment alone is worth ~25% on BeeGFS), four blocks ahead of the
    caller. Seeks fall back to direct reads. `set_prefetch(false)` turns it
    off; `release_buffers()` frees it, and the streams call it when they move
    past a file.
  - PSRFITS metadata (`TSUBINT`, `OFFS_SUB`, `DAT_FREQ`, `DAT_WTS`,
    `DAT_OFFS`, `DAT_SCL`) is read once per row in one request. Headers are
    parsed from file prefixes (`FilterbankHeader::read`,
    `parse_prefix`; `psrfits::detail::Source`, `locate_cell`, `read_cell`).
  - `FilterbankReader::view_bytes()`/`view_block()` keep their zero-copy
    contract by mapping the file on first use.
  - Opt-in page-cache eviction: `set_drop_cache(true)` or `PSRIO_DROP_CACHE=1`.
  - `psrio::psrio` now links `Threads::Threads`.
- **CPU**: hot loops (PSRFITS decode, digifil rescale/apply/encode) use
  `PSRIO_OMP(...)` pragmas, active only when the consumer compiles with
  OpenMP; psrio gains no OpenMP dependency. Every parallel loop keeps each
  output's arithmetic in serial order, so results do not depend on the thread
  count. Measured on BeeGFS with 8 cores (cold cache; run-to-run spread is
  large): full 3-file `psrio digifil -b 8` 258 s → 27–39 s; reading a 9.3 GB
  8-bit `.fil` 55 s → 9 s. Outputs are byte-identical to 0.5.0 at 1, 3 and 8
  threads.
- **`digifil::Source`**: digifil as an in-memory `BlockReader` over any
  `BlockSource` (including a multi-file PSRFITS stream), so a pipeline can
  read 8-bit redigitised data without writing a `.fil`. It reads the next
  input block on another thread while the pipeline runs. `digifil::run` now
  drives it. Also `digifil::header_of(BlockSource)` and `default_gulp`.
- **Robust rescaling (opt-in)**: `Config::rescale_clip` / `psrio digifil
  --rescale-clip k` leaves values beyond k σ of a channel's statistics out of
  its next statistics. The default stays dspsr's plain mean and variance.
- **`BlockSource::refdm()` and `BlockSource::sigproc_header()`**, detected
  from the wrapped reader (the `BlockReader` concept is unchanged).
  `refdm()` on `FilterbankReader`, `SearchFile`, and `SearchStream`.
  `SearchFile::frontend()`, `backend()`, `obs_mode()`; `SearchStream::nfiles()`
  and `file(i)`.
- **Fixes**:
  - `ftop()`/`fbottom()` (SIGPROC and `psrio::Header`) returned swapped band
    edges when `foff > 0`; a frequency table now gives edges too.
  - `psrio dedisperse --dm 0,411 in.fil` no longer takes `in.fil` as a DM;
    `--dm` may be repeated.
  - An installed package built with `PSRIO_WITH_HDF5=ON` and a fetched
    HighFive failed `find_package(psrio)` on `HighFive::HighFive`; the
    HighFive headers are now installed and `psrio::fbh5` links `HDF5::HDF5`.
- **CLI**:
  - `inspect` prints the ISO 8601 UTC start, and for PSRFITS the frontend,
    backend, observing mode, polarisation type and count, `NSBLK`, subints
    per file, and `NSUBOFFS`.
  - `-j,--threads` on `slice`, `splice`, `decimate`, `convert`, `digifil`. The
    defaults follow `OMP_NUM_THREADS`/`OMP_THREAD_LIMIT` (HTCondor sets them);
    `dedisperse` no longer defaults to every core of the host.
  - With OpenMP and no wait policy chosen, `psrio` re-executes itself with
    `OMP_WAIT_POLICY=passive`, so idle threads do not spin through file
    waits (this cut charged CPU time about threefold).
  - `PSRIO_CLI_WITH_DMT=OFF` builds the CLI without dmt/FFTW (no
    `dedisperse`). `psrio` is installed with `INSTALL_RPATH_USE_LINK_PATH`.

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
