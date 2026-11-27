# psrio

A high-performance, header-only C++20 library for reading, parsing, and streaming time-domain radio astronomy data (pulsar baseband, filterbank, and single-channel time series).

`psrio` provides zero-overhead, type-safe ingestion and writing of observational data formats (SIGPROC `.fil` / `.tim`, search-mode PSRFITS `.sf` / `.fits`, PRESTO `.dat` / `.inf`, GUPPI RAW, and DADA files), plus an optional `psrio` CLI (`inspect`, `slice`, `splice`, `decimate`, `convert`, `digifil`, `dedisperse`). Intensity formats share one block-streaming front. Baseband formats share another.

## Key Features

- **Pure Header-Only C++20**: Zero compiled library binaries. Simply link `psrio::psrio` and `#include <psrio/psrio.hpp>`.
- **Fast streaming on parallel file systems**: SIGPROC and search-mode PSRFITS readers read with `pread` and, once reads run forward, stream the rest of the file on a background thread in large aligned blocks (BeeGFS: ~1.1–1.4 GB/s, against ~0.17 GB/s through a memory mapping). Memory stays flat whatever the file size. Readers unpack directly into caller-managed buffers (`std::span<float>`, pinned CUDA host memory, etc.); `view_bytes()` still offers a zero-copy mapped view. Hot loops use OpenMP when the consumer builds with it (`-fopenmp`), and give the same results on any thread count.
- **Unified Metadata Representation**: A single flat, idiomatic [`psrio::Header`](include/psrio/header.hpp) struct.
- **One Header-Only Dependency**: Base formats require the ISO C++ standard library, POSIX system APIs (`pread`, `<sys/mman.h>`), threads, and [fmt](https://github.com/fmtlib/fmt) 12.1 or newer. A build fetches fmt 12.1.0 with CPM unless `fmt::fmt-header-only` already exists. An installed psrio does not ship fmt. FBH5 and PSRDADA are optional and stay out of `psrio.hpp` unless `PSRIO_WITH_HDF5` or `PSRIO_WITH_PSRDADA` is turned on.
- **Block streaming**: `psrio::concepts::BlockReader` and the move-only `psrio::BlockSource` read packed bytes (`read_block`) or unpacked floats (`read_samples`), and move the sample cursor with `seek`, `rewind`, and `skip`.

---

## Requirements

- **C++ Standard**: C++20 (`-std=c++20`)
- **Compilers**: GCC 11+, Clang 15+
- **Build System**: CMake 3.21+
- **Dependencies**: fmt 12.1 or newer (header-only). Fetched at configure time when `fmt::fmt-header-only` is missing. An installed psrio finds it with `find_package(fmt 12.1)`.

---

## Adding `psrio` to Your CMake Project

`psrio` exports the interface target `psrio::psrio`. Linking to this target automatically configures include directories and enables C++20 compilation features on your target.

### Method 1: Using `CPM.cmake` (Recommended)

Add `psrio` directly from GitHub using [CPM.cmake](https://github.com/cpm-cmake/CPM.cmake):

```cmake
CPMAddPackage(
  NAME psrio
  GITHUB_REPOSITORY pravirkr/psrio
  GIT_TAG v0.6.0  # Or a specific branch/commit
)

target_link_libraries(my_downstream_target PRIVATE psrio::psrio)
```

Define `fmt::fmt-header-only` before adding psrio when the parent already builds fmt 12.1 or newer. Otherwise psrio fetches its own copy.

### Method 2: Using CMake `FetchContent`

Using built-in CMake 3.21+ `FetchContent`:

```cmake
include(FetchContent)

FetchContent_Declare(
  psrio
  GIT_REPOSITORY https://github.com/pravirkr/psrio.git
  GIT_TAG v0.6.0  # Or a specific branch/commit
)

FetchContent_MakeAvailable(psrio)

target_link_libraries(my_downstream_target PRIVATE psrio::psrio)
```

### Method 3: Using `find_package` (Pre-installed)

Install fmt 12.1 or newer first, then install `psrio` into the same prefix or any prefix CMake already searches. psrio does not install fmt.

```bash
cmake -B build -DPSRIO_BUILD_TESTS=OFF -DPSRIO_BUILD_APPLICATIONS=OFF -DCMAKE_INSTALL_PREFIX=/usr/local
sudo cmake --install build
```

Then in your project's `CMakeLists.txt`:

```cmake
find_package(psrio CONFIG REQUIRED)

target_link_libraries(my_downstream_target PRIVATE psrio::psrio)
```


> [!NOTE]
> `PSRIO_BUILD_TESTS` and `PSRIO_BUILD_APPLICATIONS` both default to `OFF` automatically whenever `psrio` is consumed as a dependency via `CPM`, `FetchContent`, or `add_subdirectory`.

---

## Usage Examples

### 1. Zero-Allocation Streaming with `TimeSeriesReader`

Ideal for real-time processing and GPU pipelines where you want to unpack data directly into caller-owned buffers (e.g. pinned host memory) without heap allocations:

```cpp
#include <psrio/timeseries.hpp> // Or umbrella <psrio/psrio.hpp>

#include <iostream>
#include <span>
#include <vector>

int main() {
    // 1. Open file (.tim or .dat/.inf) - only parses metadata, 0 data array allocation
    psrio::TimeSeriesReader reader("observation.tim"); // or "observation.dat"

    std::cout << "Source:   " << reader.header().source << "\n"
              << "Samples:  " << reader.nsamples() << "\n"
              << "Sampling: " << reader.dt() * 1e6 << " us\n"
              << "DM:       " << reader.header().dm << "\n";

    // 2. Stream directly into a caller-managed buffer
    std::vector<float> buffer(reader.nsamples());
    reader.read_data(buffer);

    // 3. Or chunked reading / streaming
    reader.rewind();
    std::vector<float> chunk(1024);
    while (reader.tell() < reader.nsamples()) {
        std::uint64_t n_read = reader.read_samples(1024, chunk);
        // Process chunk...
    }

    return 0;
}
```

### 2. In-Memory TimeSeries Processing (`sigpyproc3` Style)

Loads full time series data into a contiguous `std::vector<float>` alongside the observational `Header`:

```cpp
#include <psrio/timeseries.hpp>

#include <iostream>

int main() {
    // Load from SIGPROC .tim or PRESTO .dat (auto-detects paired .inf)
    auto ts = psrio::TimeSeries::from_file("observation.dat");

    std::cout << "Observation: " << ts.header().source << "\n"
              << "Telescope:   " << ts.header().telescope << "\n"
              << "Duration:    " << ts.tobs() << " s\n";

    // Direct access to float samples
    std::span<float> samples = ts.data();
    for (float& val : samples) {
        val -= 100.0f; // In-place processing
    }

    // Export modified series back to disk
    ts.to_tim("processed.tim");
    ts.to_dat("processed"); // Generates processed.dat and processed.inf

    return 0;
}
```

### 3. Streaming Multi-Channel Filterbank Data (`.fil`)

Stream spectrometer/filterbank channels chunk by chunk:

```cpp
#include <psrio/formats/sigproc.hpp>

#include <iostream>
#include <vector>

int main() {
    psrio::formats::sigproc::FilterbankReader reader("observation.fil");

    const auto& hdr = reader.header();
    std::cout << "Channels: " << hdr.nchans << "\n"
              << "Bits:     " << hdr.nbits << "\n"
              << "F_top:    " << hdr.fch1.value_or(0.0) << " MHz\n";

    // Read 256 time spectra at a time
    const std::uint64_t chunk_samples = 256;
    std::vector<float> spectra(chunk_samples * hdr.nchans);

    while (reader.tell() < reader.nsamples()) {
        std::uint64_t read = reader.read_samples(chunk_samples, spectra);
        // spectra contains read * hdr.nchans unpacked floats
    }

    return 0;
}
```

### 4. One intensity front (`BlockSource`)

`BlockSource` erases any reader that models `concepts::BlockReader` (filterbank, time series, or an in-memory block). `read_block` copies packed bytes. `skip` moves by time samples, including backward.

```cpp
#include <psrio/psrio.hpp>

#include <vector>

int main() {
    psrio::MemoryInfo info;
    info.nchans = 4;
    info.nbits  = 8;
    info.tsamp  = 64.0e-6;
    info.fch1   = 1400.0;
    info.foff   = -0.5;

    std::vector<std::byte> packed(info.nchans * 2);
    psrio::BlockSource source(
        psrio::MemoryBlock(info, std::move(packed)));

    std::vector<std::byte> block(source.bytes_per_sample() * 2);
    const auto nread = source.read_block(2, block);
    if (nread > 0) {
        source.skip(-1); // one sample back, still inside this block
    }
    return 0;
}
```

### 5. Redigitise in memory (`digifil::Source`)

`digifil::Source` runs dspsr's `digifil` (rescale and 8-bit digitisation, optional `tscrunch`/`fscrunch`) as samples are read, over any `BlockSource`, including several contiguous PSRFITS files. It is itself a `BlockReader`, so no intermediate `.fil` is written. `header()` is the SIGPROC header `digifil` would write, and `refdm()` keeps the coherent DM (`CHAN_DM`).

```cpp
#include <psrio/psrio.hpp>

#include <filesystem>
#include <vector>

int main() {
    const std::vector<std::filesystem::path> files{"obs_0001.sf", "obs_0002.sf"};
    psrio::digifil::Config config; // 8 bits, 10 s rescale interval
    config.rescale_clip = 3.5;     // optional: robust to bursts of RFI
    psrio::digifil::Source stream(psrio::open_block_source(files), config);

    std::vector<std::uint8_t> block(4096 * stream.nchans());
    while (stream.read_samples(4096, std::span<std::uint8_t>(block)) > 0) {
        // channels descend in frequency; stream.refdm() is CHAN_DM
    }
    return 0;
}
```

### 6. Clean RFI, then redigitise, in memory (`stats`, `clean::Source`)

`stats::compute` gives per-channel mean, variance, skewness, and kurtosis in one pass. `clean::channel_mask` flags outlying channels as `spp_clean rfi` does, and `clean::Source` replaces them while reading. Like `digifil::Source`, it is a `BlockReader`, so the two chain without a file in between.

```cpp
#include <psrio/psrio.hpp>

#include <filesystem>
#include <utility>
#include <vector>

int main() {
    const std::vector<std::filesystem::path> files{"obs_0001.sf", "obs_0002.sf"};
    auto input = psrio::open_block_source(files);

    // Statistics from the first 2^18 samples, then a double-MAD mask.
    const auto stats = psrio::stats::compute(
        input, psrio::stats::ComputeOptions{.nsamps = std::uint64_t{1} << 18U});
    auto mask = psrio::clean::channel_mask(stats, psrio::clean::Method::kMad, 3.0);
    const auto fill = static_cast<float>(psrio::clean::fill_value(stats, mask));

    input.rewind();
    psrio::clean::Source cleaned(
        input, psrio::clean::Steps{.mask = std::move(mask), .fill = fill});
    psrio::digifil::Source stream(psrio::BlockSource(std::move(cleaned)),
                                  psrio::digifil::Config{});
    std::vector<std::uint8_t> block(4096 * stream.nchans());
    while (stream.read_samples(4096, std::span<std::uint8_t>(block)) > 0) {
    }
    return 0;
}
```

### 7. GUPPI RAW bytes

`RawReader` walks one header and the payload that follows it. It does not model the block fronts.

```cpp
#include <psrio/formats/guppi.hpp>

#include <iostream>

int main() {
    psrio::formats::guppi::RawReader reader("observation.raw");
    const auto header = reader.read_header();
    const auto bytes  = header.get<std::int64_t>("BLOCSIZE");
    const auto payload =
        reader.view_bytes(static_cast<std::uint64_t>(bytes));
    std::cout << "payload bytes: " << payload.size() << "\n";
    return 0;
}
```

### 8. Baseband block front

GUPPI RAW and DADA files model `concepts::BasebandReader`. `BasebandSource` erases either one. `read_block` copies canonical time-major voltages (`[time][antenna][channel][polarization][component]`). `read_samples` unpacks them into a caller buffer.

`GuppiReader::open` takes every file in a Breakthrough Listen set. Files that share `OBSFREQ` are one band continued in time. Files with different `OBSFREQ` are stitched into one channel axis, ordered by the lowest channel frequency. Each file keeps its own channel order.

```cpp
#include <psrio/psrio.hpp>

#include <complex>
#include <filesystem>
#include <iostream>
#include <vector>

int main() {
    const std::filesystem::path files[] = {
        "blc00_guppi_obs.0000.raw",
        "blc01_guppi_obs.0000.raw",
    };
    psrio::GuppiSet guppi = psrio::GuppiReader::open(files);
    std::cout << "channels " << guppi.nchan()
              << " stride " << guppi.bytes_per_sample() << "\n";

    std::vector<std::complex<float>> spectra(guppi.nchan() * guppi.npol());
    guppi.read_samples(1, spectra);

    psrio::DadaReader dada("observation.dada");
    psrio::BasebandSource source(std::move(dada));
    auto volts = source.read_samples<std::complex<float>>(1024);

    // CohFDMT wants FTPRI and channel 0 at the low edge of the band.
    // 4-bit and 2-bit bytes stay least-significant-field first (`msb_first = false`).
    std::vector<std::byte> coherent(guppi.bytes_per_sample() * 1024);
    guppi.seek(0);
    guppi.read_voltages(1024, {}, coherent);
    return 0;
}
```

---

## The `psrio` command

Built with `-DPSRIO_BUILD_APPLICATIONS=ON` (enabled by default for top-level builds). Every command accepts one or more `.fil` or search-mode `.sf` files (automatically sorted chronologically by start time), quoted globs (e.g. `"/data/obs/*.sf"`), `@filelist.txt` / `--file-list`, and `-j` threads.

### Installing the CLI

```bash
# Clone and build an optimised Release binary (defaults to -march=native)
git clone https://github.com/pravirkr/psrio.git
cd psrio
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DPSRIO_BUILD_APPLICATIONS=ON \
  -DPSRIO_BUILD_TESTS=OFF \
  -DCMAKE_INSTALL_PREFIX="$HOME/.local"
cmake --build build -j
cmake --install build
```

Make sure `$HOME/.local/bin` is on your `PATH`, then verify with `psrio --help`.

Useful CLI CMake options:

- `-DPSRIO_CLI_CPU_ARCH=<arch>` — `-march` target for the CLI (defaults to `native`). On heterogeneous HPC clusters where login and compute nodes have different CPUs, set a portable target such as `-DPSRIO_CLI_CPU_ARCH=x86-64-v3` (AVX2) to avoid illegal-instruction crashes on older nodes.
- `-DPSRIO_CLI_WITH_DMT=OFF` — skip fetching `dmt` and requiring FFTW if you do not need `psrio dedisperse` (faster build for `digifil`, `inspect`, `slice`, `clean`, etc.).
- `-DPSRIO_DMT_GPU=AUTO|CUDA|HIP|OFF` — GPU backend for `psrio dedisperse --device`. When building without a visible GPU, also pass `-DDMT_CUDA_ARCHITECTURES=<compute capability>` (e.g. `89` for an L40); `dmt`'s default `native` falls back to `sm_52`.
- `-DPSRIO_WITH_HDF5=ON` — enable reading `.h5` (FBH5) files in the CLI.

### Commands

| Command | sigpyproc / dspsr equivalent |
|---|---|
| `inspect` (alias `header`), `header FILE --set key=value` | `spp_header print/get/update` |
| `slice`, `splice`, `decimate`, `convert` | `spp_io extract/splice`, `spp_decimate` |
| `digifil` (`--scloffs`, `-P`, `-U`, `--rescale-clip`, `--progress`, `--io-buffer-mb`) | `spp_digifil`, dspsr `digifil` |
| `stats`, `bandpass` | `spp_stats`, `spp_bandpass` |
| `clean rfi/zap-channels/zap-times/zerodm` | `spp_clean` |
| `dedisperse` (`--device cpu\|cuda[:N]`, needs dmt) | `spp_dedisperse` |

### HPC / Large-Observation Example

```bash
# Convert 500+ search-mode PSRFITS files into one filterbank on an HPC node
psrio -v digifil \
  "/data/obs/*.sf" \
  -o /scratch/output.fil \
  -b 8 -t 4 -f 2 -I 10 \
  -j 8 --io-buffer-mb 256 --drop-cache
```

- `--progress auto` (default) shows an interactive progress bar on a terminal and switches to periodic newline log lines when redirected to a Slurm/PBS log (`-q` or `--progress none` disables it).
- `--io-buffer-mb` tunes background read-ahead memory (default `128` MiB per file) without altering `-U` rescaling block boundaries.
- `--drop-cache` (Linux) asks the kernel to evict consumed input and written output pages from the page cache to stay within strict Slurm/cgroup memory limits.

---

## Building and Running Tests Locally

To build and run the Catch2 test suite:

```bash
cmake -S . -B build -DPSRIO_BUILD_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

---

## License

`psrio` is licensed under the [MIT License](LICENSE).
