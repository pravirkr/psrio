# psrio

A high-performance, header-only C++20 library for reading, parsing, and streaming time-domain radio astronomy data (pulsar baseband, filterbank, and single-channel time series).

`psrio` provides zero-overhead, type-safe ingestion of observational data formats (SIGPROC `.fil` / `.tim`, PRESTO `.dat` / `.inf`, and GUPPI RAW). Intensity formats share one block-streaming front. GUPPI stays on its own baseband front.

## Key Features

- **Pure Header-Only C++20**: Zero compiled library binaries. Simply link `psrio::psrio` and `#include <psrio/psrio.hpp>`.
- **Zero-Copy & Zero-Allocation**: Direct memory-mapped access (`mmap`) to file payloads. Readers unpack directly into caller-managed buffers (`std::span<float>`, pinned CUDA host memory, etc.) without intermediate heap allocations.
- **Unified Metadata Representation**: A single flat, idiomatic [`psrio::Header`](include/psrio/header.hpp) struct.
- **Zero Third-Party Dependencies by Default**: Base formats require only the ISO C++ standard library and POSIX system APIs (`<sys/mman.h>`). FBH5 and PSRDADA are optional and stay out of `psrio.hpp` unless `PSRIO_WITH_HDF5` or `PSRIO_WITH_PSRDADA` is turned on.
- **Block streaming**: `psrio::concepts::BlockReader` and the move-only `psrio::BlockSource` read packed bytes (`read_block`) or unpacked floats (`read_samples`), and move the sample cursor with `seek`, `rewind`, and `skip`.

---

## Requirements

- **C++ Standard**: C++20 (`-std=c++20`)
- **Compilers**: GCC 13+, Clang 18+
- **Build System**: CMake 3.21+

---

## Adding `psrio` to Your CMake Project

`psrio` exports the interface target `psrio::psrio`. Linking to this target automatically configures include directories and enables C++20 compilation features on your target.

### Method 1: Using `CPM.cmake` (Recommended)

Add `psrio` directly from GitHub using [CPM.cmake](https://github.com/cpm-cmake/CPM.cmake):

```cmake
CPMAddPackage(
  NAME psrio
  GITHUB_REPOSITORY pravirkr/psrio
  GIT_TAG v0.2.0  # Or a specific branch/commit
)

target_link_libraries(my_downstream_target PRIVATE psrio::psrio)
```

### Method 2: Using CMake `FetchContent`

Using built-in CMake 3.21+ `FetchContent`:

```cmake
include(FetchContent)

FetchContent_Declare(
  psrio
  GIT_REPOSITORY https://github.com/pravirkr/psrio.git
  GIT_TAG v0.2.0  # Or a specific branch/commit
)

FetchContent_MakeAvailable(psrio)

target_link_libraries(my_downstream_target PRIVATE psrio::psrio)
```

### Method 3: Using `find_package` (Pre-installed)

First install `psrio` into your system or a custom prefix:

```bash
cmake -B build -DPSRIO_BUILD_TESTS=OFF -DCMAKE_INSTALL_PREFIX=/usr/local
sudo cmake --install build
```

Then in your project's `CMakeLists.txt`:

```cmake
find_package(psrio CONFIG REQUIRED)

target_link_libraries(my_downstream_target PRIVATE psrio::psrio)
```


> [!NOTE]
> `PSRIO_BUILD_TESTS` defaults to `OFF` automatically whenever `psrio` is consumed as a dependency via `CPM`, `FetchContent`, or `add_subdirectory`.

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

### 5. GUPPI RAW (baseband, not a block source)

GUPPI complex voltages do not model `BlockReader`. Read a header, then a payload span:

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
