# psrio

Header-only C++20 library for reading and streaming time-domain radio astronomy data (SIGPROC filterbank/time series, PRESTO, GUPPI, and optional PSRFITS/FBH5).

## Requirements

- CMake 3.18+
- GCC 13+ or Clang 18+

## Build and test

```bash
cmake -S . -B build -DPSRIO_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

## Use in another CMake project

```cmake
add_subdirectory(path/to/psrio)
target_link_libraries(your_target PRIVATE psrio::psrio)
```

When `psrio` is not the top-level project, tests are off by default; enable with `-DPSRIO_BUILD_TESTS=ON`.

## License

MIT — see [LICENSE](LICENSE).
