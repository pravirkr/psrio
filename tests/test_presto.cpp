#include "psrio/detail/mmap.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>

TEST_CASE("MappedFile exposes file bytes", "[presto][mmap]") {
    constexpr std::size_t kPayloadSize = 5;
    const auto path =
        std::filesystem::temp_directory_path() / "psrio_mmap_smoke.bin";
    {
        std::ofstream out(path, std::ios::binary);
        REQUIRE(out);
        out.write("psrio", static_cast<std::streamsize>(kPayloadSize));
    }

    const psrio::detail::MappedFile mapped(path);
    REQUIRE(mapped.size() == kPayloadSize);
    const auto bytes = mapped.bytes();
    REQUIRE(bytes.size() == kPayloadSize);
    REQUIRE(bytes.front() == std::byte{'p'});

    std::filesystem::remove(path);
}
