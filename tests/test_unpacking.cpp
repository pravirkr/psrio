#include <catch2/catch_test_macros.hpp>

#include <psrio/detail/endian.hpp>
#include <psrio/detail/unpack.hpp>

#include <array>
#include <cstdint>

TEST_CASE("endian swap round-trips 32-bit values", "[endian]")
{
    constexpr std::uint32_t kValue = 0x01020304U;
    REQUIRE(psrio::detail::swap_endian(psrio::detail::swap_endian(kValue)) == kValue);
    REQUIRE(psrio::detail::swap_endian(kValue) == 0x04030201U);
}

TEST_CASE("8-bit unpack maps bytes to floats", "[unpack]")
{
    const std::array<std::uint8_t, 4> input{0, 1, 127, 255};
    std::array<float, 4> output{};
    psrio::detail::unpack_8bit_to_float(input, output);
    REQUIRE(output.at(0) == 0.0F);
    REQUIRE(output.at(1) == 1.0F);
    REQUIRE(output.at(2) == 127.0F);
    REQUIRE(output.at(3) == 255.0F);
}
