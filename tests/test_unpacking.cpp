#include "psrio/detail/endian.hpp"
#include "psrio/detail/unpack.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

TEST_CASE("endian swap round-trips 32-bit values", "[endian]") {
    constexpr std::uint32_t kValue = 0x01020304U;
    REQUIRE(psrio::detail::swap_endian(psrio::detail::swap_endian(kValue)) ==
            kValue);
    REQUIRE(psrio::detail::swap_endian(kValue) == 0x04030201U);
}

TEST_CASE("little-endian loads keep host values on this machine", "[endian]") {
    if (psrio::detail::host_is_little_endian()) {
        REQUIRE(psrio::detail::from_little_endian(std::uint32_t{0x01020304U}) ==
                0x01020304U);
        REQUIRE(psrio::detail::from_little_endian(1.5F) == 1.5F);
    }
}

TEST_CASE("8-bit unpack maps bytes to floats", "[unpack]") {
    const std::array<std::uint8_t, 4> input{0, 1, 127, 255};
    std::array<float, 4> output{};
    psrio::detail::unpack_8bit_to_float(input, output);
    REQUIRE(output.at(0) == 0.0F);
    REQUIRE(output.at(1) == 1.0F);
    REQUIRE(output.at(2) == 127.0F);
    REQUIRE(output.at(3) == 255.0F);
}

TEST_CASE("8-bit unpack rejects a short destination", "[unpack]") {
    const std::array<std::uint8_t, 2> input{1, 2};
    std::array<float, 1> output{};
    REQUIRE_THROWS_AS(psrio::detail::unpack_8bit_to_float(input, output),
                      std::invalid_argument);
}

TEST_CASE("sub-byte unpack is least-significant field first", "[unpack]") {
    const std::array<std::byte, 1> packed{std::byte{0x81}};
    std::array<std::uint8_t, 8> bits{};
    psrio::detail::unpack_lsb<std::uint8_t>(packed, bits, 1);
    REQUIRE(bits == std::array<std::uint8_t, 8>{1, 0, 0, 0, 0, 0, 0, 1});

    const std::array<std::byte, 1> two_bit{std::byte{0xE4}};
    std::array<std::uint8_t, 4> pairs{};
    psrio::detail::unpack_lsb<std::uint8_t>(two_bit, pairs, 2);
    REQUIRE(pairs == std::array<std::uint8_t, 4>{0, 1, 2, 3});

    const std::array<std::byte, 1> four_bit{std::byte{0xAB}};
    std::array<float, 2> nibbles{};
    psrio::detail::unpack_lsb<float>(four_bit, nibbles, 4);
    REQUIRE(nibbles.at(0) == 11.0F);
    REQUIRE(nibbles.at(1) == 10.0F);
}

TEST_CASE("multi-byte unpack honours little-endian layout", "[unpack]") {
    const std::array<std::byte, 2> raw16{std::byte{0x02}, std::byte{0x01}};
    std::array<std::uint16_t, 1> as_int{};
    std::array<float, 1> as_float{};
    psrio::detail::unpack_16le<std::uint16_t>(raw16, as_int);
    psrio::detail::unpack_16le<float>(raw16, as_float);
    REQUIRE(as_int.at(0) == 0x0102U);
    REQUIRE(as_float.at(0) == 258.0F);

    const float sample = -2.5F;
    std::array<std::byte, 4> raw32{};
    const auto bits =
        psrio::detail::from_little_endian(std::bit_cast<std::uint32_t>(sample));
    raw32.at(0) = std::byte{static_cast<unsigned char>(bits & 0xFFU)};
    raw32.at(1) = std::byte{static_cast<unsigned char>((bits >> 8U) & 0xFFU)};
    raw32.at(2) = std::byte{static_cast<unsigned char>((bits >> 16U) & 0xFFU)};
    raw32.at(3) = std::byte{static_cast<unsigned char>((bits >> 24U) & 0xFFU)};
    std::array<float, 1> floats{};
    psrio::detail::unpack_32le(raw32, floats);
    REQUIRE(floats.at(0) == -2.5F);
}

TEST_CASE("sub-byte unpack supports most-significant field first (MSB)",
          "[unpack]") {
    // 0x80 = 1000 0000 in binary. MSB first: first sample is bit 7 (1), rest 0.
    const std::array<std::byte, 1> packed{std::byte{0x80}};
    std::array<std::uint8_t, 8> bits{};
    psrio::detail::unpack_sub_byte<std::uint8_t>(packed, bits, 1,
                                                 psrio::BitOrder::kMsbFirst);
    REQUIRE(bits == std::array<std::uint8_t, 8>{1, 0, 0, 0, 0, 0, 0, 0});

    // 0xE4 = 11 10 01 00 in binary. MSB first: samples are 3, 2, 1, 0.
    const std::array<std::byte, 1> two_bit{std::byte{0xE4}};
    std::array<std::uint8_t, 4> pairs{};
    psrio::detail::unpack_sub_byte<std::uint8_t>(two_bit, pairs, 2,
                                                 psrio::BitOrder::kMsbFirst);
    REQUIRE(pairs == std::array<std::uint8_t, 4>{3, 2, 1, 0});

    // 0xAB = high nibble 10, low nibble 11. MSB first: 10, then 11.
    const std::array<std::byte, 1> four_bit{std::byte{0xAB}};
    std::array<float, 2> nibbles{};
    psrio::detail::unpack_sub_byte<float>(four_bit, nibbles, 4,
                                          psrio::BitOrder::kMsbFirst);
    REQUIRE(nibbles.at(0) == 10.0F);
    REQUIRE(nibbles.at(1) == 11.0F);
}

TEST_CASE("sub-byte pack and unpack round-trip across all bit orders",
          "[unpack]") {
    for (const auto order :
         {psrio::BitOrder::kLsbFirst, psrio::BitOrder::kMsbFirst}) {
        // 1-bit test
        const std::array<std::uint8_t, 8> input1{1, 0, 1, 1, 0, 0, 1, 0};
        std::array<std::byte, 1> packed1{};
        psrio::detail::pack_sub_byte(input1, packed1, 1, order);
        std::array<std::uint8_t, 8> unpacked1{};
        psrio::detail::unpack_sub_byte<std::uint8_t>(packed1, unpacked1, 1,
                                                     order);
        REQUIRE(input1 == unpacked1);

        // 2-bit test
        const std::array<std::uint8_t, 4> input2{0, 3, 1, 2};
        std::array<std::byte, 1> packed2{};
        psrio::detail::pack_sub_byte(input2, packed2, 2, order);
        std::array<std::uint8_t, 4> unpacked2{};
        psrio::detail::unpack_sub_byte<std::uint8_t>(packed2, unpacked2, 2,
                                                     order);
        REQUIRE(input2 == unpacked2);

        // 4-bit test
        const std::array<std::uint8_t, 2> input4{14, 5};
        std::array<std::byte, 1> packed4{};
        psrio::detail::pack_sub_byte(input4, packed4, 4, order);
        std::array<std::uint8_t, 2> unpacked4{};
        psrio::detail::unpack_sub_byte<std::uint8_t>(packed4, unpacked4, 4,
                                                     order);
        REQUIRE(input4 == unpacked4);
    }
}

TEST_CASE("signed 8-bit unpack casts through int8", "[unpack]") {
    const std::array<std::byte, 2> packed{std::byte{0xFF}, std::byte{0x05}};
    std::array<float, 2> output{};
    psrio::detail::unpack_8bit<float>(packed, output, true);
    REQUIRE(output.at(0) == -1.0F);
    REQUIRE(output.at(1) == 5.0F);

    std::array<std::uint8_t, 2> raw{};
    psrio::detail::unpack_8bit<std::uint8_t>(packed, raw, true);
    REQUIRE(raw.at(0) == 255U);
    REQUIRE(raw.at(1) == 5U);
}
