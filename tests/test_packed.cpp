#include "psrio/detail/exceptions.hpp"
#include "psrio/packed.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstdint>
#include <vector>

TEST_CASE("sample_type predicates and helpers", "[sample_type]") {
    REQUIRE(psrio::sample_type_from_nbits(1, false) == psrio::SampleType::kUInt1);
    REQUIRE(psrio::sample_type_from_nbits(2, false) == psrio::SampleType::kUInt2);
    REQUIRE(psrio::sample_type_from_nbits(4, false) == psrio::SampleType::kUInt4);
    REQUIRE(psrio::sample_type_from_nbits(8, false) == psrio::SampleType::kUInt8);
    REQUIRE(psrio::sample_type_from_nbits(8, true) == psrio::SampleType::kInt8);
    REQUIRE(psrio::sample_type_from_nbits(16, false) ==
            psrio::SampleType::kUInt16);
    REQUIRE(psrio::sample_type_from_nbits(32, false) ==
            psrio::SampleType::kFloat32);
    REQUIRE(psrio::bits_of(psrio::SampleType::kUInt32) == 32);
    REQUIRE(psrio::bits_of(psrio::SampleType::kFloat32) == 32);

    REQUIRE(psrio::is_floating(psrio::SampleType::kFloat32));
    REQUIRE_FALSE(psrio::is_floating(psrio::SampleType::kUInt32));
    REQUIRE_FALSE(psrio::is_floating(psrio::SampleType::kUInt8));

    REQUIRE(psrio::is_integral(psrio::SampleType::kUInt8));
    REQUIRE(psrio::is_integral(psrio::SampleType::kInt8));
    REQUIRE_FALSE(psrio::is_integral(psrio::SampleType::kFloat32));

    REQUIRE(psrio::is_signed(psrio::SampleType::kInt8));
    REQUIRE(psrio::is_signed(psrio::SampleType::kFloat32));
    REQUIRE_FALSE(psrio::is_signed(psrio::SampleType::kUInt8));

    REQUIRE(psrio::to_nbits(psrio::SampleType::kFloat32) == -32);
    REQUIRE(psrio::to_nbits(psrio::SampleType::kUInt8) == 8);
    REQUIRE(psrio::to_nbits(psrio::SampleType::kUInt1) == 1);
    REQUIRE(psrio::to_nbits(psrio::SampleType::kUInt2) == 2);
    REQUIRE(psrio::to_nbits(psrio::SampleType::kUInt4) == 4);
    REQUIRE(psrio::to_nbits(psrio::SampleType::kUInt16) == 16);
    REQUIRE(psrio::to_nbits(psrio::SampleType::kUInt32) == 32);

    REQUIRE_THROWS_AS(psrio::sample_type_from_nbits(32, true),
                      psrio::ValidationError);
    REQUIRE_THROWS_AS(psrio::sample_type_from_nbits(12, false),
                      psrio::ValidationError);
}

TEST_CASE("reverse_channels reverses 8-bit channel order", "[packed]") {
    std::vector<std::byte> gulp{
        std::byte{10}, std::byte{20}, std::byte{30}, std::byte{40}};
    psrio::reverse_channels(gulp, 1, 4, 8);
    REQUIRE(gulp == std::vector<std::byte>{std::byte{40}, std::byte{30},
                                           std::byte{20}, std::byte{10}});
}

TEST_CASE("reverse_channels reverses 1-bit, 2-bit, and 4-bit rows accurately",
          "[packed]") {
    // 4-bit: 8 channels across 4 bytes.
    // In nibbles (low to high per byte):
    // byte 0: 1, 2 (0x21)
    // byte 1: 3, 4 (0x43)
    // byte 2: 5, 6 (0x65)
    // byte 3: 7, 8 (0x87)
    // Reversing all 8 channels yields: 8, 7, 6, 5, 4, 3, 2, 1
    // byte 0: 8, 7 -> 0x78
    // byte 1: 6, 5 -> 0x56
    // byte 2: 4, 3 -> 0x34
    // byte 3: 2, 1 -> 0x12
    const std::vector<std::byte> four_bit{
        std::byte{0x21}, std::byte{0x43}, std::byte{0x65}, std::byte{0x87}};
    auto flipped4 = four_bit;
    psrio::reverse_channels(flipped4, 1, 8, 4);
    REQUIRE(flipped4 == std::vector<std::byte>{std::byte{0x78}, std::byte{0x56},
                                               std::byte{0x34}, std::byte{0x12}});

    // 2-bit: 8 channels across 2 bytes.
    // byte 0: ch0=1, ch1=0, ch2=2, ch3=0 -> 0b00100001 = 0x21
    // byte 1: ch4=3, ch5=0, ch6=0, ch7=1 -> 0b01000011 = 0x43
    // Reversed (ch7..ch0):
    // byte 0: ch7=1, ch6=0, ch5=0, ch4=3 -> 0b11000001 = 0xC1
    // byte 1: ch3=0, ch2=2, ch1=0, ch0=1 -> 0b01001000 = 0x48
    const std::vector<std::byte> two_bit{std::byte{0x21}, std::byte{0x43}};
    auto flipped2 = two_bit;
    psrio::reverse_channels(flipped2, 1, 8, 2);
    REQUIRE(flipped2 ==
            std::vector<std::byte>{std::byte{0xC1}, std::byte{0x48}});

    // 1-bit: 8 channels in 1 byte (odd row length of 1 byte!)
    // bits 0..7: 1, 0, 1, 0, 0, 0, 1, 1 -> 0b11000101 = 0xC5
    // reversed:  1, 1, 0, 0, 0, 1, 0, 1 -> 0b10100011 = 0xA3
    std::vector<std::byte> one_bit{std::byte{0xC5}};
    psrio::reverse_channels(one_bit, 1, 8, 1);
    REQUIRE(one_bit == std::vector<std::byte>{std::byte{0xA3}});
}

TEST_CASE("reverse_channels reverses odd-length rows including middle byte",
          "[packed]") {
    // 4-bit with 6 channels = 3 bytes (odd row_bytes = 3).
    // byte 0: ch0=1, ch1=2 -> 0x21
    // byte 1 (middle): ch2=3, ch3=4 -> 0x43 (must be swapped to 0x34!)
    // byte 2: ch4=5, ch5=6 -> 0x65
    // Reversed channels:
    // byte 0: ch5=6, ch4=5 -> 0x56
    // byte 1: ch3=4, ch2=3 -> 0x34
    // byte 2: ch1=2, ch0=1 -> 0x12
    std::vector<std::byte> row3{std::byte{0x21}, std::byte{0x43}, std::byte{0x65}};
    psrio::reverse_channels(row3, 1, 6, 4);
    REQUIRE(row3 == std::vector<std::byte>{std::byte{0x56}, std::byte{0x34},
                                           std::byte{0x12}});
}

TEST_CASE("reverse_channels supports 16-bit and 32-bit", "[packed]") {
    std::vector<std::byte> row16{
        std::byte{1}, std::byte{0}, std::byte{2}, std::byte{0},
        std::byte{3}, std::byte{0}};
    psrio::reverse_channels(row16, 1, 3, 16);
    REQUIRE(row16 == std::vector<std::byte>{
        std::byte{3}, std::byte{0}, std::byte{2}, std::byte{0},
        std::byte{1}, std::byte{0}});

    std::vector<std::byte> row32{
        std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0},
        std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}};
    psrio::reverse_channels(row32, 1, 2, 32);
    REQUIRE(row32 == std::vector<std::byte>{
        std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0},
        std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}});
}

TEST_CASE("reverse_channels is its own inverse across bit depths",
          "[packed]") {
    for (const int nbits : {1, 2, 4, 8, 16, 32}) {
        const auto nchans = 32U;
        const auto row_bytes = (nchans * static_cast<std::uint64_t>(nbits)) / 8U;
        std::vector<std::byte> original(row_bytes);
        for (std::size_t i = 0; i < original.size(); ++i) {
            original[i] = static_cast<std::byte>((i * 37U + 11U) & 0xFFU);
        }
        auto gulp = original;
        psrio::reverse_channels(gulp, 1, nchans, nbits);
        REQUIRE(gulp != original);
        psrio::reverse_channels(gulp, 1, nchans, nbits);
        REQUIRE(gulp == original);
    }
}

TEST_CASE("DataExtractor extracts samples and spectra across bit depths",
          "[packed]") {
    // 4-bit data: 2 samples, 4 channels.
    // sample 0: channels 1, 2, 3, 4 -> byte 0: 0x21, byte 1: 0x43
    // sample 1: channels 5, 6, 7, 8 -> byte 2: 0x65, byte 3: 0x87
    const std::vector<std::byte> bytes4{
        std::byte{0x21}, std::byte{0x43}, std::byte{0x65}, std::byte{0x87}};
    psrio::DataExtractor ext4(bytes4, 2, 4, psrio::SampleType::kUInt4);
    REQUIRE(ext4.extract_sample(0, 0) == 1.0F);
    REQUIRE(ext4.extract_sample(0, 1) == 2.0F);
    REQUIRE(ext4.extract_sample(0, 2) == 3.0F);
    REQUIRE(ext4.extract_sample(0, 3) == 4.0F);
    REQUIRE(ext4.extract_sample(1, 0) == 5.0F);
    REQUIRE(ext4.extract_sample(1, 3) == 8.0F);

    std::vector<float> spectrum(4);
    ext4.extract_spectrum(1, spectrum);
    REQUIRE(spectrum == std::vector<float>{5.0F, 6.0F, 7.0F, 8.0F});

    // 8-bit signed
    const std::vector<std::byte> bytes_s8{std::byte{0xFE}, std::byte{0x05}};
    psrio::DataExtractor ext_s8(bytes_s8, 1, 2, psrio::SampleType::kInt8);
    REQUIRE(ext_s8.extract_sample(0, 0) == -2.0F);
    REQUIRE(ext_s8.extract_sample(0, 1) == 5.0F);

    // Float32
    float floats[2] = {12.5F, -3.25F};
    std::span<const std::byte> float_bytes(
        reinterpret_cast<const std::byte*>(floats), sizeof(floats));
    psrio::DataExtractor ext_f(float_bytes, 1, 2, psrio::SampleType::kFloat32);
    REQUIRE(ext_f.extract_sample(0, 0) == 12.5F);
    REQUIRE(ext_f.extract_sample(0, 1) == -3.25F);

    REQUIRE_THROWS_AS(ext4.extract_sample(2, 0), psrio::ValidationError);
    REQUIRE_THROWS_AS(ext4.extract_sample(0, 4), psrio::ValidationError);
}

TEST_CASE("write_sample modifies packed and unpacked samples in place",
          "[packed]") {
    std::vector<std::byte> buffer(4, std::byte{0});

    // 4-bit write
    psrio::write_sample(buffer, 4, 0, 0, psrio::SampleType::kUInt4, 5.0F);
    psrio::write_sample(buffer, 4, 0, 1, psrio::SampleType::kUInt4, 9.0F);
    // byte 0 has ch0=5 (low nibble), ch1=9 (high nibble) -> 0x95
    REQUIRE(buffer[0] == std::byte{0x95});

    // 8-bit write
    psrio::write_sample(buffer, 4, 0, 2, psrio::SampleType::kUInt8, 42.0F);
    REQUIRE(buffer[2] == std::byte{42});

    // Verify round-trip with DataExtractor
    psrio::DataExtractor ext(buffer, 1, 4, psrio::SampleType::kUInt4);
    REQUIRE(ext.extract_sample(0, 0) == 5.0F);
    REQUIRE(ext.extract_sample(0, 1) == 9.0F);
}
