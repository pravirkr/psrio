#include "psrio/detail/exceptions.hpp"
#include "psrio/packed.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstdint>
#include <vector>

TEST_CASE("sample_type predicates and helpers", "[sample_type]") {
    REQUIRE(psrio::sample_type_from_nbits(1, false) ==
            psrio::SampleType::kUInt1);
    REQUIRE(psrio::sample_type_from_nbits(2, false) ==
            psrio::SampleType::kUInt2);
    REQUIRE(psrio::sample_type_from_nbits(4, false) ==
            psrio::SampleType::kUInt4);
    REQUIRE(psrio::sample_type_from_nbits(8, false) ==
            psrio::SampleType::kUInt8);
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
    REQUIRE(psrio::is_signed(psrio::SampleType::kInt16));
    REQUIRE(psrio::is_signed(psrio::SampleType::kFloat32));
    REQUIRE_FALSE(psrio::is_signed(psrio::SampleType::kUInt8));
    REQUIRE(psrio::bits_of(psrio::SampleType::kInt16) == 16);
    REQUIRE(psrio::to_nbits(psrio::SampleType::kInt16) == 16);
    REQUIRE(psrio::sample_type_name(psrio::SampleType::kInt16) == "int16");

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
    std::vector<std::byte> gulp{std::byte{10}, std::byte{20}, std::byte{30},
                                std::byte{40}};
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
    const std::vector<std::byte> four_bit{std::byte{0x21}, std::byte{0x43},
                                          std::byte{0x65}, std::byte{0x87}};
    auto flipped4 = four_bit;
    psrio::reverse_channels(flipped4, 1, 8, 4);
    REQUIRE(flipped4 == std::vector<std::byte>{std::byte{0x78}, std::byte{0x56},
                                               std::byte{0x34},
                                               std::byte{0x12}});

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
    std::vector<std::byte> row3{std::byte{0x21}, std::byte{0x43},
                                std::byte{0x65}};
    psrio::reverse_channels(row3, 1, 6, 4);
    REQUIRE(row3 == std::vector<std::byte>{std::byte{0x56}, std::byte{0x34},
                                           std::byte{0x12}});
}

TEST_CASE("reverse_channels supports 16-bit and 32-bit", "[packed]") {
    std::vector<std::byte> row16{std::byte{1}, std::byte{0}, std::byte{2},
                                 std::byte{0}, std::byte{3}, std::byte{0}};
    psrio::reverse_channels(row16, 1, 3, 16);
    REQUIRE(row16 == std::vector<std::byte>{std::byte{3}, std::byte{0},
                                            std::byte{2}, std::byte{0},
                                            std::byte{1}, std::byte{0}});

    std::vector<std::byte> row32{std::byte{1}, std::byte{0}, std::byte{0},
                                 std::byte{0}, std::byte{2}, std::byte{0},
                                 std::byte{0}, std::byte{0}};
    psrio::reverse_channels(row32, 1, 2, 32);
    REQUIRE(row32 == std::vector<std::byte>{std::byte{2}, std::byte{0},
                                            std::byte{0}, std::byte{0},
                                            std::byte{1}, std::byte{0},
                                            std::byte{0}, std::byte{0}});
}

TEST_CASE("reverse_channels is its own inverse across bit depths", "[packed]") {
    for (const int nbits : {1, 2, 4, 8, 16, 32}) {
        const auto nchans = 32U;
        const auto row_bytes =
            (nchans * static_cast<std::uint64_t>(nbits)) / 8U;
        std::vector<std::byte> original(row_bytes);
        for (std::size_t i = 0; i < original.size(); ++i) {
            original[i] = static_cast<std::byte>(((i * 37U) + 11U) & 0xFFU);
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
    const std::vector<std::byte> bytes4{std::byte{0x21}, std::byte{0x43},
                                        std::byte{0x65}, std::byte{0x87}};
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
    std::array<float, 2> floats = {12.5F, -3.25F};
    std::span<const std::byte> float_bytes(
        reinterpret_cast<const std::byte*>(floats.data()),
        floats.size() * sizeof(float));
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

TEST_CASE(
    "BitsInfo and DigitizationInfo provide standard digitization constants",
    "[packed][bits_info]") {
    // Compile-time constexpr validation
    constexpr psrio::BitsInfo kConstexpr1Bit(1);
    static_assert(kConstexpr1Bit.nbits() == 1);
    static_assert(kConstexpr1Bit.itemsize() == 1);
    static_assert(kConstexpr1Bit.can_pack_unpack());
    static_assert(kConstexpr1Bit.bit_factor() == 8);
    static_assert(psrio::BitsInfo::digi_min() == 0);
    static_assert(kConstexpr1Bit.digi_max() == 1);

    // 1-bit
    const psrio::BitsInfo b1(1);
    CHECK(b1.nbits() == 1);
    CHECK(b1.get_nbits() == 1);
    CHECK(b1.itemsize() == 1U);
    CHECK(b1.get_itemsize() == 1U);
    CHECK(b1.can_pack_unpack());
    CHECK(b1.get_can_pack_unpack());
    CHECK(b1.bit_factor() == 8U);
    CHECK(b1.bitfact() == 8U);
    CHECK(b1.get_bitfact() == 8U);
    CHECK(b1.digi_min() == 0U);
    CHECK(b1.get_digi_min() == 0U);
    CHECK(b1.digi_max() == 1U);
    CHECK(b1.get_digi_max() == 1U);
    CHECK(b1.digi_mean() == 0.5F);
    CHECK(b1.get_digi_mean() == 0.5F);
    CHECK(b1.digi_sigma() == 0.5F);
    CHECK(b1.get_digi_sigma() == 0.5F);
    CHECK(b1.digi_scale() == 1.0F);
    CHECK(b1.get_digi_scale() == 1.0F);

    // 2-bit
    const psrio::DigitizationInfo b2(2);
    CHECK(b2.nbits() == 2);
    CHECK(b2.itemsize() == 1U);
    CHECK(b2.can_pack_unpack());
    CHECK(b2.bit_factor() == 4U);
    CHECK(b2.digi_min() == 0U);
    CHECK(b2.digi_max() == 3U);
    CHECK(b2.digi_mean() == 1.5F);
    CHECK(b2.digi_sigma() == 1.5F);
    CHECK(b2.digi_scale() == 1.0F);

    // 4-bit
    const psrio::BitsInfo b4(4);
    CHECK(b4.nbits() == 4);
    CHECK(b4.itemsize() == 1U);
    CHECK(b4.can_pack_unpack());
    CHECK(b4.bit_factor() == 2U);
    CHECK(b4.digi_max() == 15U);
    CHECK(b4.digi_mean() == 7.5F);
    CHECK(b4.digi_sigma() == 6.0F);
    CHECK(b4.digi_scale() == 1.25F);

    // 8-bit
    const psrio::BitsInfo b8(8);
    CHECK(b8.nbits() == 8);
    CHECK(b8.itemsize() == 1U);
    CHECK_FALSE(b8.can_pack_unpack());
    CHECK(b8.bit_factor() == 1U);
    CHECK(b8.digi_max() == 255U);
    CHECK(b8.digi_mean() == 127.5F);
    CHECK(b8.digi_sigma() == 6.0F);
    CHECK(b8.digi_scale() == 127.5F / 6.0F);

    // 16-bit
    const psrio::BitsInfo b16(16);
    CHECK(b16.nbits() == 16);
    CHECK(b16.itemsize() == 2U);
    CHECK_FALSE(b16.can_pack_unpack());
    CHECK(b16.bit_factor() == 1U);
    CHECK(b16.digi_max() == 65535U);
    CHECK(b16.digi_mean() == 32767.5F);
    CHECK(b16.digi_sigma() == 6.0F);
    CHECK(b16.digi_scale() == 32767.5F / 6.0F);

    // 32-bit
    const psrio::BitsInfo b32(32);
    CHECK(b32.nbits() == 32);
    CHECK(b32.itemsize() == 4U);
    CHECK_FALSE(b32.can_pack_unpack());
    CHECK(b32.bit_factor() == 1U);
    CHECK(b32.digi_max() == 4294967295ULL);
    CHECK(b32.digi_mean() == 2147483647.5F);
    CHECK(b32.digi_sigma() == 6.0F);

    // Invalid bit depths
    CHECK_THROWS_AS(psrio::BitsInfo(0), psrio::ValidationError);
    CHECK_THROWS_AS(psrio::BitsInfo(3), psrio::ValidationError);
    CHECK_THROWS_AS(psrio::BitsInfo(12), psrio::ValidationError);
    CHECK_THROWS_AS(psrio::BitsInfo(64), psrio::ValidationError);
}

TEST_CASE("pack_sub_byte and unpack_sub_byte convenience allocating overloads",
          "[packed]") {
    for (const auto order :
         {psrio::BitOrder::kLsbFirst, psrio::BitOrder::kMsbFirst}) {
        // 1-bit
        const std::vector<std::uint8_t> in1{1, 0, 1, 1, 0, 0, 1, 0};
        const auto packed1 = psrio::pack_sub_byte(in1, 1, order);
        REQUIRE(packed1.size() == 1U);
        const auto unpacked1 =
            psrio::unpack_sub_byte<std::uint8_t>(packed1, 1, order);
        CHECK(unpacked1 == in1);

        // 2-bit
        const std::vector<std::uint8_t> in2{0, 3, 1, 2, 2, 1, 3, 0};
        const auto packed2 = psrio::pack_sub_byte(in2, 2, order);
        REQUIRE(packed2.size() == 2U);
        const auto unpacked2 =
            psrio::unpack_sub_byte<std::uint8_t>(packed2, 2, order);
        CHECK(unpacked2 == in2);

        // 4-bit
        const std::vector<std::uint8_t> in4{15, 0, 4, 11};
        const auto packed4 = psrio::pack_sub_byte(in4, 4, order);
        REQUIRE(packed4.size() == 2U);
        const auto unpacked4 =
            psrio::unpack_sub_byte<std::uint8_t>(packed4, 4, order);
        CHECK(unpacked4 == in4);
    }
}

TEST_CASE("pack_inplace and unpack_inplace round-trip within a single buffer",
          "[packed]") {
    for (const auto order :
         {psrio::BitOrder::kLsbFirst, psrio::BitOrder::kMsbFirst}) {
        // 1-bit: 16 samples pack down to 2 bytes, then unpack back to 16
        // samples
        const std::vector<std::uint8_t> original1 = {
            1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0,
        };
        std::vector<std::uint8_t> buf1 = original1;
        const auto packed_span1 =
            psrio::pack_inplace<std::uint8_t>(buf1, 1, order);
        REQUIRE(packed_span1.size() == 2U);
        // Now unpack in-place backwards
        psrio::unpack_inplace<std::uint8_t>(buf1, 1, order);
        CHECK(buf1 == original1);

        // 2-bit: 8 samples pack down to 2 bytes, then unpack back to 8 samples
        const std::vector<std::uint8_t> original2 = {0, 3, 1, 2, 2, 1, 0, 3};
        std::vector<std::uint8_t> buf2            = original2;
        const auto packed_span2 =
            psrio::pack_inplace<std::uint8_t>(buf2, 2, order);
        REQUIRE(packed_span2.size() == 2U);
        psrio::unpack_inplace<std::uint8_t>(buf2, 2, order);
        CHECK(buf2 == original2);

        // 4-bit: 4 samples pack down to 2 bytes, then unpack back to 4 samples
        const std::vector<std::uint8_t> original4 = {14, 5, 2, 9};
        std::vector<std::uint8_t> buf4            = original4;
        const auto packed_span4 =
            psrio::pack_inplace<std::uint8_t>(buf4, 4, order);
        REQUIRE(packed_span4.size() == 2U);
        psrio::unpack_inplace<std::uint8_t>(buf4, 4, order);
        CHECK(buf4 == original4);

        // Also verify std::byte overloads
        std::vector<std::byte> byte_buf(original4.size());
        for (std::size_t i = 0; i < original4.size(); ++i) {
            byte_buf[i] = static_cast<std::byte>(original4[i]);
        }
        const auto packed_bytes =
            psrio::pack_inplace<std::byte>(byte_buf, 4, order);
        REQUIRE(packed_bytes.size() == 2U);
        psrio::unpack_inplace<std::byte>(byte_buf, 4, order);
        for (std::size_t i = 0; i < original4.size(); ++i) {
            CHECK(byte_buf[i] == static_cast<std::byte>(original4[i]));
        }
    }

    // Invalid arguments
    std::vector<std::uint8_t> invalid_size{1, 2,
                                           3}; // Not multiple of 8 for 1-bit
    CHECK_THROWS_AS(psrio::pack_inplace<std::uint8_t>(invalid_size, 1),
                    psrio::ValidationError);
    CHECK_THROWS_AS(psrio::unpack_inplace<std::uint8_t>(invalid_size, 1),
                    psrio::ValidationError);
    CHECK_THROWS_AS(psrio::pack_inplace<std::uint8_t>(invalid_size, 3),
                    psrio::ValidationError);
    CHECK_THROWS_AS(psrio::unpack_inplace<std::uint8_t>(invalid_size, 3),
                    psrio::ValidationError);
}
