#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/unpack.hpp"
#include "psrio/formats/sigproc/header.hpp"
#include "psrio/formats/sigproc/reader.hpp"
#include "sigproc_bytes.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ios>
#include <span>
#include <string>
#include <system_error>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using psrio::formats::sigproc::ByteCount;
using psrio::formats::sigproc::FilterbankReader;
using psrio::formats::sigproc::SampleCount;
using psrio::test::SigprocBytes;

namespace {

class TempFile {
public:
    explicit TempFile(std::span<const std::byte> bytes) {
        m_path =
            std::filesystem::temp_directory_path() /
            ("psrio-fil-" + std::to_string(++m_sequence) + "-" +
             std::to_string(
                 std::chrono::steady_clock::now().time_since_epoch().count()));
        std::ofstream stream(m_path, std::ios::binary);
        REQUIRE(stream.good());
        if (!bytes.empty()) {
            stream.write(reinterpret_cast<const char*>(bytes.data()),
                         static_cast<std::streamsize>(bytes.size()));
        }
        REQUIRE(stream.good());
    }

    TempFile(const TempFile&)            = delete;
    TempFile& operator=(const TempFile&) = delete;
    TempFile(TempFile&&)                 = delete;
    TempFile& operator=(TempFile&&)      = delete;

    ~TempFile() {
        std::error_code error;
        std::filesystem::remove(m_path, error);
    }

    [[nodiscard]] const std::filesystem::path& path() const { return m_path; }

private:
    std::filesystem::path m_path;
    static inline int m_sequence{0};
};

SigprocBytes file_8bit(std::span<const std::uint8_t> samples,
                       std::int32_t nchans    = 1,
                       std::uint32_t declared = 0xFFFFFFFFU) {
    SigprocBytes bytes;
    bytes.string("HEADER_START");
    bytes.key_i32("nchans", nchans);
    bytes.key_i32("nbits", 8);
    bytes.key_f64("tsamp", 1.0);
    bytes.key_f64("fch1", 1400.0);
    bytes.key_f64("foff", -1.0);
    if (declared != 0xFFFFFFFFU) {
        bytes.key_u32("nsamples", declared);
    }
    bytes.string("HEADER_END");
    for (const std::uint8_t sample : samples) {
        bytes.u8(sample);
    }
    return bytes;
}

} // namespace

TEST_CASE("reader streams samples and stops at the end", "[sigproc][reader]") {
    const std::array<std::uint8_t, 8> samples{10, 11, 12, 13, 14, 15, 16, 17};
    const TempFile file(file_8bit(samples).out);
    FilterbankReader reader(file.path());

    std::array<std::uint8_t, 3> first{0, 0, 0};
    REQUIRE(reader.read(SampleCount{3}, std::span<std::uint8_t>{first}) == 3U);
    REQUIRE(first == std::array<std::uint8_t, 3>{10, 11, 12});
    REQUIRE(reader.tell() == 3U);

    std::array<std::uint8_t, 3> second{};
    REQUIRE(reader.read(SampleCount{3}, std::span<std::uint8_t>{second}) == 3U);
    REQUIRE(second == std::array<std::uint8_t, 3>{13, 14, 15});

    std::array<std::uint8_t, 3> tail{9, 9, 9};
    REQUIRE(reader.read(SampleCount{3}, std::span<std::uint8_t>{tail}) == 2U);
    REQUIRE(tail.at(0) == 16);
    REQUIRE(tail.at(1) == 17);
    REQUIRE(tail.at(2) == 9);
    REQUIRE(reader.read(SampleCount{3}, std::span<std::uint8_t>{tail}) == 0U);

    reader.seek(reader.header().nsamples());
    REQUIRE(reader.tell() == 8U);
    REQUIRE_THROWS_AS(reader.seek(9), psrio::ValidationError);
    reader.rewind();
    REQUIRE(reader.tell() == 0U);
}

TEST_CASE("reader rejects a destination of the wrong length",
          "[sigproc][reader]") {
    const std::array<std::uint8_t, 2> samples{1, 2};
    const TempFile file(file_8bit(samples).out);
    FilterbankReader reader(file.path());
    std::array<std::uint8_t, 1> dest{};
    REQUIRE_THROWS_AS(
        reader.read(SampleCount{2}, std::span<std::uint8_t>{dest}),
        psrio::ValidationError);
    std::array<std::uint16_t, 2> wrong_type{};
    REQUIRE_THROWS_AS(
        reader.read(SampleCount{2}, std::span<std::uint16_t>{wrong_type}),
        psrio::ValidationError);
}

TEST_CASE("raw byte reads match a zero-copy view", "[sigproc][reader]") {
    const std::array<std::uint8_t, 4> samples{1, 2, 3, 4};
    const TempFile file(file_8bit(samples, 4).out);
    FilterbankReader reader(file.path());
    REQUIRE(reader.header().bytes_per_sample() == 4U);
    REQUIRE(reader.header().nsamples() == 1U);

    const auto viewed = reader.view(ByteCount{4});
    REQUIRE(viewed.size() == 4U);
    REQUIRE(std::to_integer<unsigned char>(viewed[0]) == 1);
    REQUIRE(std::to_integer<unsigned char>(viewed[3]) == 4);
    REQUIRE(reader.tell() == 1U);

    reader.rewind();
    std::array<std::byte, 4> copied{};
    REQUIRE(reader.read(ByteCount{4}, copied) == 4U);
    REQUIRE(std::equal(viewed.begin(), viewed.end(), copied.begin()));

    reader.rewind();
    REQUIRE_THROWS_AS(
        reader.read(ByteCount{3}, std::span<std::byte>{copied}.first(3)),
        psrio::ValidationError);
    REQUIRE_THROWS_AS(reader.view(ByteCount{8}), psrio::ValidationError);
    REQUIRE(reader.tell() == 0U);
}

TEST_CASE("4-bit byte 0xAB unpacks low nibble first", "[sigproc][reader]") {
    SigprocBytes bytes;
    bytes.string("HEADER_START");
    bytes.key_i32("nchans", 2);
    bytes.key_i32("nbits", 4);
    bytes.key_f64("tsamp", 1.0);
    bytes.key_f64("fch1", 1400.0);
    bytes.key_f64("foff", -1.0);
    bytes.string("HEADER_END");
    bytes.u8(0xAB);
    const TempFile file(bytes.out);
    FilterbankReader reader(file.path());
    std::array<float, 2> samples{};
    REQUIRE(reader.read(SampleCount{1}, std::span<float>{samples}) == 1U);
    REQUIRE(samples.at(0) == 11.0F);
    REQUIRE(samples.at(1) == 10.0F);
}

TEST_CASE("signed 8-bit samples cast through int8", "[sigproc][reader]") {
    SigprocBytes bytes;
    bytes.string("HEADER_START");
    bytes.key_i32("nchans", 2);
    bytes.key_i32("nbits", 8);
    bytes.key_i8("signed", -1);
    bytes.key_f64("tsamp", 1.0);
    bytes.key_f64("fch1", 1400.0);
    bytes.key_f64("foff", -1.0);
    bytes.string("HEADER_END");
    bytes.u8(0xFF);
    bytes.u8(0x05);
    const TempFile file(bytes.out);
    FilterbankReader reader(file.path());
    REQUIRE(reader.header().samples_are_signed());
    std::array<float, 2> samples{};
    REQUIRE(reader.read(SampleCount{1}, std::span<float>{samples}) == 1U);
    REQUIRE(samples.at(0) == -1.0F);
    REQUIRE(samples.at(1) == 5.0F);
}

TEST_CASE("16-bit and 32-bit samples unpack little-endian",
          "[sigproc][reader]") {
    SigprocBytes bits16;
    bits16.string("HEADER_START");
    bits16.key_i32("nchans", 1);
    bits16.key_i32("nbits", 16);
    bits16.key_f64("tsamp", 1.0);
    bits16.key_f64("fch1", 1400.0);
    bits16.key_f64("foff", -1.0);
    bits16.string("HEADER_END");
    bits16.u8(0x02);
    bits16.u8(0x01);
    const TempFile file16(bits16.out);
    FilterbankReader reader16(file16.path());
    std::array<std::uint16_t, 1> as_int{};
    REQUIRE(reader16.read(SampleCount{1}, std::span<std::uint16_t>{as_int}) ==
            1U);
    REQUIRE(as_int.at(0) == 0x0102U);

    SigprocBytes bits32;
    bits32.string("HEADER_START");
    bits32.key_i32("nchans", 1);
    bits32.key_i32("nbits", 32);
    bits32.key_f64("tsamp", 1.0);
    bits32.key_f64("fch1", 1400.0);
    bits32.key_f64("foff", -1.0);
    bits32.string("HEADER_END");
    const float sample = 1.5F;
    const auto sample_bits =
        psrio::detail::from_little_endian(std::bit_cast<std::uint32_t>(sample));
    bits32.u8(static_cast<std::uint8_t>(sample_bits & 0xFFU));
    bits32.u8(static_cast<std::uint8_t>((sample_bits >> 8U) & 0xFFU));
    bits32.u8(static_cast<std::uint8_t>((sample_bits >> 16U) & 0xFFU));
    bits32.u8(static_cast<std::uint8_t>((sample_bits >> 24U) & 0xFFU));
    const TempFile file32(bits32.out);
    FilterbankReader reader32(file32.path());
    std::array<float, 1> as_float{};
    REQUIRE(reader32.read(SampleCount{1}, std::span<float>{as_float}) == 1U);
    REQUIRE(as_float.at(0) == 1.5F);
    reader32.rewind();
    std::array<std::uint8_t, 1> as_byte{};
    REQUIRE_THROWS_AS(
        reader32.read(SampleCount{1}, std::span<std::uint8_t>{as_byte}),
        psrio::ValidationError);
}

TEST_CASE("multi-if samples stay in file order", "[sigproc][reader]") {
    SigprocBytes bytes;
    bytes.string("HEADER_START");
    bytes.key_i32("nchans", 2);
    bytes.key_i32("nifs", 2);
    bytes.key_i32("nbits", 8);
    bytes.key_f64("tsamp", 1.0);
    bytes.key_f64("fch1", 1400.0);
    bytes.key_f64("foff", -1.0);
    bytes.string("HEADER_END");
    for (std::uint8_t value = 1; value <= 8; ++value) {
        bytes.u8(value);
    }
    const TempFile file(bytes.out);
    FilterbankReader reader(file.path());
    REQUIRE(reader.header().bytes_per_sample() == 4U);
    REQUIRE(reader.header().nsamples() == 2U);
    std::array<float, 4> first{};
    REQUIRE(reader.read(SampleCount{1}, std::span<float>{first}) == 1U);
    REQUIRE(first == std::array<float, 4>{1.0F, 2.0F, 3.0F, 4.0F});
}

TEST_CASE("a smaller declared nsamples limits the reader",
          "[sigproc][reader]") {
    const std::array<std::uint8_t, 10> samples{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    const TempFile file(file_8bit(samples, 1, 3).out);
    FilterbankReader reader(file.path());
    REQUIRE(reader.header().nsamples() == 3U);
    REQUIRE(reader.header().warnings().size() == 1U);
    REQUIRE_THAT(reader.header().warnings().front(),
                 ContainsSubstring("disagrees"));

    std::array<std::uint8_t, 10> dest{};
    dest.fill(0);
    REQUIRE(reader.read(SampleCount{10}, std::span<std::uint8_t>{dest}) == 3U);
    REQUIRE(dest.at(0) == 1);
    REQUIRE(dest.at(1) == 2);
    REQUIRE(dest.at(2) == 3);
    REQUIRE(dest.at(3) == 0);
    REQUIRE_THROWS_AS(reader.seek(4), psrio::ValidationError);
}

TEST_CASE("missing and empty files fail with the matching error",
          "[sigproc][reader]") {
    REQUIRE_THROWS_AS(FilterbankReader("/no/such/psrio-filterbank.fil"),
                      psrio::IoError);
    const TempFile empty({});
    REQUIRE_THROWS_WITH(FilterbankReader(empty.path()),
                        ContainsSubstring("empty"));
}

TEST_CASE("a time-series data_type is rejected by the filterbank reader",
          "[sigproc][reader]") {
    SigprocBytes bytes;
    bytes.string("HEADER_START");
    bytes.key_i32("data_type", 2);
    bytes.key_i32("nchans", 1);
    bytes.key_i32("nbits", 32);
    bytes.key_f64("tsamp", 1.0);
    bytes.key_f64("fch1", 1400.0);
    bytes.key_f64("foff", -1.0);
    bytes.string("HEADER_END");
    const TempFile file(bytes.out);
    REQUIRE_THROWS_AS(FilterbankReader(file.path()), psrio::ValidationError);
}

TEST_CASE("reader respects BitOrder in constructor, setter, and read overload",
          "[sigproc][reader]") {
    // 2 channels, 4-bit, 1 sample. 2 chans * 4-bit = 8 bits = 1 byte.
    // Byte value 0xAB: high nibble 10, low nibble 11.
    SigprocBytes bytes;
    bytes.string("HEADER_START");
    bytes.key_i32("nchans", 2);
    bytes.key_i32("nbits", 4);
    bytes.key_f64("tsamp", 1.0);
    bytes.key_f64("fch1", 1400.0);
    bytes.key_f64("foff", -1.0);
    bytes.string("HEADER_END");
    bytes.u8(0xAB);
    const TempFile file(bytes.out);

    // Default reader is LSB first: {11, 10}
    FilterbankReader reader(file.path());
    REQUIRE(reader.bit_order() == psrio::BitOrder::kLsbFirst);
    std::array<std::uint8_t, 2> lsb_out{};
    REQUIRE(reader.read(SampleCount{1}, std::span<std::uint8_t>{lsb_out}) ==
            1U);
    REQUIRE(lsb_out == std::array<std::uint8_t, 2>{11, 10});

    // Rewind and use set_bit_order(kMsbFirst): {10, 11}
    reader.rewind();
    reader.set_bit_order(psrio::BitOrder::kMsbFirst);
    REQUIRE(reader.bit_order() == psrio::BitOrder::kMsbFirst);
    std::array<std::uint8_t, 2> msb_out{};
    REQUIRE(reader.read(SampleCount{1}, std::span<std::uint8_t>{msb_out}) ==
            1U);
    REQUIRE(msb_out == std::array<std::uint8_t, 2>{10, 11});

    // Rewind and use per-read overload override
    reader.rewind();
    std::array<float, 2> lsb_floats{};
    REQUIRE(reader.read(SampleCount{1}, std::span<float>{lsb_floats},
                        psrio::BitOrder::kLsbFirst) == 1U);
    REQUIRE(lsb_floats == std::array<float, 2>{11.0F, 10.0F});
}

TEST_CASE("reader convenience allocating methods read_samples and read_bytes",
          "[sigproc][reader]") {
    const std::array<std::uint8_t, 6> samples{10, 20, 30, 40, 50, 60};
    const TempFile file(
        file_8bit(samples, 2).out); // 2 chans, 8-bit => 3 samples
    FilterbankReader reader(file.path());
    REQUIRE(reader.header().nsamples() == 3U);

    // Convenience read_samples<float>
    const auto floats = reader.read_samples<float>(SampleCount{2});
    REQUIRE(floats.size() == 4U); // 2 samples * 2 chans
    REQUIRE(floats == std::vector<float>{10.0F, 20.0F, 30.0F, 40.0F});
    REQUIRE(reader.tell() == 2U);

    // Convenience read_bytes
    reader.rewind();
    const auto raw_bytes = reader.read_bytes(ByteCount{4});
    REQUIRE(raw_bytes.size() == 4U);
    REQUIRE(std::to_integer<unsigned char>(raw_bytes[0]) == 10);
    REQUIRE(std::to_integer<unsigned char>(raw_bytes[3]) == 40);
    REQUIRE(reader.tell() == 2U);
}
