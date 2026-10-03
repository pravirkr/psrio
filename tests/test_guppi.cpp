#include "psrio/block_source.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/formats/guppi.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <chrono>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iosfwd>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

using psrio::formats::guppi::GuppiHeader;
using psrio::formats::guppi::RawReader;

namespace {

class TempFile {
public:
    explicit TempFile(std::span<const std::byte> bytes) {
        static int sequence = 0;
        m_path =
            std::filesystem::temp_directory_path() /
            ("psrio-guppi-" + std::to_string(++sequence) + "-" +
             std::to_string(
                 std::chrono::steady_clock::now().time_since_epoch().count()));
        std::ofstream stream(m_path, std::ios::binary);
        REQUIRE(stream.good());
        stream.write(reinterpret_cast<const char*>(bytes.data()),
                     static_cast<std::streamsize>(bytes.size()));
        REQUIRE(stream.good());
    }

    ~TempFile() noexcept {
        std::error_code ec;
        std::filesystem::remove(m_path, ec);
    }

    TempFile(const TempFile&)            = delete;
    TempFile& operator=(const TempFile&) = delete;
    TempFile(TempFile&&)                 = delete;
    TempFile& operator=(TempFile&&)      = delete;

    [[nodiscard]] const std::filesystem::path& path() const { return m_path; }

private:
    std::filesystem::path m_path;
};

std::string make_record(std::string_view text) {
    std::string record{text};
    if (record.size() < psrio::formats::guppi::kRecordBytes) {
        record.resize(psrio::formats::guppi::kRecordBytes, ' ');
    }
    return record;
}

void append_text(std::vector<std::byte>& out, std::string_view text) {
    const auto offset = out.size();
    out.resize(offset + text.size());
    std::memcpy(out.data() + offset, text.data(), text.size());
}

void append_record(std::vector<std::byte>& out, std::string_view text) {
    append_text(out, make_record(text));
}

std::vector<std::byte> standard_header_bytes() {
    std::vector<std::byte> bytes;
    append_record(bytes, "obsnchan= 4");
    append_record(bytes, "npol    = 1");
    append_record(bytes, "nbits   = 8");
    append_record(bytes, "blocsize= 128");
    append_record(bytes, "obsfreq = 1400.0");
    append_record(bytes, "chan_bw = 1.0");
    append_record(bytes, "END");
    return bytes;
}

} // namespace

static_assert(!psrio::concepts::BlockReader<RawReader>);
static_assert(!std::is_constructible_v<psrio::BlockSource, RawReader>);

TEST_CASE("GuppiHeader parses records with case-insensitive keys", "[guppi]") {
    const auto bytes = standard_header_bytes();
    GuppiHeader header;
    REQUIRE(header.parse(bytes) == bytes.size());

    CHECK(header.has_key("OBSNCHAN"));
    CHECK(header.has_key("obsnchan"));
    CHECK(header.get<std::int64_t>("OBSNCHAN") == 4);
    CHECK(header.get<std::int64_t>("ObsNChan") == 4);
    CHECK(header.block_type() == "Array<Complex<Int8>, 3>");
}

TEST_CASE("GuppiHeader derives ntime, blocksize, and channel frequency",
          "[guppi]") {
    const auto bytes = standard_header_bytes();
    GuppiHeader header;
    header.parse(bytes);

    // blocsize=128, obsnchan=4, npol=1, nbits=8 ->
    // ntime = 8*128 / (2*4*1*8) = 16
    CHECK(header.ntime() == 16);
    CHECK(header.blocksize() == 64);
    CHECK(header.antnchan() == 4);
    // channel 1 = 1400 - 4*1/2 + 1*0.5 = 1398.5
    CHECK(header.channel_frequency(1) == 1398.5);

    const auto freqs = header.channel_frequencies(1);
    REQUIRE(freqs.size() == 4U);
    CHECK(freqs.front() == 1398.5);
    CHECK(freqs.back() == 1401.5);
}

TEST_CASE("GuppiHeader::parse throws on a truncated header", "[guppi]") {
    const auto record = make_record("obsnchan= 4");
    std::vector<std::byte> bytes(record.size());
    std::memcpy(bytes.data(), record.data(), record.size());
    GuppiHeader header;
    CHECK_THROWS_AS(header.parse(bytes), psrio::FormatError);
}

TEST_CASE("GuppiHeader stores quoted strings and preserves raw NPOL",
          "[guppi]") {
    std::vector<std::byte> bytes;
    append_record(bytes, "SRC_NAME = 'B1937+21'");
    append_record(bytes, "NPOL     = 4");
    append_record(bytes, "END");
    GuppiHeader header;
    header.parse(bytes);
    CHECK(header.get<std::string>("src_name") == "B1937+21");
    CHECK(header.get<std::int64_t>("NPOL") == 4);
    CHECK_THROWS_AS(header.get<std::int64_t>("OBSNCHAN"), psrio::FormatError);
    CHECK_THROWS_AS(header.channel_frequencies(3, 1), psrio::ValidationError);
}

TEST_CASE("GuppiHeader strips comments from records", "[guppi]") {
    std::vector<std::byte> bytes;
    append_record(bytes, "SRC_NAME= 'B1937+21' / Target pulsar");
    append_record(bytes, "TELESCOP= 'GBT' / Green Bank Telescope");
    append_record(bytes, "OBSFREQ = 1400.0 / Centre frequency in MHz");
    append_record(bytes, "OBSNCHAN= 64 / Channel count");
    append_record(bytes, "NPOL    = 2 / Complex polarisations");
    append_record(bytes, "NBITS   = 8 / Bits per component");
    append_record(bytes, "END");
    GuppiHeader header;
    header.parse(bytes);
    CHECK(header.get<std::string>("SRC_NAME") == "B1937+21");
    CHECK(header.get<std::string>("TELESCOP") == "GBT");
    CHECK(header.get<double>("OBSFREQ") == 1400.0);
    CHECK(header.get<std::int64_t>("OBSNCHAN") == 64);
    CHECK(header.get<std::int64_t>("NPOL") == 2);
    CHECK(header.get<std::int64_t>("NBITS") == 8);
}

TEST_CASE("DIRECTIO padding reaches the next 512-byte boundary", "[guppi]") {
    std::vector<std::byte> unaligned;
    append_record(unaligned, "DIRECTIO = 1");
    append_record(unaligned, "OBSNCHAN= 4");
    append_record(unaligned, "NPOL    = 1");
    append_record(unaligned, "NBITS   = 8");
    append_record(unaligned, "BLOCSIZE= 16");
    append_record(unaligned, "OBSFREQ = 1400.0");
    append_record(unaligned, "CHAN_BW = 1.0");
    append_record(unaligned, "END");
    REQUIRE(unaligned.size() == 640U);
    unaligned.resize(1024, std::byte{0xA5});
    GuppiHeader header;
    CHECK(header.parse(unaligned) == 1024U);

    std::vector<std::byte> no_pad;
    append_record(no_pad, "DIRECTIO = 0");
    append_record(no_pad, "OBSNCHAN= 4");
    append_record(no_pad, "NPOL    = 1");
    append_record(no_pad, "NBITS   = 8");
    append_record(no_pad, "BLOCSIZE= 16");
    append_record(no_pad, "OBSFREQ = 1400.0");
    append_record(no_pad, "CHAN_BW = 1.0");
    append_record(no_pad, "END");
    GuppiHeader plain;
    CHECK(plain.parse(no_pad) == no_pad.size());

    std::vector<std::byte> aligned;
    append_record(aligned, "DIRECTIO = 1");
    append_record(aligned, "OBSNCHAN= 1");
    append_record(aligned, "NPOL    = 1");
    append_record(aligned, "NBITS   = 8");
    append_record(aligned, "BLOCSIZE= 8");
    append_record(aligned, "OBSFREQ = 100.0");
    append_record(aligned, "CHAN_BW = 1.0");
    for (int extra = 0; extra < 24; ++extra) {
        append_record(aligned, "COMMENT = pad");
    }
    append_record(aligned, "END");
    REQUIRE(aligned.size() == 2560U);
    REQUIRE(aligned.size() % 512U == 0U);
    GuppiHeader aligned_header;
    CHECK(aligned_header.parse(aligned) == 2560U);
}

TEST_CASE("unpack_complex reads little-endian real and imaginary pairs",
          "[guppi]") {
    const std::array<std::byte, 4> packed_i8{
        std::byte{1},
        std::byte{2},
        std::byte{3},
        std::byte{4},
    };
    std::array<std::complex<std::int8_t>, 2> out8{};
    psrio::formats::guppi::unpack_complex<std::int8_t>(packed_i8, out8);
    CHECK(out8[0] == std::complex<std::int8_t>(1, 2));
    CHECK(out8[1] == std::complex<std::int8_t>(3, 4));

    const std::array<std::byte, 4> packed_i16{
        std::byte{0x02},
        std::byte{0x01},
        std::byte{0x04},
        std::byte{0x03},
    };
    std::array<std::complex<std::int16_t>, 1> out16{};
    psrio::formats::guppi::unpack_complex<std::int16_t>(packed_i16, out16);
    CHECK(out16[0] == std::complex<std::int16_t>(0x0102, 0x0304));

    std::array<std::complex<std::int8_t>, 1> mismatch{};
    const std::span<const std::byte> packed_span{packed_i8};
    const std::span<std::complex<std::int8_t>> mismatch_span{mismatch};
    CHECK_THROWS_AS(psrio::formats::guppi::unpack_complex<std::int8_t>(
                        packed_span, mismatch_span),
                    psrio::ValidationError);
}

TEST_CASE("RawReader walks headers and payload bytes", "[guppi]") {
    std::vector<std::byte> bytes;
    append_record(bytes, "obsnchan= 4");
    append_record(bytes, "npol    = 1");
    append_record(bytes, "nbits   = 8");
    append_record(bytes, "blocsize= 16");
    append_record(bytes, "obsfreq = 1400.0");
    append_record(bytes, "chan_bw = 1.0");
    append_record(bytes, "END");
    const auto first_header = bytes.size();
    const std::array<std::uint8_t, 16> payload{
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    };
    for (const auto value : payload) {
        bytes.push_back(std::byte{value});
    }
    append_record(bytes, "SRC_NAME = 'NEXT'");
    append_record(bytes, "obsnchan= 4");
    append_record(bytes, "npol    = 1");
    append_record(bytes, "nbits   = 8");
    append_record(bytes, "blocsize= 16");
    append_record(bytes, "obsfreq = 1400.0");
    append_record(bytes, "chan_bw = 1.0");
    append_record(bytes, "END");
    for (const auto value : payload) {
        bytes.push_back(std::byte{value});
    }

    const TempFile file(bytes);
    RawReader reader(file.path());
    const auto header = reader.read_header();
    REQUIRE(header.ntime() == 2);
    REQUIRE(header.blocksize() == 8);
    REQUIRE(reader.tell() == first_header);

    const auto view = reader.view_bytes(4);
    REQUIRE(view.size() == 4U);
    CHECK(view[0] == std::byte{1});
    CHECK(view[3] == std::byte{4});
    REQUIRE(reader.tell() == first_header + 4U);

    reader.seek(first_header + 2U);
    std::vector<std::byte> copied(4);
    REQUIRE(reader.read_bytes(4, copied) == 4U);
    CHECK(copied[0] == std::byte{3});
    CHECK(copied[3] == std::byte{6});

    auto rest = reader.read_bytes(16);
    REQUIRE(rest.size() == 10U);
    CHECK(rest.front() == std::byte{7});
    CHECK(rest.back() == std::byte{16});
    REQUIRE(reader.tell() == first_header + 16U);
    CHECK_THROWS_AS(reader.read_bytes(1), psrio::ValidationError);

    const auto second = reader.read_header();
    CHECK(second.get<std::string>("SRC_NAME") == "NEXT");
    const auto second_payload = reader.tell();
    auto again                = reader.read_bytes(16);
    REQUIRE(again.size() == 16U);

    reader.seek(second_payload);
    std::array<std::complex<std::int8_t>, 8> samples{};
    const auto packed = reader.view_bytes(16);
    psrio::formats::guppi::unpack_complex<std::int8_t>(packed, samples);
    CHECK(samples[0] == std::complex<std::int8_t>(1, 2));
    CHECK(samples[7] == std::complex<std::int8_t>(15, 16));

    reader.rewind();
    REQUIRE(reader.tell() == 0U);
    const auto rewound_header = reader.read_header();
    REQUIRE(rewound_header.get<std::int64_t>("OBSNCHAN") == 4);

    reader.seek(reader.tell() + 1U);
    CHECK_THROWS_AS(reader.read_header(), psrio::ValidationError);
    REQUIRE_THROWS_AS(reader.seek(bytes.size() + 1U), psrio::ValidationError);
}

TEST_CASE("RawReader reports a short file and a missing path", "[guppi]") {
    auto bytes = standard_header_bytes();
    // BLOCSIZE is 128, but no payload bytes follow the header.
    const TempFile file(bytes);
    RawReader reader(file.path());
    CHECK_THROWS_AS(reader.read_header(), psrio::FormatError);
    CHECK(reader.tell() == 0U);
    CHECK_THROWS_AS(RawReader{"/no/such/psrio-guppi.raw"}, psrio::IoError);
}
