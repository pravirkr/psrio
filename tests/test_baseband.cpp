#include "psrio/baseband.hpp"
#include "psrio/block_source.hpp"
#include "psrio/formats/dada.hpp"
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
#include <utility>
#include <vector>

using psrio::BasebandSource;
using psrio::FrequencyStitch;
using psrio::MemoryBaseband;
using psrio::formats::dada::DadaReader;
using psrio::formats::guppi::GuppiReader;
using psrio::formats::guppi::GuppiSet;

namespace {

class TempFile {
public:
    explicit TempFile(std::span<const std::byte> bytes) {
        static int sequence = 0;
        m_path =
            std::filesystem::temp_directory_path() /
            ("psrio-baseband-" + std::to_string(++sequence) + "-" +
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

void append_text(std::vector<std::byte>& out, std::string_view text) {
    const auto offset = out.size();
    out.resize(offset + text.size());
    std::memcpy(out.data() + offset, text.data(), text.size());
}

void append_record(std::vector<std::byte>& out, std::string_view text) {
    std::string record{text};
    if (record.size() < psrio::formats::guppi::kRecordBytes) {
        record.resize(psrio::formats::guppi::kRecordBytes, ' ');
    }
    append_text(out, record);
}

std::vector<std::byte> guppi_block(std::span<const std::string_view> cards,
                                   std::span<const std::byte> payload) {
    std::vector<std::byte> out;
    for (const auto card : cards) {
        append_record(out, card);
    }
    append_record(out, "END");
    out.insert(out.end(), payload.begin(), payload.end());
    return out;
}

std::vector<std::byte> make_dada(std::span<const std::string_view> lines,
                                 std::span<const std::byte> payload,
                                 std::uint64_t header_bytes = 256) {
    std::string text = "HEADER DADA\nHDR_SIZE " + std::to_string(header_bytes) + "\n";
    for (const auto line : lines) {
        text.append(line);
        text.push_back('\n');
    }
    text += "# end of header\n";
    REQUIRE(text.size() <= header_bytes);
    std::vector<std::byte> out(static_cast<std::size_t>(header_bytes + payload.size()));
    std::memcpy(out.data(), text.data(), text.size());
    if (!payload.empty()) {
        std::memcpy(out.data() + header_bytes, payload.data(), payload.size());
    }
    return out;
}

template <std::size_t N>
std::vector<std::byte> guppi_block(const std::array<std::string_view, N>& cards,
                                   std::span<const std::byte> payload) {
    return guppi_block(std::span<const std::string_view>(cards.data(), cards.size()),
                       payload);
}

template <std::size_t N>
std::vector<std::byte> make_dada(const std::array<std::string_view, N>& lines,
                                 std::span<const std::byte> payload,
                                 std::uint64_t header_bytes = 256) {
    return make_dada(std::span<const std::string_view>(lines.data(), lines.size()),
                     payload, header_bytes);
}

psrio::BasebandHeader complex_header(std::uint64_t nchan, int nbit) {
    psrio::BasebandHeader header;
    header.npol           = 1;
    header.nchan          = nchan;
    header.nants          = 1;
    header.nbit           = nbit;
    header.ndim           = 2;
    header.samples_signed = true;
    header.tsamp          = 1.0e-6;
    header.foff           = 1.0;
    header.fch1           = 1000.0;
    return header;
}

} // namespace

static_assert(psrio::concepts::BasebandReader<MemoryBaseband>);
static_assert(psrio::concepts::BasebandReader<FrequencyStitch>);
static_assert(psrio::concepts::BasebandReader<GuppiReader>);
static_assert(psrio::concepts::BasebandReader<GuppiSet>);
static_assert(psrio::concepts::BasebandReader<DadaReader>);
static_assert(psrio::concepts::BasebandReader<BasebandSource>);
static_assert(!psrio::concepts::BlockReader<GuppiReader>);
static_assert(!psrio::concepts::BlockReader<DadaReader>);
static_assert(!std::is_constructible_v<psrio::BlockSource, GuppiReader>);

TEST_CASE("MemoryBaseband unpacks canonical complex bytes", "[baseband]") {
    const std::vector<std::byte> packed{std::byte{1}, std::byte{2}, std::byte{3},
                                        std::byte{4},};
    MemoryBaseband memory(complex_header(2, 8), packed);
    BasebandSource source(std::move(memory));
    CHECK(source.nchan() == 2U);
    CHECK(source.ndim() == 2);
    CHECK(source.bytes_per_sample() == 4U);
    auto samples = source.read_samples<std::complex<std::int8_t>>(1);
    REQUIRE(samples.size() == 2U);
    CHECK(samples[0] == std::complex<std::int8_t>(1, 2));
    CHECK(samples[1] == std::complex<std::int8_t>(3, 4));
    CHECK(source.read_block(1).empty());
}

TEST_CASE("FrequencyStitch orders bands by lowest channel frequency", "[baseband]") {
    psrio::BasebandHeader low = complex_header(1, 8);
    low.fch1                   = 1000.0;
    psrio::BasebandHeader high = complex_header(1, 8);
    high.fch1                  = 2000.0;
    std::vector<psrio::BasebandSource> bands;
    bands.emplace_back(MemoryBaseband(
        high, std::vector<std::byte>{std::byte{3}, std::byte{4}}));
    bands.emplace_back(MemoryBaseband(
        low, std::vector<std::byte>{std::byte{1}, std::byte{2}}));
    FrequencyStitch stitch(std::move(bands));
    CHECK(stitch.nchan() == 2U);
    CHECK(stitch.fch1() == 1000.0);
    auto samples = stitch.read_samples<std::complex<std::int8_t>>(1);
    REQUIRE(samples.size() == 2U);
    CHECK(samples[0] == std::complex<std::int8_t>(1, 2));
    CHECK(samples[1] == std::complex<std::int8_t>(3, 4));
}

TEST_CASE("GuppiReader transposes channel-major blocks", "[baseband][guppi]") {
    const std::vector<std::byte> payload{std::byte{1}, std::byte{2}, std::byte{3},
                                         std::byte{4}, std::byte{5}, std::byte{6},
                                         std::byte{7}, std::byte{8}};
    const std::string_view cards[] = {
        "OBSNCHAN= 2", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 8",
        "OBSFREQ = 1400.0", "CHAN_BW = 1.0", "TBIN    = 1.0e-6", "PKTIDX  = 0"};
    const TempFile file(guppi_block(cards, payload));
    GuppiReader reader(file.path());
    CHECK(reader.nsamples() == 2U);
    CHECK(std::vector<std::byte>(reader.view_native(0).begin(),
                                 reader.view_native(0).end()) == payload);
    CHECK(reader.header().order == psrio::BasebandOrder::kChannelMajor);
    std::vector<std::byte> block(8);
    REQUIRE(reader.read_block(2, block) == 2U);
    const std::vector<std::byte> expected{std::byte{1}, std::byte{2}, std::byte{5},
                                          std::byte{6}, std::byte{3}, std::byte{4},
                                          std::byte{7}, std::byte{8}};
    CHECK(block == expected);
    reader.seek(0);
    auto samples = reader.read_samples<std::complex<std::int8_t>>(2);
    REQUIRE(samples.size() == 4U);
    CHECK(samples[0] == std::complex<std::int8_t>(1, 2));
    CHECK(samples[1] == std::complex<std::int8_t>(5, 6));
    CHECK(samples[2] == std::complex<std::int8_t>(3, 4));
    CHECK(samples[3] == std::complex<std::int8_t>(7, 8));
    CHECK_THROWS_AS(reader.read_samples<float>(1), psrio::ValidationError);
}

TEST_CASE("GuppiReader drops overlap and continues at the next block", "[baseband][guppi]") {
    const std::vector<std::byte> first{std::byte{1}, std::byte{2}, std::byte{3},
                                       std::byte{4},};
    const std::vector<std::byte> second{std::byte{5}, std::byte{6}, std::byte{7},
                                        std::byte{8},};
    const std::array<std::string_view, 8> head = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 4",
        "OVERLAP = 1", "OBSFREQ = 1400.0", "CHAN_BW = 1.0", "PKTIDX  = 0",};
    const std::array<std::string_view, 8> next = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 4",
        "OVERLAP = 1", "OBSFREQ = 1400.0", "CHAN_BW = 1.0", "PKTIDX  = 1",};
    auto bytes = guppi_block(head, first);
    const auto more = guppi_block(next, second);
    bytes.insert(bytes.end(), more.begin(), more.end());
    const TempFile file(bytes);
    GuppiReader reader(file.path());
    CHECK(reader.nsamples() == 2U);
    auto samples = reader.read_samples<std::complex<std::int8_t>>(2);
    REQUIRE(samples.size() == 2U);
    CHECK(samples[0] == std::complex<std::int8_t>(1, 2));
    CHECK(samples[1] == std::complex<std::int8_t>(5, 6));
}

TEST_CASE("GuppiReader unpacks 16-bit, 4-bit, and 2-bit voltages", "[baseband][guppi]") {
    const std::vector<std::byte> wide{std::byte{0x02}, std::byte{0x01},
                                      std::byte{0x04}, std::byte{0x03},};
    const std::array<std::string_view, 7> wide_cards = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 16", "BLOCSIZE= 4",
        "OBSFREQ = 100.0", "CHAN_BW = 1.0", "PKTIDX  = 0",};
    const TempFile wide_file(guppi_block(wide_cards, wide));
    auto wide_samples =
        GuppiReader(wide_file.path()).read_samples<std::complex<std::int16_t>>(1);
    REQUIRE(wide_samples.size() == 1U);
    CHECK(wide_samples[0] == std::complex<std::int16_t>(0x0102, 0x0304));

    const std::vector<std::byte> nibbles{std::byte{0x21},};
    const std::array<std::string_view, 7> nibble_cards = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 4", "BLOCSIZE= 1",
        "OBSFREQ = 100.0", "CHAN_BW = 1.0", "PKTIDX  = 0",};
    const TempFile nibble_file(guppi_block(nibble_cards, nibbles));
    auto nibble_samples =
        GuppiReader(nibble_file.path()).read_samples<std::complex<std::int8_t>>(1);
    REQUIRE(nibble_samples.size() == 1U);
    CHECK(nibble_samples[0] == std::complex<std::int8_t>(2, 1));

    const std::vector<std::byte> packed{std::byte{0x14}, std::byte{0xBE},};
    const std::array<std::string_view, 7> packed_cards = {
        "OBSNCHAN= 2", "NPOL    = 1", "NBITS   = 2", "BLOCSIZE= 2",
        "OBSFREQ = 100.0", "CHAN_BW = 1.0", "PKTIDX  = 0",};
    const TempFile packed_file(guppi_block(packed_cards, packed));
    auto packed_samples =
        GuppiReader(packed_file.path()).read_samples<std::complex<std::int8_t>>(2);
    REQUIRE(packed_samples.size() == 4U);
    CHECK(packed_samples[0] == std::complex<std::int8_t>(0, 1));
    CHECK(packed_samples[1] == std::complex<std::int8_t>(-2, -1));
    CHECK(packed_samples[2] == std::complex<std::int8_t>(1, 0));
    CHECK(packed_samples[3] == std::complex<std::int8_t>(-1, -2));
}

TEST_CASE("GuppiReader applies VDIF offset-binary bytes", "[baseband][guppi]") {
    const std::vector<std::byte> payload{std::byte{0x00}, std::byte{0x01},};
    const std::array<std::string_view, 8> cards = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 2",
        "OBSFREQ = 100.0", "CHAN_BW = 1.0", "PKTFMT  = 'VDIF'", "PKTIDX  = 0",};
    const TempFile file(guppi_block(cards, payload));
    auto samples =
        GuppiReader(file.path()).read_samples<std::complex<std::int8_t>>(1);
    REQUIRE(samples.size() == 1U);
    CHECK(samples[0] == std::complex<std::int8_t>(0, static_cast<std::int8_t>(0x01U ^ 0x80U)));
}

TEST_CASE("NPOL 4 and NPOL 2 describe the same complex block", "[baseband][guppi]") {
    const std::vector<std::byte> payload{std::byte{1}, std::byte{2}, std::byte{3},
                                         std::byte{4},};
    const std::array<std::string_view, 7> two = {
        "OBSNCHAN= 1", "NPOL    = 2", "NBITS   = 8", "BLOCSIZE= 4",
        "OBSFREQ = 100.0", "CHAN_BW = 1.0", "PKTIDX  = 0",};
    const std::array<std::string_view, 7> four = {
        "OBSNCHAN= 1", "NPOL    = 4", "NBITS   = 8", "BLOCSIZE= 4",
        "OBSFREQ = 100.0", "CHAN_BW = 1.0", "PKTIDX  = 0",};
    const TempFile two_file(guppi_block(two, payload));
    const TempFile four_file(guppi_block(four, payload));
    CHECK(GuppiReader(two_file.path()).nsamples() ==
          GuppiReader(four_file.path()).nsamples());
    CHECK(GuppiReader(four_file.path()).npol() == 2U);
}

TEST_CASE("GuppiReader joins time-split files and handles PKTIDX jumps with zero-fill",
          "[baseband][guppi]") {
    const std::vector<std::byte> first{std::byte{1}, std::byte{2}};
    const std::vector<std::byte> second{std::byte{3}, std::byte{4}};
    const std::array<std::string_view, 7> head = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 2",
        "OBSFREQ = 1400.0", "CHAN_BW = 1.0", "PKTIDX  = 0",};
    const std::array<std::string_view, 7> next = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 2",
        "OBSFREQ = 1400.0", "CHAN_BW = 1.0", "PKTIDX  = 1",};
    const TempFile file0(guppi_block(head, first));
    const TempFile file1(guppi_block(next, second));
    const std::array<std::filesystem::path, 2> paths = {file0.path(), file1.path()};
    GuppiReader reader(paths);
    auto samples = reader.read_samples<std::complex<std::int8_t>>(2);
    REQUIRE(samples.size() == 2U);
    CHECK(samples[0] == std::complex<std::int8_t>(1, 2));
    CHECK(samples[1] == std::complex<std::int8_t>(3, 4));

    const std::array<std::string_view, 7> jumped = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 2",
        "OBSFREQ = 1400.0", "CHAN_BW = 1.0", "PKTIDX  = 5",};
    const TempFile gap(guppi_block(jumped, second));
    const std::array<std::filesystem::path, 2> broken = {file0.path(), gap.path()};
    GuppiReader gap_reader(broken);
    CHECK(gap_reader.dropped_packets() == 4U);
    CHECK(gap_reader.dropped_samples() == 4U);
    CHECK(gap_reader.nsamples() == 6U);
    auto all_samps = gap_reader.read_samples<std::complex<std::int8_t>>(6);
    REQUIRE(all_samps.size() == 6U);
    CHECK(all_samps[0] == std::complex<std::int8_t>(1, 2));
    CHECK(all_samps[1] == std::complex<std::int8_t>(0, 0));
    CHECK(all_samps[2] == std::complex<std::int8_t>(0, 0));
    CHECK(all_samps[3] == std::complex<std::int8_t>(0, 0));
    CHECK(all_samps[4] == std::complex<std::int8_t>(0, 0));
    CHECK(all_samps[5] == std::complex<std::int8_t>(3, 4));

    // Non-monotonic backwards PKTIDX should throw FormatError
    const std::array<std::string_view, 7> backwards = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 2",
        "OBSFREQ = 1400.0", "CHAN_BW = 1.0", "PKTIDX  = 0",};
    const TempFile back(guppi_block(backwards, second));
    const std::array<std::filesystem::path, 2> bad_order = {file1.path(), back.path()};
    CHECK_THROWS_AS(GuppiReader(bad_order), psrio::FormatError);
}

TEST_CASE("GuppiReader stitches subbands in frequency order", "[baseband][guppi]") {
    const std::vector<std::byte> low_payload{std::byte{1}, std::byte{2}};
    const std::vector<std::byte> high_payload{std::byte{3}, std::byte{4}};
    const std::array<std::string_view, 8> low_cards = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 2",
        "OBSFREQ = 1000.0", "CHAN_BW = 1.0", "TBIN    = 1.0e-6", "PKTIDX  = 0",};
    const std::array<std::string_view, 8> high_cards = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 2",
        "OBSFREQ = 2000.0", "CHAN_BW = 1.0", "TBIN    = 1.0e-6", "PKTIDX  = 0",};
    const TempFile high(guppi_block(high_cards, high_payload));
    const TempFile low(guppi_block(low_cards, low_payload));
    const std::array<std::filesystem::path, 2> paths = {high.path(), low.path()};
    auto opened = GuppiReader::open(paths);
    CHECK(opened.nchan() == 2U);
    CHECK(opened.fch1() == 1000.0);
    CHECK(opened.foff() == 1.0);
    auto samples = opened.read_samples<std::complex<std::int8_t>>(1);
    REQUIRE(samples.size() == 2U);
    CHECK(samples[0] == std::complex<std::int8_t>(1, 2));
    CHECK(samples[1] == std::complex<std::int8_t>(3, 4));

    const std::vector<std::byte> low_band{std::byte{11}, std::byte{1}, std::byte{22},
                                          std::byte{2},};
    const std::vector<std::byte> high_band{std::byte{33}, std::byte{3}, std::byte{44},
                                           std::byte{4},};
    const std::array<std::string_view, 8> low_neg = {
        "OBSNCHAN= 2", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 4",
        "OBSFREQ = 1000.0", "CHAN_BW = -1.0", "TBIN    = 1.0e-6", "PKTIDX  = 0",};
    const std::array<std::string_view, 8> high_neg = {
        "OBSNCHAN= 2", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 4",
        "OBSFREQ = 1004.0", "CHAN_BW = -1.0", "TBIN    = 1.0e-6", "PKTIDX  = 0",};
    const TempFile low_file(guppi_block(low_neg, low_band));
    const TempFile high_file(guppi_block(high_neg, high_band));
    const std::array<std::filesystem::path, 2> negative = {high_file.path(), low_file.path()};
    auto reversed = GuppiReader::open(negative);
    CHECK(reversed.foff() == -1.0);
    auto native = reversed.read_samples<std::complex<std::int8_t>>(1);
    REQUIRE(native.size() == 4U);
    CHECK(native[0] == std::complex<std::int8_t>(33, 3));
    CHECK(native[1] == std::complex<std::int8_t>(44, 4));
    CHECK(native[2] == std::complex<std::int8_t>(11, 1));
    CHECK(native[3] == std::complex<std::int8_t>(22, 2));

    const std::array<std::string_view, 8> mismatched = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 2",
        "OBSFREQ = 3000.0", "CHAN_BW = 1.0", "TBIN    = 2.0e-6", "PKTIDX  = 0",};
    const TempFile other(guppi_block(mismatched, low_payload));
    const std::array<std::filesystem::path, 2> bad = {low.path(), other.path()};
    CHECK_THROWS_AS(GuppiReader::open(bad), psrio::ValidationError);
}

TEST_CASE("DadaReader parses a text header and unpacks real and complex samples",
          "[baseband][dada]") {
    const std::array<std::string_view, 11> lines = {
        "NBIT 8", "NDIM 1", "NPOL 1", "NCHAN 2", "TSAMP 1.0", "FREQ 1400.0",
        "BW 2.0", "SOURCE J1234", "UTC_START 2020-01-01-00:00:00",
        "# a comment", "TELESCOPE GBT # ignored",};
    const std::vector<std::byte> payload{std::byte{10}, std::byte{200},};
    const TempFile file(make_dada(lines, payload));
    DadaReader reader(file.path());
    CHECK(reader.tsamp() == 1.0e-6);
    CHECK(reader.fch1() == 1399.5);
    CHECK(reader.foff() == 1.0);
    CHECK(reader.source_name() == "J1234");
    CHECK(reader.telescope() == "GBT");
    CHECK(reader.sample_type() == psrio::SampleType::kUInt8);
    auto samples = reader.read_samples<float>(1);
    REQUIRE(samples.size() == 2U);
    CHECK(samples[0] == 10.0F);
    CHECK(samples[1] == 200.0F);
    CHECK_THROWS_AS(reader.read_samples<std::complex<float>>(1),
                    psrio::ValidationError);

    const std::array<std::string_view, 5> complex_lines = {
        "NBIT 8", "NDIM 2", "NPOL 1", "NCHAN 1", "TSAMP 64.0",};
    const std::vector<std::byte> complex_payload{std::byte{1}, std::byte{2},};
    const TempFile complex_file(make_dada(complex_lines, complex_payload));
    auto complex_samples =
        DadaReader(complex_file.path()).read_samples<std::complex<std::int8_t>>(1);
    REQUIRE(complex_samples.size() == 1U);
    CHECK(complex_samples[0] == std::complex<std::int8_t>(1, 2));

    const std::vector<std::byte> wide{std::byte{0x02}, std::byte{0x01},
                                      std::byte{0x04}, std::byte{0x03},};
    const std::array<std::string_view, 5> wide_lines = {
        "NBIT 16", "NDIM 2", "NPOL 1", "NCHAN 1", "TSAMP 1.0",};
    const TempFile wide_file(make_dada(wide_lines, wide));
    auto wide_samples =
        DadaReader(wide_file.path()).read_samples<std::complex<std::int16_t>>(1);
    REQUIRE(wide_samples.size() == 1U);
    CHECK(wide_samples[0] == std::complex<std::int16_t>(0x0102, 0x0304));

    const std::vector<std::byte> two_bit{std::byte{0xE4}};
    const std::array<std::string_view, 5> two_lines = {
        "NBIT 2", "NDIM 1", "NPOL 1", "NCHAN 4", "TSAMP 1.0",};
    const TempFile two_file(make_dada(two_lines, two_bit));
    auto two_samples = DadaReader(two_file.path()).read_samples<std::int8_t>(1);
    REQUIRE(two_samples.size() == 4U);
    CHECK(two_samples[0] == 0);
    CHECK(two_samples[1] == 1);
    CHECK(two_samples[2] == 2);
    CHECK(two_samples[3] == 3);

    const std::vector<std::byte> four_bit{std::byte{0x21},};
    const std::array<std::string_view, 5> four_lines = {
        "NBIT 4", "NDIM 1", "NPOL 1", "NCHAN 2", "TSAMP 1.0",};
    const TempFile four_file(make_dada(four_lines, four_bit));
    auto four_samples = DadaReader(four_file.path()).read_samples<std::int8_t>(1);
    REQUIRE(four_samples.size() == 2U);
    CHECK(four_samples[0] == 1);
    CHECK(four_samples[1] == 2);

    const std::vector<std::byte> one_bit{std::byte{0xB2},};
    const std::array<std::string_view, 5> one_lines = {
        "NBIT 1", "NDIM 1", "NPOL 1", "NCHAN 8", "TSAMP 1.0",};
    const TempFile one_file(make_dada(one_lines, one_bit));
    auto one_samples = DadaReader(one_file.path()).read_samples<std::int8_t>(1);
    REQUIRE(one_samples.size() == 8U);
    CHECK(one_samples[0] == 0);
    CHECK(one_samples[1] == 1);
    CHECK(one_samples[4] == 1);
    CHECK(one_samples[7] == 1);
}

TEST_CASE("DadaReader transposes polarization and channel", "[baseband][dada]") {
    const std::vector<std::byte> payload{std::byte{1}, std::byte{2}, std::byte{3},
                                         std::byte{4},};
    const std::array<std::string_view, 6> lines = {
        "NBIT 8", "NDIM 1", "NPOL 2", "NCHAN 2", "TSAMP 1.0", "ORDER TF",};
    const TempFile file(make_dada(lines, payload));
    std::vector<std::byte> block(4);
    DadaReader reader(file.path());
    REQUIRE(reader.read_block(1, block) == 1U);
    const std::vector<std::byte> expected{std::byte{1}, std::byte{3}, std::byte{2},
                                          std::byte{4},};
    CHECK(block == expected);
}

TEST_CASE("DadaReader joins files on OBS_OFFSET and skips OBS_OVERLAP", "[baseband][dada]") {
    const std::array<std::string_view, 7> first_lines = {
        "NBIT 8", "NDIM 1", "NPOL 1", "NCHAN 1", "TSAMP 1.0", "OBS_OFFSET 0",
        "FILE_SIZE 4",};
    const std::array<std::string_view, 8> second_lines = {
        "NBIT 8", "NDIM 1", "NPOL 1", "NCHAN 1", "TSAMP 1.0", "OBS_OFFSET 2",
        "OBS_OVERLAP 2", "FILE_SIZE 4",};
    const std::vector<std::byte> first{std::byte{1}, std::byte{2}, std::byte{3},
                                       std::byte{4},};
    const std::vector<std::byte> second{std::byte{3}, std::byte{4}, std::byte{5},
                                        std::byte{6},};
    const TempFile file0(make_dada(first_lines, first));
    const TempFile file1(make_dada(second_lines, second));
    const std::array<std::filesystem::path, 2> paths = {file1.path(), file0.path()};
    DadaReader reader(paths);
    auto samples = reader.read_samples<std::int8_t>(6);
    REQUIRE(samples.size() == 6U);
    CHECK(samples[0] == 1);
    CHECK(samples[3] == 4);
    CHECK(samples[4] == 5);
    CHECK(samples[5] == 6);
}

TEST_CASE("DadaReader defaults optional NDIM and NCHAN and supports total FILE_SIZE",
          "[baseband][dada]") {
    // Missing NDIM and NCHAN should default to 1.
    // In make_dada, header_bytes = 256, payload = 4 bytes -> total file size = 260.
    const std::array<std::string_view, 4> lines = {
        "NBIT 8", "NPOL 1", "TSAMP 1.0", "FILE_SIZE 260",};
    const std::vector<std::byte> payload{std::byte{10}, std::byte{20}, std::byte{30},
                                         std::byte{40},};
    const TempFile file(make_dada(lines, payload));
    DadaReader reader(file.path());
    CHECK(reader.nchan() == 1U);
    CHECK(reader.ndim() == 1);
    CHECK(reader.nsamples() == 4U);
    auto samples = reader.read_samples<std::int8_t>(4);
    REQUIRE(samples.size() == 4U);
    CHECK(samples[0] == 10);
    CHECK(samples[3] == 40);
}

TEST_CASE("DadaReader rejects a bad header", "[baseband][dada]") {
    const std::vector<std::byte> tiny(16);
    const TempFile truncated(tiny);
    CHECK_THROWS_AS(DadaReader(truncated.path()), psrio::FormatError);

    std::string text = "HEADER FOO\nHDR_SIZE 64\n# end of header\n";
    std::vector<std::byte> bytes(64);
    std::memcpy(bytes.data(), text.data(), text.size());
    const TempFile wrong(bytes);
    CHECK_THROWS_AS(DadaReader(wrong.path()), psrio::FormatError);

    const std::array<std::string_view, 6> lines = {
        "NBIT 8", "NDIM 1", "NPOL 1", "NCHAN 1", "TSAMP 1.0", "ORDER FT",};
    const std::vector<std::byte> payload{std::byte{1}};
    const TempFile ft(make_dada(lines, payload));
    CHECK_THROWS_AS(DadaReader(ft.path()), psrio::FormatError);
}

TEST_CASE("BasebandSource erases GUPPI and DADA readers", "[baseband]") {
    const std::vector<std::byte> guppi_payload{std::byte{9}, std::byte{8}};
    const std::array<std::string_view, 7> cards = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 2",
        "OBSFREQ = 100.0", "CHAN_BW = 1.0", "PKTIDX  = 0",};
    const TempFile guppi_file(guppi_block(cards, guppi_payload));
    BasebandSource guppi(GuppiReader(guppi_file.path()));
    auto guppi_samples = guppi.read_samples<std::complex<std::int8_t>>(1);
    REQUIRE(guppi_samples.size() == 1U);
    CHECK(guppi_samples[0] == std::complex<std::int8_t>(9, 8));

    const std::array<std::string_view, 5> lines = {
        "NBIT 8", "NDIM 2", "NPOL 1", "NCHAN 1", "TSAMP 1.0",};
    const std::vector<std::byte> dada_payload{std::byte{4}, std::byte{5}};
    const TempFile dada_path(make_dada(lines, dada_payload));
    BasebandSource dada(DadaReader(dada_path.path()));
    auto dada_samples = dada.read_samples<std::complex<float>>(1);
    REQUIRE(dada_samples.size() == 1U);
    CHECK(dada_samples[0] == std::complex<float>(4.0F, 5.0F));
}

TEST_CASE("read_voltages matches FTPRI strides and ascending frequency",
          "[baseband][guppi]") {
    // Two channels, two times, two polarisations. On disk is FTPRI.
    // chan0 times: 10, 11. chan1 times: 20, 21. Each sample is 4 bytes.
    const std::vector<std::byte> payload{
        std::byte{10}, std::byte{1}, std::byte{2}, std::byte{3},
        std::byte{11}, std::byte{1}, std::byte{2}, std::byte{3},
        std::byte{20}, std::byte{1}, std::byte{2}, std::byte{3},
        std::byte{21}, std::byte{1}, std::byte{2}, std::byte{3},};
    const std::array<std::string_view, 7> rising = {
        "OBSNCHAN= 2", "NPOL    = 2", "NBITS   = 8", "BLOCSIZE= 16",
        "OBSFREQ = 1000.0", "CHAN_BW = 1.0", "PKTIDX  = 0",};
    const TempFile rising_file(guppi_block(rising, payload));
    GuppiReader rising_reader(rising_file.path());
    std::vector<std::byte> time_major(16);
    std::vector<std::byte> freq_major(16);
    REQUIRE(rising_reader.read_block(2, time_major) == 2U);
    rising_reader.seek(0);
    REQUIRE(rising_reader.read_voltages(
                2, {.order = psrio::VoltageOrder::kTimeMajor}, time_major) == 2U);
    const std::vector<std::byte> time_expected{
        std::byte{10}, std::byte{1}, std::byte{2}, std::byte{3},
        std::byte{20}, std::byte{1}, std::byte{2}, std::byte{3},
        std::byte{11}, std::byte{1}, std::byte{2}, std::byte{3},
        std::byte{21}, std::byte{1}, std::byte{2}, std::byte{3},};
    CHECK(time_major == time_expected);
    rising_reader.seek(0);
    REQUIRE(rising_reader.read_voltages(2, {}, freq_major) == 2U);
    CHECK(freq_major == payload);
    // FTPRI: channel stride is 4 elements * nsamps, time stride is 4.
    CHECK(freq_major[0] == std::byte{10});
    CHECK(freq_major[4] == std::byte{11});
    CHECK(freq_major[8] == std::byte{20});

    const std::array<std::string_view, 7> falling = {
        "OBSNCHAN= 2", "NPOL    = 2", "NBITS   = 8", "BLOCSIZE= 16",
        "OBSFREQ = 1000.0", "CHAN_BW = -1.0", "PKTIDX  = 0",};
    const TempFile falling_file(guppi_block(falling, payload));
    GuppiReader falling_reader(falling_file.path());
    std::vector<std::byte> reversed(16);
    REQUIRE(falling_reader.read_voltages(2, {}, reversed) == 2U);
    const std::vector<std::byte> reversed_expected{
        std::byte{20}, std::byte{1}, std::byte{2}, std::byte{3},
        std::byte{21}, std::byte{1}, std::byte{2}, std::byte{3},
        std::byte{10}, std::byte{1}, std::byte{2}, std::byte{3},
        std::byte{11}, std::byte{1}, std::byte{2}, std::byte{3},};
    CHECK(reversed == reversed_expected);
    CHECK(reversed[0] == std::byte{20});
    CHECK(reversed[8] == std::byte{10});

    // Two subbands, one time, negative spacing. Ascending sky frequency is
    // 999.5, 1000.5, 1003.5, 1004.5.
    const std::vector<std::byte> low_payload{
        std::byte{11}, std::byte{0}, std::byte{0}, std::byte{0},
        std::byte{22}, std::byte{0}, std::byte{0}, std::byte{0},};
    const std::vector<std::byte> high_payload{
        std::byte{33}, std::byte{0}, std::byte{0}, std::byte{0},
        std::byte{44}, std::byte{0}, std::byte{0}, std::byte{0},};
    const std::array<std::string_view, 8> low_cards = {
        "OBSNCHAN= 2", "NPOL    = 2", "NBITS   = 8", "BLOCSIZE= 8",
        "OBSFREQ = 1000.0", "CHAN_BW = -1.0", "TBIN    = 1.0e-6", "PKTIDX  = 0",};
    const std::array<std::string_view, 8> high_cards = {
        "OBSNCHAN= 2", "NPOL    = 2", "NBITS   = 8", "BLOCSIZE= 8",
        "OBSFREQ = 1004.0", "CHAN_BW = -1.0", "TBIN    = 1.0e-6", "PKTIDX  = 0",};
    const TempFile low_file(guppi_block(low_cards, low_payload));
    const TempFile high_file(guppi_block(high_cards, high_payload));
    const std::array<std::filesystem::path, 2> paths = {high_file.path(), low_file.path()};
    auto stitched = GuppiReader::open(paths);
    std::vector<std::byte> joined(16);
    REQUIRE(stitched.read_voltages(1, {}, joined) == 1U);
    const std::vector<std::byte> joined_expected{
        std::byte{22}, std::byte{0}, std::byte{0}, std::byte{0},
        std::byte{11}, std::byte{0}, std::byte{0}, std::byte{0},
        std::byte{44}, std::byte{0}, std::byte{0}, std::byte{0},
        std::byte{33}, std::byte{0}, std::byte{0}, std::byte{0},};
    CHECK(joined == joined_expected);

    const std::array<std::string_view, 7> single = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 2",
        "OBSFREQ = 100.0", "CHAN_BW = 1.0", "PKTIDX  = 0",};
    const TempFile one_pol(guppi_block(single, std::vector<std::byte>{std::byte{1}, std::byte{2}}));
    std::vector<std::byte> one_dest(2);
    CHECK_THROWS_AS(GuppiReader(one_pol.path()).read_voltages(1, {}, one_dest),
                    psrio::ValidationError);

    const std::array<std::string_view, 7> wide = {
        "OBSNCHAN= 1", "NPOL    = 2", "NBITS   = 16", "BLOCSIZE= 8",
        "OBSFREQ = 100.0", "CHAN_BW = 1.0", "PKTIDX  = 0",};
    const std::vector<std::byte> wide_payload(8);
    const TempFile wide_file(guppi_block(wide, wide_payload));
    std::vector<std::byte> wide_dest(8);
    CHECK_THROWS_AS(GuppiReader(wide_file.path()).read_voltages(1, {}, wide_dest),
                    psrio::ValidationError);

    const std::array<std::string_view, 8> antennas = {
        "OBSNCHAN= 2", "NANTS   = 2", "NPOL    = 2", "NBITS   = 8",
        "BLOCSIZE= 8", "OBSFREQ = 100.0", "CHAN_BW = 1.0", "PKTIDX  = 0",};
    const TempFile antenna_file(guppi_block(antennas, wide_payload));
    CHECK_THROWS_AS(GuppiReader(antenna_file.path()).read_voltages(1, {}, wide_dest),
                    psrio::ValidationError);
}

TEST_CASE("GuppiReader and FrequencyStitch stream native FTPRI layouts and subband groups",
          "[baseband][guppi]") {
    // 2 channels, 1 pol (complex -> 2 components per time), 8-bit, 2 time samples
    // Total blocksize = 2 chans * 2 times * 2 components * 1 byte = 8 bytes.
    // In GUPPI on-disk FTPRI:
    // chan 0: t0(r,i), t1(r,i) -> 4 bytes
    // chan 1: t0(r,i), t1(r,i) -> 4 bytes
    const std::vector<std::byte> payload_low{
        std::byte{10}, std::byte{11}, std::byte{12}, std::byte{13},
        std::byte{20}, std::byte{21}, std::byte{22}, std::byte{23},
    };
    const std::vector<std::byte> payload_high{
        std::byte{30}, std::byte{31}, std::byte{32}, std::byte{33},
        std::byte{40}, std::byte{41}, std::byte{42}, std::byte{43},
    };

    const std::array<std::string_view, 8> low_cards = {
        "OBSNCHAN= 2", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 8",
        "OBSFREQ = 1000.0", "CHAN_BW = 1.0", "TBIN    = 1.0e-6", "PKTIDX  = 0",};
    const std::array<std::string_view, 8> high_cards = {
        "OBSNCHAN= 2", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 8",
        "OBSFREQ = 1002.0", "CHAN_BW = 1.0", "TBIN    = 1.0e-6", "PKTIDX  = 0",};

    const TempFile low_file(guppi_block(low_cards, payload_low));
    const TempFile high_file(guppi_block(high_cards, payload_high));

    // 1. Single reader read_native
    GuppiReader single_reader(low_file.path());
    CHECK(single_reader.num_groups() == 1U);
    std::vector<std::byte> native_buf(8);
    REQUIRE(single_reader.read_native(2, native_buf) == 2U);
    CHECK(native_buf == payload_low);

    // 2. Multi-band FrequencyStitch read_groups
    const std::array<std::filesystem::path, 2> paths = {high_file.path(), low_file.path(),};
    auto stitched = GuppiReader::open(paths);
    CHECK(stitched.num_groups() == 2U);
    CHECK(stitched.nchan() == 4U);

    std::vector<std::byte> group0(8);
    std::vector<std::byte> group1(8);
    std::array<std::span<std::byte>, 2> group_spans = {group0, group1};
    REQUIRE(stitched.read_groups(2, group_spans) == 2U);
    CHECK(group0 == payload_low);
    CHECK(group1 == payload_high);

    // Test convenience allocating overload
    stitched.seek(0);
    auto alloc_groups = stitched.read_groups(2);
    REQUIRE(alloc_groups.size() == 2U);
    CHECK(alloc_groups[0] == payload_low);
    CHECK(alloc_groups[1] == payload_high);
}

TEST_CASE("FrequencyStitch orders subband groups in ascending frequency for negative foff",
          "[baseband][guppi]") {
    // low freq subband: OBSFREQ = 1000.0, CHAN_BW = -1.0, NCHAN = 2
    // channels: 1000.5 MHz (chan 0), 999.5 MHz (chan 1) -> sky freq 999.0 - 1001.0
    // high freq subband: OBSFREQ = 1004.0, CHAN_BW = -1.0, NCHAN = 2
    // channels: 1004.5 MHz (chan 0), 1003.5 MHz (chan 1) -> sky freq 1003.0 - 1005.0
    const std::vector<std::byte> payload_low{
        std::byte{1}, std::byte{2},
        std::byte{3}, std::byte{4},
    };
    const std::vector<std::byte> payload_high{
        std::byte{5}, std::byte{6},
        std::byte{7}, std::byte{8},
    };
    const std::array<std::string_view, 8> low_cards = {
        "OBSNCHAN= 2", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 4",
        "OBSFREQ = 1000.0", "CHAN_BW = -1.0", "TBIN    = 1.0e-6", "PKTIDX  = 0",};
    const std::array<std::string_view, 8> high_cards = {
        "OBSNCHAN= 2", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 4",
        "OBSFREQ = 1004.0", "CHAN_BW = -1.0", "TBIN    = 1.0e-6", "PKTIDX  = 0",};

    const TempFile low_file(guppi_block(low_cards, payload_low));
    const TempFile high_file(guppi_block(high_cards, payload_high));

    const std::array<std::filesystem::path, 2> paths = {high_file.path(), low_file.path(),};
    auto stitched = GuppiReader::open(paths);
    CHECK(stitched.num_groups() == 2U);

    // read_groups should order groups in ascending sky frequency order:
    // group 0 is the 1000 MHz subband, group 1 is the 1004 MHz subband.
    // Within each group, channels must also be in ascending frequency order:
    // payload_low has chan 0 (1000.5 MHz) as {1, 2} and chan 1 (999.5 MHz) as {3, 4}.
    // Ascending order places 999.5 MHz first ({3, 4}) then 1000.5 MHz ({1, 2}).
    const std::vector<std::byte> expected_group0{
        std::byte{3}, std::byte{4},
        std::byte{1}, std::byte{2},
    };
    const std::vector<std::byte> expected_group1{
        std::byte{7}, std::byte{8},
        std::byte{5}, std::byte{6},
    };
    auto groups = stitched.read_groups(1);
    REQUIRE(groups.size() == 2U);
    CHECK(groups[0] == expected_group0);
    CHECK(groups[1] == expected_group1);
}

TEST_CASE("Baseband readers unpack full range of signed 4-bit and 2-bit voltages",
          "[baseband][signed]") {
    // 4-bit two's complement range is -8 (0x8) to +7 (0x7).
    // GUPPI 4-bit is MSB-first: high nibble is first component, low nibble is second component.
    const std::vector<std::byte> payload4{
        std::byte{0x89}, std::byte{0xFE}, std::byte{0x07}, std::byte{0x16},};
    const std::array<std::string_view, 7> cards4 = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 4", "BLOCSIZE= 4",
        "OBSFREQ = 100.0", "CHAN_BW = 1.0", "PKTIDX  = 0",};
    const TempFile file4(guppi_block(cards4, payload4));
    GuppiReader reader4(file4.path());
    CHECK(reader4.msb_first());
    auto samps4 = reader4.read_samples<std::complex<std::int8_t>>(4);
    REQUIRE(samps4.size() == 4U);
    CHECK(samps4[0] == std::complex<std::int8_t>(-8, -7));
    CHECK(samps4[1] == std::complex<std::int8_t>(-1, -2));
    CHECK(samps4[2] == std::complex<std::int8_t>(0, 7));
    CHECK(samps4[3] == std::complex<std::int8_t>(1, 6));

    // 2-bit two's complement:
    // 00 -> 0, 01 -> 1, 10 -> -2, 11 -> -1
    // With 2 channels and 2 time samples, blocsize is 2 bytes (1 byte per channel):
    // Channel 0 byte (0x4E = 0b01001110):
    //   t0: bits 1..0: 10 (-2), bits 3..2: 11 (-1) -> (-2, -1)
    //   t1: bits 5..4: 00 (0),  bits 7..6: 01 (+1) -> (0, 1)
    // Channel 1 byte (0x39 = 0b00111001):
    //   t0: bits 1..0: 01 (+1), bits 3..2: 10 (-2) -> (1, -2)
    //   t1: bits 5..4: 11 (-1), bits 7..6: 00 (0)  -> (-1, 0)
    const std::vector<std::byte> payload2{std::byte{0x4E}, std::byte{0x39},};
    const std::array<std::string_view, 7> cards2 = {
        "OBSNCHAN= 2", "NPOL    = 1", "NBITS   = 2", "BLOCSIZE= 2",
        "OBSFREQ = 100.0", "CHAN_BW = 1.0", "PKTIDX  = 0",};
    const TempFile file2(guppi_block(cards2, payload2));
    GuppiReader reader2(file2.path());
    CHECK_FALSE(reader2.msb_first());
    auto samps2 = reader2.read_samples<std::complex<std::int8_t>>(2);
    REQUIRE(samps2.size() == 4U);
    CHECK(samps2[0] == std::complex<std::int8_t>(-2, -1));
    CHECK(samps2[1] == std::complex<std::int8_t>(1, -2));
    CHECK(samps2[2] == std::complex<std::int8_t>(0, 1));
    CHECK(samps2[3] == std::complex<std::int8_t>(-1, 0));
}

TEST_CASE("GuppiReader read_native zero-fills dropped packets", "[baseband][guppi]") {
    const std::vector<std::byte> first{std::byte{0xAA}, std::byte{0xBB},};
    const std::vector<std::byte> second{std::byte{0xCC}, std::byte{0xDD},};
    const std::array<std::string_view, 7> head = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 2",
        "OBSFREQ = 1400.0", "CHAN_BW = 1.0", "PKTIDX  = 0",};
    const std::array<std::string_view, 7> jumped = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 2",
        "OBSFREQ = 1400.0", "CHAN_BW = 1.0", "PKTIDX  = 3",};
    const TempFile file0(guppi_block(head, first));
    const TempFile gap(guppi_block(jumped, second));
    const std::array<std::filesystem::path, 2> paths = {file0.path(), gap.path(),};
    GuppiReader reader(paths);
    CHECK(reader.dropped_packets() == 2U);
    CHECK(reader.dropped_samples() == 2U);
    CHECK(reader.nsamples() == 4U);

    std::vector<std::byte> dest(8);
    REQUIRE(reader.read_native(4, dest) == 4U);
    CHECK(dest[0] == std::byte{0xAA});
    CHECK(dest[1] == std::byte{0xBB});
    CHECK(dest[2] == std::byte{0x00});
    CHECK(dest[3] == std::byte{0x00});
    CHECK(dest[4] == std::byte{0x00});
    CHECK(dest[5] == std::byte{0x00});
    CHECK(dest[6] == std::byte{0xCC});
    CHECK(dest[7] == std::byte{0xDD});
}

TEST_CASE("GuppiReader fills VDIF packet gaps with zero voltage", "[baseband][guppi]") {
    // In VDIF offset-binary, on-disk 0x80 is 0 voltage, 0x00 is -128.
    // Gaps must be filled with signed 0 (0x00), NOT 0x80 (-128).
    const std::vector<std::byte> first{std::byte{0x80}, std::byte{0x81},};
    const std::vector<std::byte> second{std::byte{0x82}, std::byte{0x83},};
    const std::array<std::string_view, 8> head = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 2",
        "OBSFREQ = 1400.0", "CHAN_BW = 1.0", "PKTFMT  = 'VDIF'", "PKTIDX  = 0",};
    const std::array<std::string_view, 8> jumped = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 2",
        "OBSFREQ = 1400.0", "CHAN_BW = 1.0", "PKTFMT  = 'VDIF'", "PKTIDX  = 2",};
    const TempFile file0(guppi_block(head, first));
    const TempFile gap(guppi_block(jumped, second));
    const std::array<std::filesystem::path, 2> paths = {file0.path(), gap.path(),};
    GuppiReader reader(paths);
    CHECK(reader.dropped_packets() == 1U);
    CHECK(reader.nsamples() == 3U);

    auto samps = reader.read_samples<std::complex<std::int8_t>>(3);
    REQUIRE(samps.size() == 3U);
    // 0x80 ^ 0x80 = 0, 0x81 ^ 0x80 = 1
    CHECK(samps[0] == std::complex<std::int8_t>(0, 1));
    // Gap sample must be signed 0, NOT -128!
    CHECK(samps[1] == std::complex<std::int8_t>(0, 0));
    // 0x82 ^ 0x80 = 2, 0x83 ^ 0x80 = 3
    CHECK(samps[2] == std::complex<std::int8_t>(2, 3));
}

TEST_CASE("GuppiReader rejects absurdly large PKTIDX jumps", "[baseband][guppi]") {
    const std::vector<std::byte> payload{std::byte{1}, std::byte{2},};
    const std::array<std::string_view, 7> head = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 2",
        "OBSFREQ = 1400.0", "CHAN_BW = 1.0", "PKTIDX  = 0",};
    const std::array<std::string_view, 7> absurd = {
        "OBSNCHAN= 1", "NPOL    = 1", "NBITS   = 8", "BLOCSIZE= 2",
        "OBSFREQ = 1400.0", "CHAN_BW = 1.0", "PKTIDX  = 50000000",};
    const TempFile file0(guppi_block(head, payload));
    const TempFile gap(guppi_block(absurd, payload));
    const std::array<std::filesystem::path, 2> paths = {file0.path(), gap.path(),};
    CHECK_THROWS_AS(GuppiReader(paths), psrio::FormatError);
}
