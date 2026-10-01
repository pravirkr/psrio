#include "psrio/detail/exceptions.hpp"
#include "psrio/formats/sigproc/header.hpp"
#include "psrio/formats/sigproc/reader.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <span>
#include <string>
#include <vector>

// Files in tests/data are the sigpyproc3 copies. Header checks follow
// sigpyproc3 parse_header. Packed samples are least-significant field
// first, the DSPSR order. sigpyproc3 unpacks 2-bit and 4-bit the other
// way inside each byte.
#ifndef PSRIO_TEST_DATA_DIR
#error PSRIO_TEST_DATA_DIR must be the tests/data directory
#endif

using psrio::formats::sigproc::FilterbankHeader;
using psrio::formats::sigproc::FilterbankReader;

namespace {

enum class Family : std::uint8_t { kParkes, kTutorial };

struct PayloadExpect {
    const char* name;
    std::uint64_t data_bytes;
    std::uint64_t nsamples;
    std::uint64_t sample0_sum;
    std::uint64_t unpacked_sum;
    Family family;
    std::int32_t nbits;
    std::uint32_t payload_crc;
    std::uint8_t sample0_last;
    std::array<std::uint8_t, 16> first;
};

constexpr auto
    kFiles =
        std::to_array<PayloadExpect>(
            {
                {
                    .name         = "parkes_1bit.fil",
                    .data_bytes   = 425984,
                    .nsamples     = 4096,
                    .sample0_sum  = 406,
                    .unpacked_sum = 1707422,
                    .family       = Family::kParkes,
                    .nbits        = 1,
                    .payload_crc  = 330239577U,
                    .sample0_last = 1,
                    .first = {0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1},
                },
                {
                    .name         = "parkes_2bit.fil",
                    .data_bytes   = 851968,
                    .nsamples     = 4096,
                    .sample0_sum  = 1250,
                    .unpacked_sum = 5150964,
                    .family       = Family::kParkes,
                    .nbits        = 2,
                    .payload_crc  = 3580086145U,
                    .sample0_last = 2,
                    .first = {1, 2, 2, 0, 2, 2, 2, 3, 3, 2, 2, 2, 1, 3, 3, 3},
                },
                {
                    .name         = "parkes_4bit.fil",
                    .data_bytes   = 1703936,
                    .nsamples     = 4096,
                    .sample0_sum  = 6239,
                    .unpacked_sum = 25658391,
                    .family       = Family::kParkes,
                    .nbits        = 4,
                    .payload_crc  = 1750014171U,
                    .sample0_last = 8,
                    .first = {7, 8, 8, 5, 9, 9, 8, 9, 9, 9, 8, 8, 7, 9, 10, 9},
                },
                {
                    .name         = "parkes_8bit_1.fil",
                    .data_bytes   = 3407872,
                    .nsamples     = 4096,
                    .sample0_sum  = 106336,
                    .unpacked_sum = 436194585,
                    .family       = Family::kParkes,
                    .nbits        = 8,
                    .payload_crc  = 195344466U,
                    .sample0_last = 137,
                    .first =
                        {
                            113,
                            136,
                            142,
                            92,
                            146,
                            148,
                            129,
                            151,
                            155,
                            147,
                            140,
                            130,
                            120,
                            158,
                            173,
                            157,
                        },
                },
                {
                    .name         = "tutorial.fil",
                    .data_bytes   = 3000320,
                    .nsamples     = 187520,
                    .sample0_sum  = 107,
                    .unpacked_sum = 19638770,
                    .family       = Family::kTutorial,
                    .nbits        = 2,
                    .payload_crc  = 2371190182U,
                    .sample0_last = 2,
                    .first = {2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 2, 1, 1, 2, 1, 2},
                },
                {
                    .name         = "tutorial_2bit.fil",
                    .data_bytes   = 3000320,
                    .nsamples     = 187520,
                    .sample0_sum  = 107,
                    .unpacked_sum = 19638770,
                    .family       = Family::kTutorial,
                    .nbits        = 2,
                    .payload_crc  = 2371190182U,
                    .sample0_last = 2,
                    .first = {2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 2, 1, 1, 2, 1, 2},
                },
            });

[[nodiscard]] std::filesystem::path data_path(const char* name) {
    return std::filesystem::path{PSRIO_TEST_DATA_DIR} / name;
}

[[nodiscard]] std::uint32_t crc32(std::span<const std::byte> bytes) {
    std::uint32_t crc = 0xFFFFFFFFU;
    for (const std::byte value : bytes) {
        crc ^= std::to_integer<unsigned char>(value);
        for (int bit = 0; bit < 8; ++bit) {
            const auto mask = static_cast<std::uint32_t>(-(crc & 1U));
            crc             = (crc >> 1U) ^ (0xEDB88320U & mask);
        }
    }
    return ~crc;
}

[[nodiscard]] std::uint64_t sample_sum(std::span<const std::uint8_t> values) {
    std::uint64_t sum = 0;
    for (const std::uint8_t value : values) {
        sum += value;
    }
    return sum;
}

void check_common(const FilterbankHeader& header) {
    REQUIRE(header.nifs == 1);
    REQUIRE(header.data_type == 1);
    REQUIRE(header.data_type_name() == "filterbank");
    REQUIRE(header.machine_id == 0);
    REQUIRE(header.machine_name() == "FAKE");
    REQUIRE_FALSE(header.samples_are_signed());
    REQUIRE_FALSE(header.signed_data.has_value());
    REQUIRE_FALSE(header.declared_nsamples.has_value());
    REQUIRE_FALSE(header.barycentric.has_value());
    REQUIRE_FALSE(header.pulsarcentric.has_value());
    REQUIRE_FALSE(header.refdm.has_value());
    REQUIRE_FALSE(header.period.has_value());
    REQUIRE(header.frequency_table.empty());
    REQUIRE(header.extra.empty());
    REQUIRE(header.warnings().empty());
    REQUIRE(header.trailing_bytes() == 0U);
}

void check_parkes(const FilterbankHeader& header) {
    check_common(header);
    REQUIRE(header.header_bytes == 351U);
    REQUIRE(header.nchans == 832);
    REQUIRE(header.telescope_id == 4);
    REQUIRE(header.telescope_name() == "Parkes");
    REQUIRE(header.source_name == "J0534+2200");
    REQUIRE(header.rawdatafile == "unknown");
    REQUIRE(header.src_raj == 0x1.a16fccccccccdp+15);
    REQUIRE(header.src_dej == 220052.0);
    REQUIRE(header.az_start == 0.0);
    REQUIRE(header.za_start == 0.0);
    REQUIRE(header.ibeam == 0);
    REQUIRE(header.nbeams == 0);
    REQUIRE(header.fch1 == 4030.0);
    REQUIRE(header.foff == -4.0);
    REQUIRE(header.tsamp == 0x1.0c6f7a0b5ed8dp-11);
    REQUIRE(header.tstart == 0x1.c95ea92884393p+15);
    REQUIRE(header.bandwidth() == 3328.0);
    REQUIRE(header.ftop() == 4032.0);
    REQUIRE(header.fbottom() == 704.0);
    REQUIRE(header.fcenter() == 2368.0);

    std::vector<double> channels(static_cast<std::size_t>(header.nchans));
    header.channel_frequencies(channels);
    REQUIRE(channels.front() == 4030.0);
    REQUIRE(channels.at(1) == 4026.0);
    REQUIRE(channels.back() == 706.0);
}

void check_tutorial(const FilterbankHeader& header) {
    check_common(header);
    REQUIRE(header.header_bytes == 244U);
    REQUIRE(header.nchans == 64);
    REQUIRE(header.nbits == 2);
    REQUIRE(header.nsamples() == 187520U);
    REQUIRE(header.telescope_id == 0);
    REQUIRE(header.telescope_name() == "Fake");
    REQUIRE(header.source_name == "P: 250.000000000000 ms, DM: 30.000");
    REQUIRE_FALSE(header.rawdatafile.has_value());
    REQUIRE_FALSE(header.src_raj.has_value());
    REQUIRE_FALSE(header.src_dej.has_value());
    REQUIRE_FALSE(header.az_start.has_value());
    REQUIRE_FALSE(header.za_start.has_value());
    REQUIRE_FALSE(header.ibeam.has_value());
    REQUIRE_FALSE(header.nbeams.has_value());
    REQUIRE(header.fch1 == 1510.0);
    REQUIRE(header.foff == -0x1.170a3d70a3d71p+0);
    REQUIRE(header.tsamp == 0x1.4f8b588e368f1p-12);
    REQUIRE(header.tstart == 50000.0);
    REQUIRE(header.bandwidth() == 0x1.170a3d70a3d71p+6);
    REQUIRE(header.ftop() == 0x1.79a2e147ae148p+10);
    REQUIRE(header.fbottom() == 0x1.68323d70a3d71p+10);
    REQUIRE(header.fcenter() == 0x1.70ea8f5c28f5cp+10);

    std::vector<double> channels(static_cast<std::size_t>(header.nchans));
    header.channel_frequencies(channels);
    REQUIRE(channels.front() == 1510.0);
    REQUIRE(channels.at(1) == 0x1.793a3d70a3d71p+10);
    REQUIRE(channels.back() == 0x1.68551eb851eb8p+10);
}

void check_reads(FilterbankReader& reader, const PayloadExpect& expect) {
    const FilterbankHeader& header = reader.header();
    REQUIRE(header.nbits == expect.nbits);
    REQUIRE(header.data_bytes == expect.data_bytes);
    REQUIRE(header.nsamples() == expect.nsamples);
    REQUIRE(header.samples_in_file() == expect.nsamples);
    REQUIRE(header.file_bytes == header.header_bytes + header.data_bytes);
    REQUIRE(header.bytes_per_sample() == expect.data_bytes / expect.nsamples);

    const auto nchans = static_cast<std::size_t>(header.nchans);
    std::vector<std::uint8_t> sample(nchans);
    REQUIRE(reader.read(1, std::span<std::uint8_t>{sample}) == 1U);
    REQUIRE(
        std::equal(expect.first.begin(), expect.first.end(), sample.begin()));
    REQUIRE(sample.back() == expect.sample0_last);
    REQUIRE(sample_sum(sample) == expect.sample0_sum);

    reader.rewind();
    std::vector<float> as_float(nchans);
    REQUIRE(reader.read(1, std::span<float>{as_float}) == 1U);
    REQUIRE(as_float.front() == static_cast<float>(expect.first.front()));
    REQUIRE(as_float.back() == static_cast<float>(expect.sample0_last));

    reader.rewind();
    constexpr std::uint64_t kGulp = 512;
    std::vector<std::uint8_t> block(static_cast<std::size_t>(kGulp) * nchans);
    std::uint64_t seen = 0;
    std::uint64_t sum  = 0;
    while (seen < expect.nsamples) {
        const auto count = reader.read(kGulp, std::span<std::uint8_t>{block});
        REQUIRE(count > 0U);
        REQUIRE(count <= kGulp);
        if (seen + count < expect.nsamples) {
            REQUIRE(count == kGulp);
        }
        sum += sample_sum(std::span<const std::uint8_t>{
            block.data(), static_cast<std::size_t>(count) * nchans});
        seen += count;
    }
    REQUIRE(seen == expect.nsamples);
    REQUIRE(sum == expect.unpacked_sum);
    REQUIRE(reader.read(kGulp, std::span<std::uint8_t>{block}) == 0U);
    REQUIRE(reader.tell() == expect.nsamples);
    reader.seek(expect.nsamples);
    REQUIRE_THROWS_AS(reader.seek(expect.nsamples + 1U),
                      psrio::ValidationError);

    reader.rewind();
    const auto stride = header.bytes_per_sample();
    const auto viewed = reader.view(stride);
    REQUIRE(viewed.size() == stride);
    REQUIRE(reader.tell() == 1U);
    REQUIRE_THROWS_AS(reader.view(1), psrio::ValidationError);
    REQUIRE(reader.tell() == 1U);

    reader.rewind();
    std::vector<std::byte> payload(static_cast<std::size_t>(expect.data_bytes));
    REQUIRE(reader.read_bytes(expect.data_bytes, payload) == expect.data_bytes);
    REQUIRE(crc32(payload) == expect.payload_crc);
    REQUIRE(std::equal(viewed.begin(), viewed.end(), payload.begin()));
    REQUIRE(reader.tell() == expect.nsamples);
}

} // namespace

TEST_CASE("crc32 matches the ISO-HDLC check value", "[sigproc][data]") {
    const std::string digits = "123456789";
    std::vector<std::byte> bytes(digits.size());
    for (std::size_t index = 0; index < digits.size(); ++index) {
        bytes.at(index) = static_cast<std::byte>(digits.at(index));
    }
    REQUIRE(crc32(bytes) == 0xCBF43926U);
}

TEST_CASE("sigpyproc3 filterbank files round-trip through the reader",
          "[sigproc][data]") {
    for (const PayloadExpect& expect : kFiles) {
        SECTION(expect.name) {
            const auto path = data_path(expect.name);
            REQUIRE(std::filesystem::is_regular_file(path));
            FilterbankReader reader(path);
            if (expect.family == Family::kParkes) {
                check_parkes(reader.header());
            } else {
                check_tutorial(reader.header());
            }
            check_reads(reader, expect);
            REQUIRE(std::filesystem::file_size(path) ==
                    reader.header().file_bytes);
        }
    }
}
