#include "psrfits_bytes.hpp"
#include "psrio/astro.hpp"
#include "psrio/block_source.hpp"
#include "psrio/common/types.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/digifil.hpp"
#include "psrio/formats/psrfits.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <span>
#include <string>
#include <string_view>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using psrio::formats::psrfits::SearchFile;
using psrio::test::append_be;
using psrio::test::FitsFileBytes;

namespace {

class TempPath {
public:
    explicit TempPath(std::string_view name) {
        m_path = std::filesystem::temp_directory_path() / name;
        std::filesystem::remove(m_path);
    }
    ~TempPath() { std::filesystem::remove(m_path); }
    TempPath(const TempPath&)            = delete;
    TempPath& operator=(const TempPath&) = delete;
    TempPath(TempPath&&)                 = delete;
    TempPath& operator=(TempPath&&)      = delete;
    [[nodiscard]] const std::filesystem::path& path() const { return m_path; }

private:
    std::filesystem::path m_path;
};

void primary_cards(FitsFileBytes& file, std::string_view mode = "SEARCH") {
    file.begin_hdu();
    file.key("SIMPLE", "T");
    file.key("BITPIX", "8");
    file.key("NAXIS", "0");
    file.key("EXTEND", "T");
    file.key("FITSTYPE", "'PSRFITS'");
    file.key("OBS_MODE", std::string("'") + std::string(mode) + "'");
    file.key("TELESCOP", "'Parkes'");
    file.key("SRC_NAME", "'Test'");
    file.key("RA", "'01:02:03.5'");
    file.key("DEC", "'+04:05:06.7'");
    file.key("STT_IMJD", "50000");
    file.key("STT_SMJD", "3600");
    file.key("STT_OFFS", "0.5");
    file.key("CHAN_DM", "12.5");
    file.key("IBEAM", "2");
    file.end_header();
}

struct Geometry {
    int nchan{2};
    int npol{4};
    int nbits{8};
    int nsblk{2};
    int nsub{1};
    std::string pol{"AABBCRCI"};
    double tbin{1.0e-3};
    double zero{0.0};
    int signint{0};
    std::string tdim;
};

void subint_header(FitsFileBytes& file,
                   const Geometry& geo,
                   std::uint64_t row_bytes,
                   std::string_view extra_tdim = {}) {
    file.begin_hdu();
    file.key("XTENSION", "'BINTABLE'");
    file.key("BITPIX", "8");
    file.key("NAXIS", "2");
    file.key("NAXIS1", std::to_string(row_bytes));
    file.key("NAXIS2", std::to_string(geo.nsub));
    file.key("PCOUNT", "0");
    file.key("GCOUNT", "1");
    file.key("TFIELDS", "7");
    file.key("EXTNAME", "'SUBINT'");
    file.key("TTYPE1", "'TSUBINT'");
    file.key("TFORM1", "'1D'");
    file.key("TTYPE2", "'OFFS_SUB'");
    file.key("TFORM2", "'1D'");
    file.key("TTYPE3", "'DAT_FREQ'");
    file.key("TFORM3", std::to_string(geo.nchan) + "D");
    file.key("TTYPE4", "'DAT_WTS'");
    file.key("TFORM4", std::to_string(geo.nchan) + "E");
    file.key("TTYPE5", "'DAT_OFFS'");
    file.key("TFORM5", std::to_string(geo.npol * geo.nchan) + "E");
    file.key("TTYPE6", "'DAT_SCL'");
    file.key("TFORM6", std::to_string(geo.npol * geo.nchan) + "E");
    const auto elements =
        static_cast<std::uint64_t>(geo.nchan) * geo.npol * geo.nsblk;
    const auto data_count =
        geo.nbits == 16 ? elements
                        : (elements * static_cast<unsigned>(geo.nbits)) / 8U;
    const char* data_code = geo.nbits == 16 ? "I" : "B";
    file.key("TTYPE7", "'DATA'");
    file.key("TFORM7", std::to_string(data_count) + data_code);
    if (!extra_tdim.empty()) {
        file.key("TDIM7", extra_tdim);
    }
    file.key("NPOL", std::to_string(geo.npol));
    file.key("NCHAN", std::to_string(geo.nchan));
    file.key("NBITS", std::to_string(geo.nbits));
    file.key("NSBLK", std::to_string(geo.nsblk));
    file.key("TBIN", std::to_string(geo.tbin));
    file.key("POL_TYPE", std::string("'") + geo.pol + "'");
    file.key("ZERO_OFF", std::to_string(geo.zero));
    file.key("SIGNINT", std::to_string(geo.signint));
    file.end_header();
}

void push_row_prefix(std::vector<std::byte>& row,
                     const Geometry& geo,
                     std::span<const double> freqs,
                     std::span<const float> weights,
                     std::span<const float> offsets,
                     std::span<const float> scales,
                     double tsubint,
                     double offs_sub) {
    append_be(row, tsubint);
    append_be(row, offs_sub);
    for (const double freq : freqs) {
        append_be(row, freq);
    }
    for (const float weight : weights) {
        append_be(row, weight);
    }
    for (const float offset : offsets) {
        append_be(row, offset);
    }
    for (const float scale : scales) {
        append_be(row, scale);
    }
    REQUIRE(freqs.size() == static_cast<std::size_t>(geo.nchan));
    REQUIRE(weights.size() == static_cast<std::size_t>(geo.nchan));
    REQUIRE(offsets.size() == static_cast<std::size_t>(geo.npol) * geo.nchan);
    REQUIRE(scales.size() == static_cast<std::size_t>(geo.npol) * geo.nchan);
}

} // namespace

TEST_CASE("SearchFile streams calibrated Stokes I", "[psrfits]") {
    static_assert(psrio::concepts::BlockReader<SearchFile>);

    Geometry geo;
    geo.zero = 10.0;
    const std::vector<double> freqs{100.0, 100.5};
    const std::vector<float> weights{1.0F, 0.5F};
    const std::vector<float> offsets(8, 0.0F);
    const std::vector<float> scales(8, 1.0F);
    // Two samples, four products, two channels. Only the first two products
    // form Stokes I. Layout is time, pol, channel.
    const std::array<std::uint8_t, 16> row0{
        12, 14, 18, 20, 0, 0, 0, 0, // sample 0
        11, 11, 11, 11, 9, 9, 9, 9, // sample 1
    };
    const std::vector<float> row1_scales{2, 2, 2, 2, 0, 0, 0, 0};
    const std::array<std::uint8_t, 16> row1{
        15, 15, 15, 15, 7, 7, 7, 7, 15, 15, 15, 15, 7, 7, 7, 7,
    };

    std::vector<std::byte> table;
    push_row_prefix(table, geo, freqs, weights, offsets, scales, 0.002, 0.001);
    for (const auto byte : row0) {
        append_be(table, byte);
    }
    push_row_prefix(table, geo, freqs, weights, offsets, row1_scales, 0.002,
                    0.003);
    for (const auto byte : row1) {
        append_be(table, byte);
    }
    geo.nsub = 2;

    FitsFileBytes file;
    primary_cards(file);
    subint_header(file, geo, table.size() / 2U, "'(2,4,2)'");
    file.raw(table);
    file.pad_data();

    const TempPath path("psrio-psrfits-8.sf");
    file.write(path.path());
    SearchFile reader(path.path());

    CHECK(reader.nchans() == 2U);
    CHECK(reader.nifs() == 1U);
    CHECK(reader.nbits() == 32);
    CHECK(reader.input_nbits() == 8);
    CHECK(reader.sample_type() == psrio::SampleType::kFloat32);
    CHECK(reader.nsamples() == 4U);
    CHECK(reader.pol_type() == "AABBCRCI");
    CHECK(reader.source_name() == "Test");
    CHECK(reader.telescope() == "Parkes");
    CHECK(reader.beam() == 2);
    CHECK(reader.fch1() == Catch::Approx(100.0));
    CHECK(reader.foff() == Catch::Approx(0.5));
    CHECK(reader.tsamp() == Catch::Approx(1.0e-3));
    CHECK(reader.tstart() == Catch::Approx(50000.0 + (3600.5 / 86400.0)));
    CHECK(reader.raj() == Catch::Approx(psrio::astro::string_to_packed(
                              "01:02:03.5", "right ascension")));
    CHECK(reader.dej() == Catch::Approx(psrio::astro::string_to_packed(
                              "+04:05:06.7", "declination")));
    CHECK(reader.sigproc_header().refdm == Catch::Approx(12.5));
    CHECK(reader.sigproc_header().telescope_id == 4);
    CHECK(reader.bytes_per_sample() == 8U);

    const auto samples = reader.read_samples<float>(3);
    REQUIRE(samples.size() == 6U);
    // (12 - 10) + (18 - 10) = 10; channel 1 weight 0.5.
    CHECK(samples[0] == 10.0F);
    CHECK(samples[1] == 7.0F);
    CHECK(samples[2] == 2.0F);
    CHECK(samples[3] == 1.0F);
    // Row 1 scale is 2: (15 - 10) * 2, summed, then the weight.
    CHECK(samples[4] == 20.0F);
    CHECK(samples[5] == 10.0F);
    CHECK(reader.tell() == 3U);

    reader.seek(0);
    reader.set_fswap(true);
    const auto swapped = reader.read_samples<float>(1);
    CHECK(swapped[0] == 7.0F);
    CHECK(swapped[1] == 10.0F);
    reader.set_fswap(false);

    reader.seek(0);
    auto block = reader.read_block(1);
    REQUIRE(block.size() == 8U);
    float first = 0.0F;
    std::memcpy(&first, block.data(), sizeof(first));
    CHECK(first == 10.0F);

    reader.seek(4);
    CHECK(reader.read_samples<float>(2).empty());
    CHECK_THROWS_AS(reader.read_bytes(8, std::span<std::byte>{}),
                    psrio::ValidationError);
    reader.seek(3);
    std::vector<std::byte> too_far(16);
    CHECK_THROWS_AS(reader.read_bytes(16, too_far), psrio::ValidationError);
    CHECK_THROWS_AS(reader.read_samples(1, std::span<std::uint8_t>{}),
                    psrio::ValidationError);
}

TEST_CASE("SearchFile unpacks 2-bit samples with dspsr levels", "[psrfits]") {
    Geometry geo;
    geo.nchan = 4;
    geo.npol  = 2;
    geo.nbits = 2;
    geo.nsblk = 4;
    geo.pol   = "AABB";
    geo.zero  = 127.5;
    geo.tbin  = 1.0;
    const std::vector<double> freqs{400.0, 399.0, 398.0, 397.0};
    const std::vector<float> weights(4, 1.0F);
    const std::vector<float> offsets(8, 0.0F);
    const std::vector<float> scales(8, 1.0F);
    // One byte holds four 2-bit channels, MSB first: 3, 2, 1, 0.
    const std::uint8_t packed = 0xE4;
    std::vector<std::byte> table;
    push_row_prefix(table, geo, freqs, weights, offsets, scales, 4.0, 2.0);
    for (int sample = 0; sample < 4; ++sample) {
        append_be(table, packed);
    }
    for (int sample = 0; sample < 4; ++sample) {
        append_be(table, packed);
    }

    FitsFileBytes file;
    primary_cards(file);
    subint_header(file, geo, table.size(), "'(4,2,1)'");
    file.raw(table);
    file.pad_data();
    const TempPath path("psrio-psrfits-2bit.sf");
    file.write(path.path());

    SearchFile reader(path.path());
    CHECK(reader.input_nbits() == 2);
    CHECK(reader.foff() == Catch::Approx(-1.0));
    const auto samples = reader.read_samples<float>(4);
    REQUIRE(samples.size() == 16U);
    const std::array<float, 4> expected{4.0F, 1.0F, -1.0F, -4.0F};
    for (int sample = 0; sample < 4; ++sample) {
        for (int chan = 0; chan < 4; ++chan) {
            CHECK(samples[static_cast<std::size_t>((sample * 4) + chan)] ==
                  expected[static_cast<std::size_t>(chan)]);
        }
    }
}

TEST_CASE("SearchFile keeps an irregular frequency table", "[psrfits]") {
    Geometry geo;
    geo.nchan = 3;
    geo.npol  = 1;
    geo.pol   = "AA+BB";
    geo.nsblk = 1;
    const std::vector<double> freqs{100.0, 101.0, 103.0};
    const std::vector<float> weights{1.0F, 1.0F, 1.0F};
    const std::vector<float> offsets{0.0F, 0.0F, 0.0F};
    const std::vector<float> scales{1.0F, 1.0F, 1.0F};
    std::vector<std::byte> table;
    push_row_prefix(table, geo, freqs, weights, offsets, scales, 1.0e-3,
                    5.0e-4);
    append_be(table, std::uint8_t{5});
    append_be(table, std::uint8_t{6});
    append_be(table, std::uint8_t{7});

    FitsFileBytes file;
    primary_cards(file);
    subint_header(file, geo, table.size());
    file.raw(table);
    file.pad_data();
    const TempPath path("psrio-psrfits-irregular.sf");
    file.write(path.path());

    SearchFile reader(path.path());
    CHECK(reader.foff() == 0.0);
    CHECK(reader.fch1() == Catch::Approx(100.0));
    const auto& header = reader.sigproc_header();
    REQUIRE(header.frequency_table.size() == 3U);
    CHECK(header.frequency_table[2] == Catch::Approx(103.0));
    CHECK_FALSE(header.foff.has_value());
    const auto samples = reader.read_samples<float>(1);
    CHECK(samples == std::vector<float>{5.0F, 6.0F, 7.0F});
}

TEST_CASE("SearchFile honours SIGNINT and a transposed DATA axis",
          "[psrfits]") {
    Geometry geo;
    geo.nchan   = 3;
    geo.npol    = 4;
    geo.nsblk   = 5;
    geo.signint = 1;
    geo.tbin    = 1.0;
    geo.pol     = "AABBCRCI";
    const std::vector<double> freqs{10.0, 11.0, 12.0};
    const std::vector<float> weights(3, 1.0F);
    const std::vector<float> offsets(12, 0.0F);
    const std::vector<float> scales(12, 1.0F);
    std::vector<std::uint8_t> data(60, 0);
    // Fortran TDIM (5, 4, 3) is C order (chan, pol, time) with time fastest.
    const auto at = [](int chan, int pol, int sample) {
        return (((chan * 4) + pol) * 5) + sample;
    };
    data[static_cast<std::size_t>(at(0, 0, 0))] = 200;
    data[static_cast<std::size_t>(at(0, 1, 0))] = 10;
    data[static_cast<std::size_t>(at(2, 0, 1))] = 3;
    data[static_cast<std::size_t>(at(2, 1, 1))] = 4;

    std::vector<std::byte> table;
    push_row_prefix(table, geo, freqs, weights, offsets, scales, 5.0, 2.5);
    for (const auto byte : data) {
        append_be(table, byte);
    }
    FitsFileBytes file;
    primary_cards(file);
    subint_header(file, geo, table.size(), "'(5,4,3)'");
    file.raw(table);
    file.pad_data();
    const TempPath path("psrio-psrfits-transpose.sf");
    file.write(path.path());

    SearchFile reader(path.path());
    const auto samples = reader.read_samples<float>(2);
    REQUIRE(samples.size() == 6U);
    CHECK(samples[0] == (200.0F - 256.0F) + 10.0F);
    CHECK(samples[1] == 0.0F);
    CHECK(samples[2] == 0.0F);
    CHECK(samples[5] == 7.0F);
}

TEST_CASE("SearchFile trims a short final sub-integration", "[psrfits]") {
    Geometry geo;
    geo.npol  = 1;
    geo.pol   = "INTEN";
    geo.nsblk = 4;
    geo.nsub  = 2;
    geo.tbin  = 1.0;
    const std::vector<double> freqs{1.0, 2.0};
    const std::vector<float> weights{1.0F, 1.0F};
    const std::vector<float> offsets{0.0F, 0.0F};
    const std::vector<float> scales{1.0F, 1.0F};
    std::vector<std::byte> table;
    push_row_prefix(table, geo, freqs, weights, offsets, scales, 4.0, 2.0);
    for (int i = 0; i < 8; ++i) {
        append_be(table, std::uint8_t{1});
    }
    push_row_prefix(table, geo, freqs, weights, offsets, scales, 2.0, 5.0);
    for (int i = 0; i < 8; ++i) {
        append_be(table, std::uint8_t{2});
    }
    FitsFileBytes file;
    primary_cards(file);
    subint_header(file, geo, table.size() / 2U);
    file.raw(table);
    file.pad_data();
    const TempPath path("psrio-psrfits-short.sf");
    file.write(path.path());

    SearchFile reader(path.path());
    CHECK(reader.nsamples() == 6U);
    const auto samples = reader.read_samples<float>(10);
    CHECK(samples.size() == 12U);
    CHECK(samples[8] == 2.0F);
    CHECK(samples[10] == 2.0F);
}

TEST_CASE("digifil sizes a search file from the stored sub-integration",
          "[psrfits][digifil]") {
    CHECK(psrio::digifil::psrfits_block_samples(8, 4, 6656, 8192, true, 1, 1) ==
          8192U);
}

TEST_CASE("SearchFile reads DATA and DAT_FREQ from the FITS heap",
          "[psrfits]") {
    FitsFileBytes file;
    primary_cards(file);
    file.begin_hdu();
    file.key("XTENSION", "'BINTABLE'");
    file.key("BITPIX", "8");
    file.key("NAXIS", "2");
    file.key("NAXIS1", "56");
    file.key("NAXIS2", "1");
    file.key("PCOUNT", "20");
    file.key("GCOUNT", "1");
    file.key("TFIELDS", "7");
    file.key("EXTNAME", "'SUBINT'");
    file.key("TTYPE1", "'TSUBINT'");
    file.key("TFORM1", "'1D'");
    file.key("TTYPE2", "'OFFS_SUB'");
    file.key("TFORM2", "'1D'");
    file.key("TTYPE3", "'DAT_FREQ'");
    file.key("TFORM3", "'1PD'");
    file.key("TTYPE4", "'DAT_WTS'");
    file.key("TFORM4", "'2E'");
    file.key("TTYPE5", "'DAT_OFFS'");
    file.key("TFORM5", "'2E'");
    file.key("TTYPE6", "'DAT_SCL'");
    file.key("TFORM6", "'2E'");
    file.key("TTYPE7", "'DATA'");
    file.key("TFORM7", "'1PB'");
    file.key("NPOL", "1");
    file.key("NCHAN", "2");
    file.key("NBITS", "8");
    file.key("NSBLK", "2");
    file.key("TBIN", "1.0");
    file.key("POL_TYPE", "'INTEN'");
    file.key("ZERO_OFF", "0.0");
    file.key("SIGNINT", "0");
    file.end_header();

    std::vector<std::byte> row;
    append_be(row, 2.0);
    append_be(row, 1.0);
    append_be(row, std::int32_t{2});
    append_be(row, std::int32_t{0});
    append_be(row, 1.0F);
    append_be(row, 1.0F);
    append_be(row, 0.0F);
    append_be(row, 0.0F);
    append_be(row, 1.0F);
    append_be(row, 1.0F);
    append_be(row, std::int32_t{4});
    append_be(row, std::int32_t{16});
    REQUIRE(row.size() == 56U);
    file.raw(row);
    std::vector<std::byte> heap;
    append_be(heap, 100.0);
    append_be(heap, 101.0);
    heap.push_back(std::byte{10});
    heap.push_back(std::byte{20});
    heap.push_back(std::byte{30});
    heap.push_back(std::byte{40});
    REQUIRE(heap.size() == 20U);
    file.raw(heap);
    file.pad_data();
    const TempPath path("psrio-psrfits-heap.sf");
    file.write(path.path());

    SearchFile reader(path.path());
    CHECK(reader.subint_samples() == 2U);
    CHECK(reader.input_npol() == 1U);
    CHECK(reader.fch1() == Catch::Approx(100.0));
    CHECK(reader.foff() == Catch::Approx(1.0));
    CHECK_FALSE(reader.sums_polarisations());
    const auto samples = reader.read_samples<float>(2);
    CHECK(samples == std::vector<float>{10.0F, 20.0F, 30.0F, 40.0F});
}

TEST_CASE("SearchFile rejects fold mode and a missing table", "[psrfits]") {
    {
        FitsFileBytes file;
        primary_cards(file, "PSR");
        file.begin_hdu();
        file.key("XTENSION", "'BINTABLE'");
        file.key("BITPIX", "8");
        file.key("NAXIS", "2");
        file.key("NAXIS1", "1");
        file.key("NAXIS2", "1");
        file.key("PCOUNT", "0");
        file.key("GCOUNT", "1");
        file.key("TFIELDS", "1");
        file.key("EXTNAME", "'SUBINT'");
        file.key("TTYPE1", "'DATA'");
        file.key("TFORM1", "'1B'");
        file.end_header();
        file.raw(std::vector<std::byte>{std::byte{0}});
        file.pad_data();
        const TempPath path("psrio-psrfits-fold.sf");
        file.write(path.path());
        CHECK_THROWS_WITH(SearchFile(path.path()), ContainsSubstring("SEARCH"));
    }
    {
        FitsFileBytes file;
        primary_cards(file);
        const TempPath path("psrio-psrfits-empty.sf");
        file.write(path.path());
        CHECK_THROWS_WITH(SearchFile(path.path()), ContainsSubstring("SUBINT"));
    }
}

TEST_CASE("Parkes UWL search file matches sigpyproc Stokes I",
          "[psrfits][uwl]") {
    const std::filesystem::path path =
        std::filesystem::path(PSRIO_SOURCE_DIR) /
        "tmp_files/test_data/uwl_210405_055326_60.sf";
    if (!std::filesystem::exists(path)) {
        SKIP("UWL PSRFITS file is not in tmp_files/test_data");
    }
    SearchFile reader(path);
    CHECK(reader.input_nbits() == 8);
    CHECK(reader.nchans() == 6656U);
    CHECK(reader.nifs() == 1U);
    CHECK(reader.pol_type() == "AABBCRCI");
    CHECK(reader.nsamples() == std::uint64_t{9} * 8192U);
    CHECK(reader.tsamp() == Catch::Approx(6.4e-5));
    CHECK(reader.fch1() == Catch::Approx(704.25));
    CHECK(reader.foff() == Catch::Approx(0.5));
    CHECK(reader.source_name() == "FRB20201124A");
    CHECK(reader.telescope() == "Parkes");
    CHECK(reader.beam() == 1);
    CHECK(reader.sigproc_header().refdm == Catch::Approx(412.0));
    const double tstart = 59309.0 + ((22999.0 + 0.169855999999901) / 86400.0);
    CHECK(reader.tstart() == Catch::Approx(tstart).epsilon(1.0e-12));

    const auto same_bits = [](float value, std::uint32_t bits) {
        return std::bit_cast<std::uint32_t>(value) == bits;
    };
    auto first = reader.read_samples<float>(32);
    REQUIRE(first.size() == std::size_t{32} * 6656U);
    CHECK(same_bits(first[0], 0x46b5e20c));
    CHECK(same_bits(first[1], 0x46b0c942));
    CHECK(same_bits(first[100], 0x479d814f));
    CHECK(same_bits(first[6655], 0x472baf84));
    CHECK(same_bits(first[6656], 0x46bcd85a));
    CHECK(same_bits(first[10UL * 6656UL], 0x46c468bc));
    std::uint32_t fold = 0;
    for (std::size_t index = 0; index < first.size(); index += 997U) {
        const auto bits = std::bit_cast<std::uint32_t>(first[index]);
        fold            = static_cast<std::uint32_t>(
            (static_cast<std::uint64_t>(fold) * 1315423911U) + bits);
    }
    CHECK(fold == 911571271U);

    reader.seek(8191);
    auto edge = reader.read_samples<float>(2);
    REQUIRE(edge.size() == std::size_t{2} * 6656U);
    CHECK(same_bits(edge[0], 0x46abbd9e));
    CHECK(same_bits(edge[6656], 0x46aa5978));
    CHECK(same_bits(edge[6657], 0x46a6b188));
}
