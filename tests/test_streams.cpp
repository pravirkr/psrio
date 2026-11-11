#include "psrfits_bytes.hpp"
#include "psrio/block_source.hpp"
#include "psrio/formats/psrfits.hpp"
#include "psrio/formats/sigproc.hpp"
#include "psrio/open.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <unistd.h>
#include <utility>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using psrio::formats::psrfits::SearchStream;
using psrio::formats::sigproc::FilterbankHeader;
using psrio::formats::sigproc::FilterbankStream;
using psrio::formats::sigproc::FilterbankWriter;
using psrio::test::append_be;
using psrio::test::FitsFileBytes;

namespace {

class TempPath {
public:
    explicit TempPath(std::string_view name) {
        m_path = std::filesystem::temp_directory_path() /
                 ("psrio_stream_" + std::to_string(::getpid()) + "_" +
                  std::string(name));
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

FilterbankHeader fil_header(double tstart, std::string source = "stream") {
    FilterbankHeader header;
    header.nchans      = 2;
    header.nifs        = 1;
    header.nbits       = 8;
    header.tsamp       = 1.0;
    header.tstart      = tstart;
    header.fch1        = 1400.0;
    header.foff        = -1.0;
    header.source_name = std::move(source);
    return header;
}

void write_fil(const std::filesystem::path& path,
               double tstart,
               std::span<const std::uint8_t> samples,
               std::string source             = "stream",
               std::optional<std::string> raw = std::nullopt) {
    auto header        = fil_header(tstart, std::move(source));
    header.rawdatafile = std::move(raw);
    FilterbankWriter writer(path, header);
    std::vector<std::byte> bytes(samples.size());
    for (std::size_t i = 0; i < samples.size(); ++i) {
        bytes[i] = static_cast<std::byte>(samples[i]);
    }
    writer.write_block(samples.size() / 2U, bytes);
    writer.close();
}

double mjd_after(double tstart, std::uint64_t samples, double tsamp) {
    return tstart + (static_cast<double>(samples) * tsamp / 86400.0);
}

void write_search(const std::filesystem::path& path,
                  double stt_offs,
                  std::optional<std::int64_t> nsuboffs,
                  std::string_view source,
                  std::span<const std::uint8_t> data) {
    FitsFileBytes file;
    file.begin_hdu();
    file.key("SIMPLE", "T");
    file.key("BITPIX", "8");
    file.key("NAXIS", "0");
    file.key("EXTEND", "T");
    file.key("FITSTYPE", "'PSRFITS'");
    file.key("OBS_MODE", "'SEARCH'");
    file.key("TELESCOP", "'Parkes'");
    file.key("SRC_NAME", std::string("'") + std::string(source) + "'");
    file.key("STT_IMJD", "50000");
    file.key("STT_SMJD", "0");
    file.key("STT_OFFS", std::to_string(stt_offs));
    file.key("CHAN_DM", "12.5");
    file.key("IBEAM", "1");
    file.end_header();

    file.begin_hdu();
    file.key("XTENSION", "'BINTABLE'");
    file.key("BITPIX", "8");
    file.key("NAXIS", "2");
    file.key("NAXIS1", "60");
    file.key("NAXIS2", "1");
    file.key("PCOUNT", "0");
    file.key("GCOUNT", "1");
    file.key("TFIELDS", "7");
    file.key("EXTNAME", "'SUBINT'");
    file.key("TTYPE1", "'TSUBINT'");
    file.key("TFORM1", "'1D'");
    file.key("TTYPE2", "'OFFS_SUB'");
    file.key("TFORM2", "'1D'");
    file.key("TTYPE3", "'DAT_FREQ'");
    file.key("TFORM3", "'2D'");
    file.key("TTYPE4", "'DAT_WTS'");
    file.key("TFORM4", "'2E'");
    file.key("TTYPE5", "'DAT_OFFS'");
    file.key("TFORM5", "'2E'");
    file.key("TTYPE6", "'DAT_SCL'");
    file.key("TFORM6", "'2E'");
    file.key("TTYPE7", "'DATA'");
    file.key("TFORM7", "'4B'");
    file.key("NPOL", "1");
    file.key("NCHAN", "2");
    file.key("NBITS", "8");
    file.key("NSBLK", "2");
    file.key("TBIN", "1.0");
    file.key("POL_TYPE", "'INTEN'");
    file.key("ZERO_OFF", "0.0");
    file.key("SIGNINT", "0");
    if (nsuboffs.has_value()) {
        file.key("NSUBOFFS", std::to_string(*nsuboffs));
    }
    file.end_header();

    std::vector<std::byte> row;
    append_be(row, 2.0);
    append_be(row, 1.0);
    append_be(row, 1400.0);
    append_be(row, 1399.0);
    append_be(row, 1.0F);
    append_be(row, 1.0F);
    append_be(row, 0.0F);
    append_be(row, 0.0F);
    append_be(row, 1.0F);
    append_be(row, 1.0F);
    for (const auto value : data) {
        append_be(row, value);
    }
    REQUIRE(row.size() == 60U);
    REQUIRE(data.size() == 4U);
    file.raw(row);
    file.pad_data();
    file.write(path);
}

} // namespace

TEST_CASE("FilterbankStream joins contiguous SIGPROC files",
          "[sigproc][stream]") {
    static_assert(psrio::concepts::BlockReader<FilterbankStream>);

    const TempPath first("a.fil");
    const TempPath second("b.fil");
    constexpr std::array<std::uint8_t, 4> kLeft  = {1, 2, 3, 4};
    constexpr std::array<std::uint8_t, 4> kRight = {5, 6, 7, 8};
    write_fil(first.path(), 50000.0, kLeft);
    write_fil(second.path(), mjd_after(50000.0, 2, 1.0), kRight, "stream",
              "raw.dat");

    const std::array<std::filesystem::path, 2> paths = {first.path(),
                                                        second.path()};
    FilterbankStream stream(paths);
    CHECK(stream.nsamples() == 4U);
    CHECK(stream.header().nsamples() == 4U);
    CHECK(stream.tstart() == 50000.0);
    CHECK(stream.source_name() == "stream");

    std::vector<std::uint8_t> got(6);
    CHECK(stream.read_samples(3, std::span<std::uint8_t>(got)) == 3U);
    CHECK(got == std::vector<std::uint8_t>{1, 2, 3, 4, 5, 6});
    CHECK(stream.tell() == 3U);

    stream.seek(2);
    stream.set_fswap(true);
    std::vector<std::uint8_t> swapped(4);
    CHECK(stream.read_samples(2, std::span<std::uint8_t>(swapped)) == 2U);
    CHECK(swapped == std::vector<std::uint8_t>{6, 5, 8, 7});

    stream.seek(1);
    stream.set_fswap(false);
    std::vector<std::byte> raw(6);
    CHECK(stream.read_bytes(6, raw) == 6U);
    CHECK(raw[0] == std::byte{3});
    CHECK(raw[1] == std::byte{4});
    CHECK(raw[2] == std::byte{5});
    CHECK(raw[4] == std::byte{7});

    psrio::BlockSource block =
        psrio::open_block_source(std::span<const std::filesystem::path>(paths));
    CHECK(block.nsamples() == 4U);
    std::vector<float> floats(2);
    block.seek(3);
    CHECK(block.read_samples(1, std::span<float>(floats)) == 1U);
    CHECK(floats[0] == 7.0F);
    CHECK(floats[1] == 8.0F);

    stream.rewind();
    CHECK(stream.read_samples<std::uint8_t>(2) ==
          std::vector<std::uint8_t>{1, 2, 3, 4});
}

TEST_CASE("FilterbankStream rejects a gap or a different header",
          "[sigproc][stream]") {
    const TempPath first("a.fil");
    const TempPath second("b.fil");
    constexpr std::array<std::uint8_t, 4> kSamples = {1, 2, 3, 4};
    write_fil(first.path(), 50000.0, kSamples);
    write_fil(second.path(), mjd_after(50000.0, 3, 1.0), kSamples);
    const std::array<std::filesystem::path, 2> gapped = {first.path(),
                                                         second.path()};
    CHECK_THROWS_WITH(FilterbankStream(gapped),
                      ContainsSubstring("not contiguous"));

    write_fil(second.path(), mjd_after(50000.0, 2, 1.0), kSamples, "other");
    const std::array<std::filesystem::path, 2> renamed = {first.path(),
                                                          second.path()};
    CHECK_THROWS_WITH(FilterbankStream(renamed),
                      ContainsSubstring("source_name"));

    write_fil(second.path(), mjd_after(50000.0, 2, 1.0), kSamples);
    const std::array<std::filesystem::path, 2> reversed = {second.path(),
                                                           first.path()};
    CHECK_THROWS_WITH(FilterbankStream(reversed),
                      ContainsSubstring("not contiguous"));
}

TEST_CASE("SearchStream reads across the file boundary", "[psrfits][stream]") {
    static_assert(psrio::concepts::BlockReader<SearchStream>);

    const TempPath first("a.sf");
    const TempPath second("b.sf");
    constexpr std::array<std::uint8_t, 4> kLeft  = {10, 20, 30, 40};
    constexpr std::array<std::uint8_t, 4> kRight = {50, 60, 70, 80};
    write_search(first.path(), 0.0, 0, "Join", kLeft);
    write_search(second.path(), 2.0, 1, "Join", kRight);
    const std::array<std::filesystem::path, 2> paths = {first.path(),
                                                        second.path()};

    SearchStream stream(paths);
    std::vector<float> got(6);
    REQUIRE(stream.read_samples(3, std::span<float>(got)) == 3U);
    CHECK(got == std::vector<float>{10.0F, 20.0F, 30.0F, 40.0F, 50.0F, 60.0F});

    stream.seek(2);
    stream.set_fswap(true);
    std::vector<float> swapped(4);
    REQUIRE(stream.read_samples(2, std::span<float>(swapped)) == 2U);
    CHECK(swapped == std::vector<float>{60.0F, 50.0F, 80.0F, 70.0F});

    const psrio::BlockSource block =
        psrio::open_block_source(std::span<const std::filesystem::path>(paths));
    CHECK(block.nsamples() == 4U);

    stream.rewind();
    stream.set_fswap(false);
    CHECK(stream.read_samples<float>(2) ==
          std::vector<float>{10.0F, 20.0F, 30.0F, 40.0F});
}

TEST_CASE("SearchStream rejects a gap, a header change, and a bad NSUBOFFS",
          "[psrfits][stream]") {
    const TempPath first("a.sf");
    const TempPath second("b.sf");
    constexpr std::array<std::uint8_t, 4> kData = {1, 2, 3, 4};
    write_search(first.path(), 0.0, 0, "Join", kData);
    write_search(second.path(), 3.0, 1, "Join", kData);
    const std::array<std::filesystem::path, 2> gapped = {first.path(),
                                                         second.path()};
    CHECK_THROWS_WITH(SearchStream(gapped),
                      ContainsSubstring("not contiguous"));

    write_search(second.path(), 2.0, 5, "Join", kData);
    const std::array<std::filesystem::path, 2> shifted = {first.path(),
                                                          second.path()};
    CHECK_THROWS_WITH(SearchStream(shifted), ContainsSubstring("NSUBOFFS"));

    write_search(second.path(), 2.0, std::nullopt, "Join", kData);
    const std::array<std::filesystem::path, 2> partial = {first.path(),
                                                          second.path()};
    CHECK_THROWS_WITH(SearchStream(partial), ContainsSubstring("NSUBOFFS"));

    write_search(first.path(), 0.0, std::nullopt, "Join", kData);
    write_search(second.path(), 2.0, std::nullopt, "Other", kData);
    const std::array<std::filesystem::path, 2> renamed = {first.path(),
                                                          second.path()};
    CHECK_THROWS_WITH(SearchStream(renamed),
                      ContainsSubstring("observing setup"));

    write_search(second.path(), 2.0, std::nullopt, "Join", kData);
    SearchStream open{std::span<const std::filesystem::path>(renamed)};
    CHECK(open.nsamples() == 4U);
}
