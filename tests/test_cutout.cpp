#include "cutout_fixture.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/mjd.hpp"
#include "psrio/detail/rfile.hpp"
#include "psrio/formats/psrfits/archive.hpp"
#include "psrio/formats/psrfits/extract.hpp"
#include "psrio/formats/psrfits/fits_writer.hpp"
#include "psrio/formats/psrfits/reader.hpp"
#include "psrio/formats/psrfits/stream.hpp"
#include "psrio/formats/psrfits/table.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <span>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

namespace fs = std::filesystem;
namespace pf = psrio::formats::psrfits;
using psrio::detail::Mjd;
using psrio::test::BurstFixture;

class TempDir {
public:
    TempDir() {
        static std::atomic<int> counter{0};
        m_path = fs::temp_directory_path() /
                 ("psrio_cutout_" + std::to_string(::getpid()) + "_" +
                  std::to_string(counter.fetch_add(1)));
        fs::remove_all(m_path);
        fs::create_directories(m_path);
    }
    ~TempDir() { fs::remove_all(m_path); }
    TempDir(const TempDir&)            = delete;
    TempDir& operator=(const TempDir&) = delete;
    TempDir(TempDir&&)                 = delete;
    TempDir& operator=(TempDir&&)      = delete;
    [[nodiscard]] fs::path operator/(const std::string& name) const {
        return m_path / name;
    }

private:
    fs::path m_path;
};

/// Stream start of the fixture, as dspsr builds it from STT_*.
Mjd fixture_start() { return {59630, 3600, 0.5}; }

pf::DspsrBand fixture_band(const BurstFixture& fixture, bool dedisperse) {
    return pf::make_dspsr_band(BurstFixture::kObsFreq, BurstFixture::kObsBw,
                               BurstFixture::kNchan, 1.0 / BurstFixture::kTbin,
                               fixture.dm, dedisperse);
}

/// fma(raw - 127.5, 2, 1), the fixture's calibrated value.
float calibrated(std::uint8_t raw) {
    return std::fma(static_cast<float>(raw) - 127.5F, 2.0F, 1.0F);
}

std::vector<float> fold(pf::SearchStream& stream,
                        const pf::ExtractWindow& window,
                        const pf::DspsrBand& band) {
    std::vector<float> out(static_cast<std::size_t>(
        stream.input_npol() * stream.nchans() * window.nbin));
    // A small chunk so the window spans several reads.
    pf::fold_window(stream, window, band, out, 4096);
    return out;
}

} // namespace

TEST_CASE("dspsr sample delays and band", "[cutout]") {
    // Parkes UWL: dspsr -K advanced STT by 3128 samples (0.200192 s).
    const auto uwl =
        pf::make_dspsr_band(2368.0, 3328.0, 6656, 1.0 / 6.4e-05, 413.0, true);
    CHECK(uwl.frequencies.front() == 704.25);
    CHECK(uwl.frequencies.back() == 4031.75);
    CHECK(uwl.span_centre == 0);
    CHECK(uwl.zero_delay == 3128);
    CHECK(uwl.applied.back() == 0);
    CHECK(uwl.min_applied() == 0);
    CHECK(uwl.applied.front() ==
          3128 - pf::dspsr_sample_delay(413.0, 2368.0, 704.25, 15625.0));

    const auto flat =
        pf::make_dspsr_band(2368.0, 3328.0, 6656, 1.0 / 6.4e-05, 413.0, false);
    CHECK(flat.zero_delay == 0);
    CHECK(flat.max_applied() == 0);
    CHECK(flat.max_abs_delay ==
          -pf::dspsr_sample_delay(413.0, 2368.0, 704.25, 15625.0));

    const double offset = pf::dispersion_offset(413.0, 4031.75, 2368.0);
    CHECK(offset == Catch::Approx(0.2001859).margin(1e-6));
    CHECK(std::llround(offset * 15625.0) == 3128);
}

TEST_CASE("Mjd and window timing reproduce dspsr's stored times",
          "[cutout][mjd]") {
    // Parkes UWL file uwl_220220_084024_122.sf, dspsr -S 5.24288 -c 0.065536
    // -b 1024 -s: STT_SMJD 34875, STT_OFFS 0.246540000000194, and
    // OFFS_SUB of the first sub-integration as below.
    const Mjd start(59630, 34870, 0.00366000000019373);
    const double rate = 1.0 / 6.4e-05;
    const auto flat =
        pf::make_dspsr_band(2368.0, 3328.0, 6656, rate, 413.0, false);
    const pf::ExtractWindow window{
        .first = 81920, .nsamples = 1024, .nbin = 1024};
    const auto timing = pf::window_timing(start, 8192, window, flat);
    CHECK(timing.start.day() == 59630);
    CHECK(timing.start.second() == 34875);
    CHECK(pf::Card::format_double(timing.start.fraction()) ==
          "0.246540000000194");
    CHECK(timing.offs_sub == -0.018380552232550973);
    CHECK(timing.period == 0.065536);

    // The same window with -K: STT moves by 3128 samples. (dspsr writes
    // ...194 with -U 16384 and ...193 with -U 8192: its block layout changes
    // the last rounding. OFFS_SUB is the same either way.)
    const auto band =
        pf::make_dspsr_band(2368.0, 3328.0, 6656, rate, 413.0, true);
    const auto dedispersed = pf::window_timing(start, 8192, window, band);
    CHECK(pf::Card::format_double(dedispersed.start.fraction()) ==
          "0.446732000000194");
    CHECK(dedispersed.offs_sub == -0.021964050035285343);

    const auto parsed = pf::parse_mjd("59630.403700000000001");
    CHECK(parsed.day() == 59630);
    CHECK(pf::format_mjd(parsed) == "59630.403700000000001");
    CHECK(pf::format_mjd(Mjd(59630, 0, 0.0)) == "59630.000000000000000");
    CHECK_THROWS_AS(pf::parse_mjd("abc"), psrio::ValidationError);

    // Seconds carry into days and back.
    const Mjd before(59630, 86399, 0.75);
    const auto late = before + 0.5;
    CHECK(late.day() == 59631);
    CHECK(late.second() == 0);
    CHECK(late.fraction() == 0.25);
    CHECK(late.seconds_since(before) == 0.5);
    CHECK(before.seconds_since(late) == -0.5);
    CHECK((late - 0.5) == before);

    // A negative interval keeps every part non-positive.
    const Mjd back(0, 0, -1.25);
    CHECK(back.second() == -1);
    CHECK(back.fraction() == -0.25);
    CHECK(Mjd(0, 2, -0.25).fraction() == 0.75);

    // Midpoints, including an odd number of whole seconds.
    const auto mid = Mjd::midpoint(Mjd(59630, 10, 0.5), Mjd(59630, 13, 0.0));
    CHECK(mid.day() == 59630);
    CHECK(mid.second() == 11);
    CHECK(mid.fraction() == 0.75);
}

TEST_CASE("quantise_profile matches psrchive's int16 rule", "[cutout]") {
    const std::array<float, 4> amps{-1.0F, 0.0F, 1.0F, 2.0F};
    std::array<std::int16_t, 4> data{};
    const auto scale  = pf::quantise_profile(amps, data);
    const double low  = -32751.0;
    const double high = 32750.0;
    CHECK(scale.offset ==
          static_cast<float>(((-1.0 * high) - (2.0 * low)) / (high - low)));
    CHECK(scale.scale == static_cast<float>(3.0 / (high - low)));
    CHECK(data[0] == -32751);
    CHECK(data[3] == 32750);
    for (std::size_t bin = 0; bin < amps.size(); ++bin) {
        const float back =
            (static_cast<float>(data[bin]) * scale.scale) + scale.offset;
        CHECK(std::abs(back - amps[bin]) <= scale.scale);
    }

    // A flat profile stores scale 0, the value as offset, and zeros.
    const std::array<float, 3> flat{5.0F, 5.0F, 5.0F};
    std::array<std::int16_t, 3> zeros{1, 1, 1};
    const auto flat_scale = pf::quantise_profile(flat, zeros);
    CHECK(flat_scale.scale == 0.0F);
    CHECK(flat_scale.offset == 5.0F);
    CHECK(zeros == std::array<std::int16_t, 3>{0, 0, 0});

    std::array<std::int16_t, 2> wrong{};
    CHECK_THROWS_AS(pf::quantise_profile(amps, wrong), psrio::ValidationError);
}

TEST_CASE("FITS cards follow CFITSIO formatting", "[cutout][fits]") {
    CHECK(pf::Card::format_double(2368.0) == "2368.");
    CHECK(pf::Card::format_double(6.4e-05) == "6.4E-05");
    CHECK(pf::Card::format_double(1e-05) == "1.E-05");
    CHECK(pf::Card::format_double(-4554231.6) == "-4554231.6");
    const auto real = pf::Card::real("OBSFREQ", 2368.0, "[MHz] centre").image();
    CHECK(real.size() == 80U);
    CHECK(real.starts_with("OBSFREQ =                2368. / [MHz] centre"));
    const auto text = pf::Card::string("SRC_NAME", "a'b").image();
    CHECK(text.starts_with("SRC_NAME= 'a''b    '"));
    CHECK(pf::Card::integer("NBIN", 1024).image().substr(10, 20) ==
          "                1024");
}

TEST_CASE("read_products unpacks every product as dspsr", "[cutout][psrfits]") {
    const TempDir dir;
    const BurstFixture fixture;
    const auto path = dir / "burst.sf";
    fixture.write(path, 0, fixture.nrows);
    const std::vector<fs::path> paths{path};

    pf::SearchStream calibrated_stream(
        paths,
        pf::DecodeOptions{.scloffs = true, .weights = true, .product = 1});
    CHECK(calibrated_stream.product_values() == 32U);
    calibrated_stream.seek(3);
    std::vector<float> one(32);
    REQUIRE(calibrated_stream.read_products(1, one) == 1U);
    for (int pol = 0; pol < 4; ++pol) {
        for (int chan = 0; chan < 8; ++chan) {
            // Weights (0 for channel 0) and the product choice are ignored.
            CHECK(one[static_cast<std::size_t>((pol * 8) + chan)] ==
                  calibrated(fixture.raw(3, pol, chan)));
        }
    }

    pf::SearchStream raw_stream(paths, pf::DecodeOptions{.scloffs = false});
    raw_stream.seek(3);
    REQUIRE(raw_stream.read_products(1, one) == 1U);
    CHECK(one[0] == 100.0F);
    CHECK(one[16] == 50.0F);
    std::vector<float> wrong(5);
    CHECK_THROWS_AS(raw_stream.read_products(1, wrong), psrio::ValidationError);
}

TEST_CASE("fold_window dedisperses and centres the burst", "[cutout]") {
    const TempDir dir;
    const BurstFixture fixture;
    const auto path = dir / "burst.sf";
    fixture.write(path, 0, fixture.nrows);
    const std::vector<fs::path> paths{path};
    pf::SearchStream stream(paths, pf::DecodeOptions{.scloffs = true});
    const auto band = fixture_band(fixture, true);
    REQUIRE(band.zero_delay > 0);

    // TOA of the pulse at the band centre; dspsr's f0 channel of output o
    // reads input o + span + zero.
    const auto toa = fixture_start() +
                     (static_cast<double>(fixture.pulse) * BurstFixture::kTbin);
    const auto window = pf::centred_window(
        fixture_start(), toa, BurstFixture::kObsFreq, 128, 32, band);
    CHECK(window.first + band.span_centre + band.zero_delay + 64 ==
          fixture.pulse);

    // The same burst named by its arrival at the top channel.
    const double top = band.frequencies.back();
    const auto toa_top =
        toa - pf::dispersion_offset(fixture.dm, top, BurstFixture::kObsFreq);
    const auto from_top =
        pf::centred_window(fixture_start(), toa_top, top, 128, 32, band);
    CHECK(from_top.first == window.first);

    const auto amps  = fold(stream, window, band);
    const float base = calibrated(BurstFixture::kBaseline);
    const float peak =
        (calibrated(BurstFixture::kPulse) + (3.0F * base)) / 4.0F;
    for (int pol = 0; pol < 2; ++pol) {
        for (int chan = 0; chan < 8; ++chan) {
            const auto* prof =
                amps.data() +
                (static_cast<std::size_t>((pol * 8) + chan) * 32U);
            for (int bin = 0; bin < 32; ++bin) {
                INFO("pol " << pol << " chan " << chan << " bin " << bin);
                CHECK(prof[bin] == (bin == 16 ? peak : base));
            }
        }
    }
    for (int chan = 0; chan < 8; ++chan) {
        const auto* prof =
            amps.data() + (static_cast<std::size_t>(16 + chan) * 32U);
        CHECK(prof[0] == calibrated(static_cast<std::uint8_t>(
                             BurstFixture::kCross + chan)));
    }
}

TEST_CASE("fold_window keeps the sweep without -K", "[cutout]") {
    const TempDir dir;
    const BurstFixture fixture;
    const auto path = dir / "burst.sf";
    fixture.write(path, 0, fixture.nrows);
    const std::vector<fs::path> paths{path};
    pf::SearchStream stream(paths, pf::DecodeOptions{.scloffs = false});
    const auto band = fixture_band(fixture, false);
    const auto toa = fixture_start() +
                     (static_cast<double>(fixture.pulse) * BurstFixture::kTbin);
    const auto window = pf::centred_window(
        fixture_start(), toa, BurstFixture::kObsFreq, 256, 256, band);
    CHECK(window.first == fixture.pulse - 128);
    const auto amps   = fold(stream, window, band);
    const auto delays = fixture.delays();
    for (int chan = 0; chan < 8; ++chan) {
        const auto* prof =
            amps.data() + (static_cast<std::size_t>(chan) * 256U);
        const auto peak =
            std::distance(prof, std::max_element(prof, prof + 256));
        CHECK(peak == 128 - delays[static_cast<std::size_t>(chan)]);
        CHECK(prof[peak] == 200.0F);
    }
}

TEST_CASE("fold_window reads across contiguous files", "[cutout][stream]") {
    const TempDir dir;
    const BurstFixture fixture;
    fixture.write(dir / "whole.sf", 0, fixture.nrows);
    fixture.write(dir / "a.sf", 0, 8);
    fixture.write(dir / "b.sf", 8, 8);
    const std::vector<fs::path> whole{dir / "whole.sf"};
    const std::vector<fs::path> split{dir / "a.sf", dir / "b.sf"};
    pf::SearchStream one(whole, pf::DecodeOptions{.scloffs = true});
    pf::SearchStream two(split, pf::DecodeOptions{.scloffs = true});
    const auto band = fixture_band(fixture, true);
    // Output samples 450..578 read input across the boundary at 512.
    const pf::ExtractWindow window{.first = 450, .nsamples = 128, .nbin = 64};
    CHECK(fold(one, window, band) == fold(two, window, band));
}

TEST_CASE("fold_window rejects windows outside the data", "[cutout]") {
    const TempDir dir;
    const BurstFixture fixture;
    fixture.write(dir / "burst.sf", 0, fixture.nrows);
    const std::vector<fs::path> paths{dir / "burst.sf"};
    pf::SearchStream stream(paths, pf::DecodeOptions{});
    const auto band = fixture_band(fixture, true);
    std::vector<float> out(std::size_t{4} * 8U * 16U);
    CHECK_THROWS_AS(pf::fold_window(stream,
                                    {.first = 1000, .nsamples = 32, .nbin = 16},
                                    band, out),
                    psrio::ValidationError);
    CHECK_THROWS_AS(pf::fold_window(stream,
                                    {.first = 0, .nsamples = 33, .nbin = 16},
                                    band, out),
                    psrio::ValidationError);
    CHECK_THROWS_AS(pf::centred_window(fixture_start(), fixture_start(),
                                       BurstFixture::kObsFreq, 64, 16, band),
                    psrio::ValidationError);
}

TEST_CASE("write_fold_archive writes a readable fold-mode PSRFITS file",
          "[cutout][fits]") {
    const TempDir dir;
    const BurstFixture fixture;
    fixture.write(dir / "burst.sf", 0, fixture.nrows);
    const std::vector<fs::path> paths{dir / "burst.sf"};
    pf::SearchStream stream(paths, pf::DecodeOptions{.scloffs = true});
    const auto band = fixture_band(fixture, true);
    const pf::ExtractWindow window{.first = 400, .nsamples = 64, .nbin = 64};
    const auto amps   = fold(stream, window, band);
    const auto timing = pf::window_timing(fixture_start(), 64, window, band);

    pf::FoldArchiveSpec spec;
    spec.input_primary    = stream.file(0).table().primary;
    spec.input_subint     = stream.file(0).table().subint;
    spec.start            = timing.start;
    spec.period           = timing.period;
    spec.tsubint          = timing.tsubint;
    spec.offs_sub         = timing.offs_sub;
    spec.nbin             = 64;
    spec.nchan            = 8;
    spec.npol             = 4;
    spec.pol_type         = pf::archive_pol_type(stream.pol_type(), 4);
    spec.dm               = fixture.dm;
    spec.dedispersed      = true;
    spec.centre_frequency = band.centre_frequency;
    spec.bandwidth        = band.bandwidth;
    spec.frequencies      = band.frequencies;
    spec.command          = "psrio cutout test";
    spec.created          = 0;
    const auto out        = dir / "out.ar";
    pf::write_fold_archive(out, spec, amps);
    CHECK(fs::file_size(out) % 2880U == 0U);

    const psrio::detail::ReadFile file(out);
    const pf::detail::Source source(file);
    const auto table = pf::detail::parse_search(source, out.string());
    CHECK(table.primary.value("OBS_MODE") == "PSR");
    CHECK(table.primary.value("COORD_MD") == "J2000");
    CHECK(table.primary.value("RA") == "05:08:03.508");
    CHECK(table.primary.value("DEC") == "+26:03:38.504");
    CHECK(table.primary.value("IBEAM") == "1");
    CHECK(table.primary.value("SCANLEN") == "*");
    CHECK(table.primary.integer("STT_IMJD") == timing.start.day());
    CHECK(table.primary.integer("STT_SMJD") == timing.start.second());
    CHECK(table.subint.value("POL_TYPE") == "AABBCRCI");
    CHECK(table.subint.integer("NBIN") == 64);
    CHECK(table.subint.number("TBIN") == BurstFixture::kTbin);
    CHECK(table.subint.number("ZERO_OFF") == 127.5);
    CHECK(table.subint.value("NSUBOFFS") == "*");
    CHECK(table.nrows == 1U);

    std::vector<std::byte> cell;
    std::vector<double> freqs(8);
    pf::detail::load_reals(
        pf::detail::read_cell(source, table, 0,
                              pf::detail::require_column(table, "DAT_FREQ"),
                              cell),
        pf::detail::require_column(table, "DAT_FREQ"), freqs);
    CHECK(freqs == band.frequencies);
    std::vector<double> weights(8);
    pf::detail::load_reals(
        pf::detail::read_cell(source, table, 0,
                              pf::detail::require_column(table, "DAT_WTS"),
                              cell),
        pf::detail::require_column(table, "DAT_WTS"), weights);
    CHECK(std::ranges::all_of(weights, [](double w) { return w == 1.0; }));

    std::vector<double> scales(32);
    std::vector<double> offsets(32);
    pf::detail::load_reals(
        pf::detail::read_cell(source, table, 0,
                              pf::detail::require_column(table, "DAT_SCL"),
                              cell),
        pf::detail::require_column(table, "DAT_SCL"), scales);
    pf::detail::load_reals(
        pf::detail::read_cell(source, table, 0,
                              pf::detail::require_column(table, "DAT_OFFS"),
                              cell),
        pf::detail::require_column(table, "DAT_OFFS"), offsets);
    const auto& data = pf::detail::require_column(table, "DATA");
    CHECK(data.repeat == 4U * 8U * 64U);
    std::vector<double> stored(static_cast<std::size_t>(data.repeat));
    pf::detail::load_reals(pf::detail::read_cell(source, table, 0, data, cell),
                           data, stored);
    for (std::size_t prof = 0; prof < 32U; ++prof) {
        for (std::size_t bin = 0; bin < 64U; ++bin) {
            const auto value =
                (stored[(prof * 64U) + bin] * scales[prof]) + offsets[prof];
            const auto want = static_cast<double>(amps[(prof * 64U) + bin]);
            CHECK(std::abs(value - want) <= std::max(scales[prof], 1e-6));
        }
    }
}

TEST_CASE("write_fold_archive clamps HISTORY NBIN as psrchive does",
          "[cutout][fits]") {
    const TempDir dir;
    const BurstFixture fixture;
    fixture.write(dir / "burst.sf", 0, 2);
    const std::vector<fs::path> paths{dir / "burst.sf"};
    const pf::SearchStream stream(paths, pf::DecodeOptions{});
    constexpr std::uint32_t kBins = 32768;
    pf::FoldArchiveSpec spec;
    spec.input_primary    = stream.file(0).table().primary;
    spec.input_subint     = stream.file(0).table().subint;
    spec.start            = fixture_start();
    spec.period           = 1.0;
    spec.tsubint          = 1.0;
    spec.nbin             = kBins;
    spec.nchan            = 8;
    spec.npol             = 4;
    spec.pol_type         = "AABBCRCI";
    spec.dm               = 1.0;
    spec.centre_frequency = BurstFixture::kObsFreq;
    spec.bandwidth        = BurstFixture::kObsBw;
    spec.frequencies      = pf::dspsr_channel_frequencies(
        BurstFixture::kObsFreq, BurstFixture::kObsBw, BurstFixture::kNchan);
    std::vector<float> amps(std::size_t{4} * 8U * kBins, 1.0F);
    amps[7]        = 2.0F;
    const auto out = dir / "wide.ar";
    pf::write_fold_archive(out, spec, amps);

    const psrio::detail::ReadFile file(out);
    const pf::detail::Source source(file);
    std::uint64_t history_at = 0;
    (void)pf::detail::parse_header(source, 0, history_at);
    std::uint64_t row_at = 0;
    const auto history   = pf::detail::parse_header(source, history_at, row_at);
    REQUIRE(history.value("EXTNAME") == "HISTORY");
    std::vector<std::byte> row(
        static_cast<std::size_t>(history.integer("NAXIS1")));
    source.read(row_at, row);
    // DATE_PRO 24, PROC_CMD 256, SCALE 8, POL_TYPE 8, NSUB 4, NPOL 2, NBIN 2.
    constexpr std::size_t kScale = 24 + 256;
    constexpr std::size_t kNbin  = kScale + 8 + 8 + 4 + 2;
    CHECK(pf::detail::load_be<std::int16_t>(row.data() + kNbin) == 32767);
    CHECK(std::all_of(row.begin() + kScale, row.begin() + kScale + 8,
                      [](std::byte b) { return b == std::byte{0}; }));
    CHECK(std::all_of(row.begin() + kNbin + 2, row.end(),
                      [](std::byte b) { return b == std::byte{0}; }));

    const auto table = pf::detail::parse_search(source, out.string());
    CHECK(table.subint.integer("NBIN") == kBins);
}
