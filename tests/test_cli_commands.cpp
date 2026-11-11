#include "cli.hpp"
#include "psrfits_bytes.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <dmt/algorithms/ddmt.hpp>
#include <dmt/algorithms/fdmt.hpp>
#include <filesystem>
#include <fmt/format.h>
#include <initializer_list>
#include <psrio/block_source.hpp>
#include <psrio/digifil.hpp>
#include <psrio/formats/sigproc.hpp>
#include <psrio/timeseries.hpp>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {

namespace fs = std::filesystem;
using psrio::TimeSeriesReader;
using psrio::formats::sigproc::FilterbankReader;

const fs::path kData{PSRIO_TEST_DATA_DIR};
const fs::path kTutorial = kData / "tutorial.fil";
const fs::path kTim      = kData / "GBT_J1807-0847.tim";

class TempDir {
public:
    TempDir() {
        static std::atomic<int> counter{0};
        m_path = fs::temp_directory_path() /
                 ("psrio_cli_" + std::to_string(::getpid()) + "_" +
                  std::to_string(counter.fetch_add(1)));
        fs::remove_all(m_path);
        fs::create_directories(m_path);
    }
    ~TempDir() { fs::remove_all(m_path); }

    TempDir(const TempDir&)            = delete;
    TempDir& operator=(const TempDir&) = delete;
    TempDir(TempDir&&)                 = delete;
    TempDir& operator=(TempDir&&)      = delete;

    [[nodiscard]] std::string operator/(const std::string& name) const {
        return (m_path / name).string();
    }
    [[nodiscard]] const fs::path& path() const noexcept { return m_path; }

private:
    fs::path m_path;
};

struct Result {
    int code;
    std::string out;
    std::string err;
};

Result psrio_run(std::initializer_list<std::string> args) {
    const std::vector<std::string> argv(args);
    std::ostringstream out;
    std::ostringstream err;
    const int code = psrio::cli::run(argv, out, err);
    return {.code = code, .out = out.str(), .err = err.str()};
}

void require_ok(const Result& result) {
    INFO(result.err);
    REQUIRE(result.code == 0);
}

std::vector<float> read_all(FilterbankReader& reader) {
    reader.rewind();
    return reader.read_samples<float>(reader.nsamples());
}

std::vector<float> read_series(const fs::path& path) {
    TimeSeriesReader reader(path);
    return reader.read_samples(reader.nsamples());
}

std::vector<std::byte> read_payload(const fs::path& path) {
    FilterbankReader reader(path);
    return reader.read_block(reader.nsamples());
}

/// Values of channels [c0, c0 + nc) of every sample.
std::vector<float> channels(std::span<const float> data,
                            std::size_t nchans,
                            std::size_t c0,
                            std::size_t nc) {
    std::vector<float> out;
    for (std::size_t base = 0; base < data.size(); base += nchans) {
        out.insert(out.end(),
                   data.begin() + static_cast<std::ptrdiff_t>(base + c0),
                   data.begin() + static_cast<std::ptrdiff_t>(base + c0 + nc));
    }
    return out;
}

} // namespace

TEST_CASE("psrio without a command prints usage and fails", "[cli]") {
    const auto result = psrio_run({});
    CHECK(result.code != 0);
    const auto help = psrio_run({"--help"});
    CHECK(help.code == 0);
    CHECK(help.out.find("dedisperse") != std::string::npos);
}

TEST_CASE("psrio inspect prints keys, JSON, and statistics", "[cli][inspect]") {
    auto result = psrio_run({"inspect", kTutorial.string(), "-k", "nchans"});
    require_ok(result);
    CHECK(result.out == "64\n");

    result = psrio_run({"header", kTutorial.string(), "--key", "tsamp"});
    require_ok(result);
    CHECK(std::stod(result.out) == 0.00032);

    result = psrio_run({"inspect", kTutorial.string(), "-k", "source"});
    require_ok(result);
    CHECK(result.out.find("P: 250") == 0U);

    result = psrio_run({"inspect", kTutorial.string(), "--json", "--stats"});
    require_ok(result);
    CHECK(result.out.front() == '{');
    CHECK(result.out.find("\"nchans\": 64") != std::string::npos);
    CHECK(result.out.find("\"nan_fraction\": 0") != std::string::npos);

    result = psrio_run({"inspect", kTim.string(), "-k", "nsamples"});
    require_ok(result);
    CHECK(std::stoull(result.out) == TimeSeriesReader(kTim).nsamples());

    result = psrio_run(
        {"inspect", (kData / "GBT_J1807-0847.dat").string(), "-k", "dm"});
    require_ok(result);
    CHECK(std::stod(result.out) == Catch::Approx(112.3802));

    result = psrio_run({"inspect", kTutorial.string(), "-k", "nope"});
    CHECK(result.code == 1);
    CHECK(result.err.find("available keys") != std::string::npos);
}

TEST_CASE("psrio slice copies a time and channel range", "[cli][slice]") {
    const TempDir dir;
    FilterbankReader in(kTutorial);
    const auto all = read_all(in);

    SECTION("byte-aligned channels are copied raw") {
        require_ok(psrio_run({"slice", kTutorial.string(), "-o", dir / "s.fil",
                              "-s", "1000", "-n", "5000", "--start-chan", "8",
                              "--nchans", "16"}));
        FilterbankReader out(dir / "s.fil");
        CHECK(out.nsamples() == 5000U);
        CHECK(out.nchans() == 16U);
        CHECK(out.header().fch1.value() == Catch::Approx(1510.0 - (8 * 1.09)));
        CHECK(out.header().tstart.value() ==
              Catch::Approx(in.header().mjd_after_samples(1000)));
        const auto expected =
            channels(std::span<const float>(all).subspan(
                         std::size_t{1000} * 64, std::size_t{5000} * 64),
                     64, 8, 16);
        CHECK(read_all(out) == expected);
    }

    SECTION("unaligned sub-byte channels are repacked") {
        require_ok(psrio_run({"extract", kTutorial.string(), "-o",
                              dir / "u.fil", "--start-chan", "3", "--nchans",
                              "8", "--duration", "1.0"}));
        FilterbankReader out(dir / "u.fil");
        CHECK(out.nsamples() == 3125U);
        const auto expected =
            channels(std::span<const float>(all).first(std::size_t{3125} * 64),
                     64, 3, 8);
        CHECK(read_all(out) == expected);
    }

    SECTION("a whole-byte stride is required") {
        const auto result =
            psrio_run({"slice", kTutorial.string(), "-o", dir / "bad.fil",
                       "--start-chan", "3", "--nchans", "5"});
        CHECK(result.code == 1);
    }

    SECTION("time series slice to PRESTO") {
        require_ok(psrio_run({"slice", kTim.string(), "-o", dir / "t.dat",
                              "--start-sample", "10", "-n", "100"}));
        const auto full = read_series(kTim);
        const auto part = read_series(dir / "t.dat");
        CHECK(part ==
              std::vector<float>(full.begin() + 10, full.begin() + 110));
    }

    SECTION("the output may not overwrite the input") {
        const auto result =
            psrio_run({"slice", kTutorial.string(), "-o", kTutorial.string()});
        CHECK(result.code == 1);
    }
}

TEST_CASE("psrio splice joins sub-bands", "[cli][splice]") {
    const TempDir dir;
    require_ok(psrio_run({"slice", kTutorial.string(), "-o", dir / "lo.fil",
                          "--start-chan", "32", "-n", "4000"}));
    require_ok(psrio_run({"slice", kTutorial.string(), "-o", dir / "hi.fil",
                          "--nchans", "32", "-n", "4000"}));
    require_ok(psrio_run(
        {"slice", kTutorial.string(), "-o", dir / "ref.fil", "-n", "4000"}));

    SECTION("contiguous bands sort by descending frequency") {
        require_ok(psrio_run({"splice", dir / "lo.fil", dir / "hi.fil", "-o",
                              dir / "sp.fil", "--gulp", "999"}));
        const FilterbankReader out(dir / "sp.fil");
        CHECK(out.header().frequency_table.empty());
        CHECK(out.header().fch1.value() == 1510.0);
        CHECK(out.header().foff.value() == Catch::Approx(-1.09));
        CHECK(read_payload(dir / "sp.fil") == read_payload(dir / "ref.fil"));
    }

    SECTION("--keep-order writes a frequency table") {
        require_ok(psrio_run({"splice", dir / "lo.fil", dir / "hi.fil", "-o",
                              dir / "kp.fil", "--keep-order"}));
        const FilterbankReader out(dir / "kp.fil");
        REQUIRE(out.header().frequency_table.size() == 64U);
        CHECK(out.header().frequency_table.front() ==
              Catch::Approx(1510.0 - (32 * 1.09)));
    }

    SECTION("mismatched sample times are rejected") {
        require_ok(psrio_run({"decimate", dir / "lo.fil", "-o",
                              dir / "slow.fil", "-t", "2", "-b", "8"}));
        CHECK(psrio_run({"splice", dir / "slow.fil", dir / "hi.fil", "-o",
                         dir / "x.fil"})
                  .code == 1);
    }
}

TEST_CASE("psrio decimate averages samples and channels", "[cli][decimate]") {
    const TempDir dir;
    FilterbankReader in(kTutorial);
    const auto all = read_all(in);

    require_ok(psrio_run({"decimate", kTutorial.string(), "-o", dir / "d.fil",
                          "-t", "4", "-c", "8", "--gulp", "1001"}));
    FilterbankReader out(dir / "d.fil");
    CHECK(out.nbits() == 32);
    CHECK(out.nchans() == 8U);
    CHECK(out.nsamples() == 187520U / 4U);
    CHECK(out.header().tsamp.value() == Catch::Approx(0.00128));
    CHECK(out.header().foff.value() == Catch::Approx(-8.72));
    CHECK(out.header().fch1.value() == Catch::Approx(1510.0 - (3.5 * 1.09)));
    const auto got = out.read_samples<float>(3);
    for (std::size_t s = 0; s < 3; ++s) {
        for (std::size_t c = 0; c < 8; ++c) {
            float acc = 0.0F;
            for (std::size_t t = 0; t < 4; ++t) {
                for (std::size_t k = 0; k < 8; ++k) {
                    acc += all[(((s * 4) + t) * 64) + (c * 8) + k];
                }
            }
            CHECK(got[(s * 8) + c] == Catch::Approx(acc / 32.0F));
        }
    }

    SECTION("one remaining channel writes a time series") {
        require_ok(psrio_run({"decimate", kTutorial.string(), "-o",
                              dir / "d.tim", "-t", "3", "-c", "64", "--sum"}));
        const auto series = read_series(dir / "d.tim");
        CHECK(series.size() == 187520U / 3U);
        float first = 0.0F;
        for (std::size_t i = 0; i < std::size_t{3} * 64; ++i) {
            first += all[i];
        }
        CHECK(series.front() == first);
    }

    SECTION("a factor that does not divide nchans is rejected") {
        CHECK(psrio_run({"decimate", kTutorial.string(), "-o", dir / "x.fil",
                         "-c", "3"})
                  .code == 1);
    }
}

TEST_CASE("psrio convert changes bit depth, channel order, and format",
          "[cli][convert]") {
    const TempDir dir;
    FilterbankReader in(kTutorial);
    const auto all = read_all(in);

    SECTION("--fswap reverses the channels losslessly") {
        require_ok(psrio_run(
            {"convert", kTutorial.string(), "-o", dir / "f.fil", "--fswap"}));
        FilterbankReader out(dir / "f.fil");
        CHECK(out.nbits() == 2);
        CHECK(out.header().fch1.value() == Catch::Approx(1510.0 - (63 * 1.09)));
        CHECK(out.header().foff.value() == Catch::Approx(1.09));
        out.set_fswap(true);
        CHECK(read_all(out) == all);
    }

    SECTION("more bits keep the values") {
        require_ok(psrio_run(
            {"convert", kTutorial.string(), "-o", dir / "b8.fil", "-b", "8"}));
        FilterbankReader out(dir / "b8.fil");
        CHECK(out.nbits() == 8);
        CHECK(read_all(out) == all);
    }

    SECTION("fewer bits go through digifil") {
        require_ok(psrio_run({"convert", kTutorial.string(), "-o",
                              dir / "b16.fil", "-b", "16"}));
        require_ok(psrio_run(
            {"convert", dir / "b16.fil", "-o", dir / "b4.fil", "-b", "4"}));
        const FilterbankReader out(dir / "b4.fil");
        CHECK(out.nbits() == 4);
        CHECK(out.nsamples() == 187520U);
        CHECK(psrio_run({"convert", dir / "b16.fil", "-o", dir / "x.fil", "-b",
                         "8", "--fswap"})
                  .code == 1);
    }

    SECTION(".tim and .dat round trip") {
        require_ok(psrio_run({"convert", kTim.string(), "-o", dir / "a.dat"}));
        require_ok(psrio_run({"convert", dir / "a.dat", "-o", dir / "b.tim"}));
        const auto original = read_series(kTim);
        CHECK(read_series(dir / "a.dat") == original);
        CHECK(read_series(dir / "b.tim") == original);
        CHECK(psrio_run(
                  {"convert", kTim.string(), "-o", dir / "c.dat", "-b", "8"})
                  .code == 1);
    }
}

TEST_CASE("psrio digifil matches the library pipeline", "[cli][digifil]") {
    const TempDir dir;
    require_ok(psrio_run({"digifil", kTutorial.string(), "-o", dir / "d.fil",
                          "-b", "2", "-t", "3", "-f", "4", "-I", "0.5"}));

    FilterbankReader reader(kTutorial);
    const auto header = reader.header();
    psrio::BlockSource source(std::move(reader));
    psrio::digifil::Config config;
    config.nbits           = 2;
    config.tscrunch        = 3;
    config.fscrunch        = 4;
    config.rescale_seconds = 0.5;
    psrio::digifil::run(source, header, dir / "lib.fil", config);
    CHECK(read_payload(dir / "d.fil") == read_payload(dir / "lib.fil"));
    CHECK(FilterbankReader(dir / "d.fil").nchans() == 16U);
}

TEST_CASE("psrio dedisperse writes one series per DM", "[cli][dedisperse]") {
    const TempDir dir;
    FilterbankReader in(kTutorial);
    const auto all      = read_all(in);
    const auto nsamples = in.nsamples();
    const auto nchans   = in.nchans();

    std::vector<float> channel_sum(static_cast<std::size_t>(nsamples), 0.0F);
    for (std::size_t s = 0; s < nsamples; ++s) {
        for (std::size_t c = 0; c < nchans; ++c) {
            channel_sum[s] += all[(s * nchans) + c];
        }
    }

    const auto f_min = static_cast<float>(1510.0 - (63 * 1.09) - 0.545);
    const auto f_max = static_cast<float>(1510.0 + 0.545);
    const std::vector<float> dm10{10.0F};
    const dmt::algorithms::DDMT reference(f_min, f_max, nchans, 0.00032F, dm10,
                                          {}, 32, {}, 1, f_max);
    const auto nout10 = reference.get_output_nsamps(nsamples);

    SECTION("single DMs to .tim and .dat") {
        require_ok(psrio_run({"dedisperse", kTutorial.string(), "-o",
                              dir / "dm0.tim", "--dm", "0", "--gulp", "5000"}));
        CHECK(read_series(dir / "dm0.tim") == channel_sum);

        require_ok(psrio_run({"dedisperse", kTutorial.string(), "-o",
                              dir / "dm10.dat", "-d", "10", "-j", "2"}));
        TimeSeriesReader series(dir / "dm10.dat");
        CHECK(series.nsamples() == nout10);
        CHECK(series.header().dm == Catch::Approx(10.0));
        CHECK(series.header().tstart == 50000.0);
        CHECK(fs::exists(dir / "dm10.inf"));

        const auto& table = reference.get_plan().get_container().delay_table;
        const auto got    = series.read_samples(16);
        for (std::size_t k = 0; k < got.size(); ++k) {
            float acc = 0.0F;
            for (std::size_t c = 0; c < nchans; ++c) {
                const auto file_chan = nchans - 1U - c;
                acc += all[((k + table[c]) * nchans) + file_chan];
            }
            CHECK(got[k] == acc);
        }
    }

    SECTION("a DM list writes prefixed files that agree across engines") {
        require_ok(psrio_run({"dedisperse", kTutorial.string(), "-o",
                              (dir.path() / "run").string(), "--dm", "0,10"}));
        require_ok(psrio_run({"dedisperse", kTutorial.string(), "-o",
                              (dir.path() / "sd").string(), "--dm-start", "0",
                              "--dm-end", "10", "--dm-step", "10", "--algo",
                              "sdmt", "--gulp", "777"}));
        for (const auto* dm : {"0.00", "10.00"}) {
            const auto ddmt =
                read_series(dir / ("run_DM" + std::string(dm) + ".tim"));
            const auto sdmt =
                read_series(dir / ("sd_DM" + std::string(dm) + ".tim"));
            CHECK(ddmt.size() == nout10);
            CHECK(ddmt == sdmt);
        }
        CHECK(read_series(dir / "run_DM0.00.tim") ==
              std::vector<float>(channel_sum.begin(),
                                 channel_sum.begin() +
                                     static_cast<std::ptrdiff_t>(nout10)));
    }

    SECTION("fdmt at DM 0 is the channel sum") {
        require_ok(psrio_run({"dedisperse", kTutorial.string(), "-o",
                              (dir.path() / "f").string(), "-d", "0,10",
                              "--algo", "fdmt", "--gulp", "4096"}));
        const auto dm0 = read_series(dir / "f_DM0.00.tim");
        REQUIRE(dm0.size() > nsamples - 64U);
        REQUIRE(dm0.size() < nsamples);
        CHECK(dm0 ==
              std::vector<float>(channel_sum.begin(),
                                 channel_sum.begin() +
                                     static_cast<std::ptrdiff_t>(dm0.size())));
    }

    SECTION("fdmt output does not depend on the gulp") {
        for (const auto* fref : {"top", "bottom"}) {
            for (const auto* gulp : {"1000", "4096", "262144"}) {
                require_ok(psrio_run({"dedisperse", kTutorial.string(), "-o",
                                      (dir.path() / gulp).string(), "-d",
                                      "0,30", "--algo", "fdmt", "--fref", fref,
                                      "--gulp", gulp}));
            }
            const auto whole = read_series(dir / "262144_DM30.16.tim");
            REQUIRE(whole.size() > 100000U);
            CHECK(read_series(dir / "1000_DM30.16.tim") == whole);
            CHECK(read_series(dir / "4096_DM30.16.tim") == whole);
        }
    }

    SECTION("fdmt places an impulse on its own trace at the reference time") {
        const std::vector<float> grid{0.0F, 30.0F};
        const dmt::algorithms::FDMT plan_engine(
            f_min, f_max, nchans, 4096, 0.00032F, grid, false, "valid", {});
        const auto& plan    = plan_engine.get_plan();
        const auto trial    = plan.get_dm_grid_final().size() - 1U;
        const auto trace    = plan.trace_dm(trial);
        const auto lowest   = std::ranges::min(trace);
        std::int64_t dt_max = 0;
        for (std::size_t d = 0; d < plan.get_dm_grid_final().size(); ++d) {
            dt_max = std::max<std::int64_t>(
                dt_max, -std::ranges::min(plan.trace_dm(d)));
        }

        constexpr std::int64_t kOut = 3000;
        auto header                 = in.header();
        header.nbits                = 32;
        std::vector<float> data(static_cast<std::size_t>(nsamples * nchans),
                                0.0F);
        for (std::size_t c = 0; c < nchans; ++c) {
            const auto sample = static_cast<std::size_t>(kOut + trace[c]);
            data[(sample * nchans) + (nchans - 1U - c)] = 1.0F;
        }
        {
            psrio::formats::sigproc::FilterbankWriter writer(dir / "imp.fil",
                                                             header);
            writer.write_samples(nsamples, std::span<const float>(data));
            writer.close();
        }

        const auto peak = [](const std::vector<float>& series) {
            const auto it = std::ranges::max_element(series);
            return std::pair{static_cast<std::int64_t>(it - series.begin()),
                             *it};
        };
        const auto name =
            fmt::format("_DM{:.2f}.tim", plan.get_dm_grid_final()[trial]);
        for (const auto* fref : {"top", "bottom"}) {
            const auto prefix = (dir.path() / fref).string();
            require_ok(
                psrio_run({"dedisperse", dir / "imp.fil", "-o", prefix, "-d",
                           "0,30", "--algo", "fdmt", "--fref", fref}));
            const auto [index, value] = peak(read_series(prefix + name));
            CHECK(value == static_cast<float>(nchans));
            // Top: the top channel's sample. Bottom: the output time at the
            // lower band edge, less the sweep folded into tstart.
            const auto expected =
                std::string_view(fref) == "top" ? kOut + lowest : kOut - dt_max;
            CHECK(index == expected);
        }
    }

    SECTION("masked channels are left out of the sum") {
        require_ok(
            psrio_run({"dedisperse", kTutorial.string(), "-o", dir / "m.tim",
                       "-d", "0", "--ignorechans", "0,2-3"}));
        const auto masked = read_series(dir / "m.tim");
        for (std::size_t s = 0; s < 8; ++s) {
            const auto* row = &all[s * nchans];
            CHECK(masked[s] == channel_sum[s] - row[0] - row[2] - row[3]);
        }
    }

    SECTION("bad arguments are rejected") {
        CHECK(psrio_run({"dedisperse", kTutorial.string(), "-o", dir / "x.tim"})
                  .code == 1);
        CHECK(psrio_run({"dedisperse", kTutorial.string(), "-o", dir / "x.tim",
                         "-d", "1", "--dm-start", "0"})
                  .code != 0);
        CHECK(psrio_run({"dedisperse", kTutorial.string(), "-o", dir / "x.tim",
                         "-d", "1", "--ignorechans", "64"})
                  .code == 1);
    }
}

TEST_CASE("psrio inspect reads a search-mode PSRFITS header",
          "[cli][psrfits]") {
    using psrio::test::append_be;
    using psrio::test::FitsFileBytes;
    FitsFileBytes file;
    file.begin_hdu();
    file.key("SIMPLE", "T");
    file.key("BITPIX", "8");
    file.key("NAXIS", "0");
    file.key("EXTEND", "T");
    file.key("FITSTYPE", "'PSRFITS'");
    file.key("OBS_MODE", "'SEARCH'");
    file.key("TELESCOP", "'Parkes'");
    file.key("SRC_NAME", "'CLITest'");
    file.key("STT_IMJD", "50000");
    file.key("STT_SMJD", "0");
    file.key("STT_OFFS", "0.0");
    file.end_header();
    file.begin_hdu();
    file.key("XTENSION", "'BINTABLE'");
    file.key("BITPIX", "8");
    file.key("NAXIS", "2");
    file.key("NAXIS1", "58");
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
    file.key("TFORM7", "'2B'");
    file.key("NPOL", "1");
    file.key("NCHAN", "2");
    file.key("NBITS", "8");
    file.key("NSBLK", "1");
    file.key("TBIN", "0.001");
    file.key("POL_TYPE", "'INTEN'");
    file.key("ZERO_OFF", "0.0");
    file.key("SIGNINT", "0");
    file.end_header();
    std::vector<std::byte> row;
    append_be(row, 0.001);
    append_be(row, 0.0005);
    append_be(row, 1400.0);
    append_be(row, 1399.0);
    append_be(row, 1.0F);
    append_be(row, 1.0F);
    append_be(row, 0.0F);
    append_be(row, 0.0F);
    append_be(row, 1.0F);
    append_be(row, 1.0F);
    append_be(row, std::uint8_t{3});
    append_be(row, std::uint8_t{4});
    REQUIRE(row.size() == 58U);
    file.raw(row);
    file.pad_data();

    const TempDir dir;
    const auto fits = dir.path() / "tiny.sf";
    file.write(fits);
    const auto source =
        psrio_run({"inspect", fits.string(), "-k", "source_name"});
    require_ok(source);
    CHECK(source.out == "CLITest\n");
    const auto bits = psrio_run({"inspect", fits.string(), "-k", "nbits"});
    require_ok(bits);
    CHECK(bits.out == "32\n");
    const auto stored =
        psrio_run({"inspect", fits.string(), "-k", "stored_nbits"});
    require_ok(stored);
    CHECK(stored.out == "8\n");
    const auto format = psrio_run({"inspect", fits.string(), "-k", "format"});
    require_ok(format);
    CHECK(format.out == "PSRFITS\n");

    const auto sliced = dir.path() / "out.fil";
    require_ok(
        psrio_run({"slice", fits.string(), "-o", sliced.string(), "-n", "1"}));
    psrio::formats::sigproc::FilterbankReader reader(sliced);
    CHECK(reader.nchans() == 2U);
    CHECK(reader.nbits() == 32);
    const auto samples = reader.read_samples<float>(1);
    CHECK(samples[0] == 3.0F);
    CHECK(samples[1] == 4.0F);
}

TEST_CASE("psrio joins contiguous filterbanks in time", "[cli][stream]") {
    const TempDir dir;
    const auto header = [](double tstart) {
        psrio::formats::sigproc::FilterbankHeader out;
        out.nchans      = 2;
        out.nifs        = 1;
        out.nbits       = 8;
        out.tsamp       = 1.0;
        out.tstart      = tstart;
        out.fch1        = 1400.0;
        out.foff        = -1.0;
        out.source_name = "stream";
        return out;
    };
    const auto write = [](const fs::path& path,
                          const psrio::formats::sigproc::FilterbankHeader& hdr,
                          std::initializer_list<std::uint8_t> samples) {
        psrio::formats::sigproc::FilterbankWriter writer(path, hdr);
        std::vector<std::byte> bytes(samples.size());
        std::size_t index = 0;
        for (const auto sample : samples) {
            bytes[index++] = static_cast<std::byte>(sample);
        }
        writer.write_block(samples.size() / 2U, bytes);
        writer.close();
    };
    const auto first  = dir.path() / "a.fil";
    const auto second = dir.path() / "b.fil";
    write(first, header(50000.0), {1, 2, 3, 4});
    write(second, header(50000.0 + (2.0 / 86400.0)), {5, 6, 7, 8});

    const auto count = psrio_run(
        {"inspect", first.string(), second.string(), "-k", "nsamples"});
    require_ok(count);
    CHECK(count.out == "4\n");
    const auto files =
        psrio_run({"inspect", first.string(), second.string(), "-k", "nfiles"});
    require_ok(files);
    CHECK(files.out == "2\n");

    const auto sliced = dir.path() / "joined.fil";
    require_ok(psrio_run(
        {"slice", first.string(), second.string(), "-o", sliced.string()}));
    FilterbankReader out(sliced);
    CHECK(out.nsamples() == 4U);
    const auto got = out.read_samples<std::uint8_t>(4);
    CHECK(got == std::vector<std::uint8_t>{1, 2, 3, 4, 5, 6, 7, 8});

    const auto gapped = dir.path() / "gap.fil";
    write(gapped, header(50000.0 + (4.0 / 86400.0)), {9, 10, 11, 12});
    CHECK(psrio_run(
              {"inspect", first.string(), gapped.string(), "-k", "nsamples"})
              .code == 1);
}
