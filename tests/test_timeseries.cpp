#include "psrio/psrio.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cstddef>
#include <filesystem>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

[[nodiscard]] std::filesystem::path test_data_path(std::string_view filename) {
    return std::filesystem::path(PSRIO_TEST_DATA_DIR) / filename;
}

} // namespace

TEST_CASE("Header astronomical and timing conversions",
          "[timeseries][header]") {
    psrio::Header hdr;
    hdr.tsamp    = 0.0001; // 100 us
    hdr.nsamples = 100000; // 10 s
    hdr.tstart   = 55000.0;
    hdr.raj      = 120000.0;  // 12h 00m 00s
    hdr.dej      = -300000.0; // -30d 00m 00s
    hdr.fch1     = 1420.0;
    hdr.foff     = -1.0;
    hdr.nchans   = 1;

    REQUIRE_THAT(hdr.observation_duration(), WithinAbs(10.0, 1e-9));
    REQUIRE_THAT(hdr.tobs(), WithinAbs(10.0, 1e-9));
    REQUIRE(hdr.duration_string() == "10.0 seconds");
    REQUIRE(hdr.gregorian_date() == "2009-06-18");

    REQUIRE_THAT(hdr.ra_hours(), WithinAbs(12.0, 1e-6));
    REQUIRE_THAT(hdr.ra_degrees(), WithinAbs(180.0, 1e-6));
    REQUIRE_THAT(hdr.dec_degrees(), WithinAbs(-30.0, 1e-6));

    REQUIRE_THAT(hdr.fcenter(), WithinAbs(1420.0, 1e-6));
    REQUIRE_THAT(hdr.bandwidth(), WithinAbs(1.0, 1e-6));
}

TEST_CASE("TimeSeriesReader satisfies concepts::TimeSeriesReader concept",
          "[timeseries][concept]") {
    static_assert(psrio::concepts::TimeSeriesReader<psrio::TimeSeriesReader>);
}

TEST_CASE("TimeSeries loads GBT_J1807-0847.tim in-memory",
          "[timeseries][sigproc]") {
    const auto tim_file = test_data_path("GBT_J1807-0847.tim");
    REQUIRE(std::filesystem::exists(tim_file));

    const auto ts = psrio::TimeSeries::from_tim(tim_file);

    // Verify metadata matching sigpyproc3 ground truth
    REQUIRE(ts.nsamples() == 131130U);
    REQUIRE(ts.data().size() == 131130U);
    REQUIRE_THAT(ts.dt(), WithinAbs(0.00016384, 1e-10));
    REQUIRE_THAT(ts.header().tstart, WithinAbs(59313.30983797434, 1e-9));
    REQUIRE_THAT(ts.header().dm, WithinAbs(112.3802, 1e-4));
    REQUIRE(ts.header().source == "J1807-0847");
    REQUIRE(ts.header().telescope == "GBT");
    REQUIRE(ts.header().nbits == 32);
    REQUIRE(ts.header().nchans == 1);
    REQUIRE_THAT(ts.header().fch1, WithinAbs(920.0, 1e-3));
    REQUIRE_THAT(ts.tobs(), WithinRel(21.4843392, 1e-6));

    // Verify first 5 sample values against sigpyproc3 ground truth:
    // [-28.0, 110.0, -190.0, 92.0, 353.0]
    const auto samples = ts.data();
    REQUIRE_THAT(samples[0], WithinAbs(-28.0F, 1e-4F));
    REQUIRE_THAT(samples[1], WithinAbs(110.0F, 1e-4F));
    REQUIRE_THAT(samples[2], WithinAbs(-190.0F, 1e-4F));
    REQUIRE_THAT(samples[3], WithinAbs(92.0F, 1e-4F));
    REQUIRE_THAT(samples[4], WithinAbs(353.0F, 1e-4F));
}

TEST_CASE("TimeSeriesReader provides lite zero-allocation reading for .tim",
          "[timeseries][sigproc][lite]") {
    const auto tim_file = test_data_path("GBT_J1807-0847.tim");
    REQUIRE(std::filesystem::exists(tim_file));

    // 1. Lite open: parses metadata only, zero heap array allocation
    psrio::TimeSeriesReader reader(tim_file);
    REQUIRE(reader.nsamples() == 131130U);
    REQUIRE(reader.header().source == "J1807-0847");
    REQUIRE(reader.tell() == 0U);

    // 2. Caller-managed memory: unpacks directly into user's contiguous buffer
    std::vector<float> caller_buffer(reader.nsamples());
    const auto read_count = reader.read_data(caller_buffer);
    REQUIRE(read_count == 131130U);
    REQUIRE(reader.tell() == 131130U);

    // Verify samples match ground truth
    REQUIRE_THAT(caller_buffer[0], WithinAbs(-28.0F, 1e-4F));
    REQUIRE_THAT(caller_buffer[1], WithinAbs(110.0F, 1e-4F));
    REQUIRE_THAT(caller_buffer[2], WithinAbs(-190.0F, 1e-4F));
    REQUIRE_THAT(caller_buffer[3], WithinAbs(92.0F, 1e-4F));
    REQUIRE_THAT(caller_buffer[4], WithinAbs(353.0F, 1e-4F));

    // 3. Navigation: rewind and chunked streaming
    reader.rewind();
    REQUIRE(reader.tell() == 0U);

    std::array<float, 3> chunk{};
    const auto n_chunk = reader.read_samples(3, chunk);
    REQUIRE(n_chunk == 3U);
    REQUIRE(reader.tell() == 3U);
    REQUIRE_THAT(chunk[0], WithinAbs(-28.0F, 1e-4F));
    REQUIRE_THAT(chunk[1], WithinAbs(110.0F, 1e-4F));
    REQUIRE_THAT(chunk[2], WithinAbs(-190.0F, 1e-4F));

    // 4. Seeking
    reader.seek(4);
    REQUIRE(reader.tell() == 4U);
    std::array<float, 1> sample4{};
    REQUIRE(reader.read_samples(1, sample4) == 1U);
    REQUIRE_THAT(sample4[0], WithinAbs(353.0F, 1e-4F));
}

TEST_CASE("TimeSeries loads GBT_J1807-0847.dat and .inf in-memory",
          "[timeseries][presto]") {
    const auto dat_file = test_data_path("GBT_J1807-0847.dat");
    const auto inf_file = test_data_path("GBT_J1807-0847.inf");
    REQUIRE(std::filesystem::exists(dat_file));
    REQUIRE(std::filesystem::exists(inf_file));

    const auto ts = psrio::TimeSeries::from_dat(dat_file, inf_file);

    // Verify metadata matching sigpyproc3 ground truth
    REQUIRE(ts.nsamples() == 131072U);
    REQUIRE(ts.data().size() == 131072U);
    REQUIRE_THAT(ts.dt(), WithinAbs(0.00016384, 1e-10));
    REQUIRE_THAT(ts.header().tstart, WithinAbs(59313.30983797434, 1e-9));
    REQUIRE_THAT(ts.header().dm, WithinAbs(112.3802, 1e-4));
    REQUIRE(ts.header().source == "J1807-0847");
    REQUIRE(ts.header().telescope == "GBT");
    REQUIRE(ts.header().backend == "VEGAS");
    REQUIRE(ts.header().nbits == 32);
    REQUIRE_THAT(ts.tobs(), WithinRel(21.47483648, 1e-6));

    // Verify first 5 sample values against sigpyproc3 ground truth:
    // [444259.0, 445709.0, 445747.0, 446036.0, 446026.0]
    const auto samples = ts.data();
    REQUIRE_THAT(samples[0], WithinAbs(444259.0F, 1e-2F));
    REQUIRE_THAT(samples[1], WithinAbs(445709.0F, 1e-2F));
    REQUIRE_THAT(samples[2], WithinAbs(445747.0F, 1e-2F));
    REQUIRE_THAT(samples[3], WithinAbs(446036.0F, 1e-2F));
    REQUIRE_THAT(samples[4], WithinAbs(446026.0F, 1e-2F));
}

TEST_CASE(
    "TimeSeriesReader provides lite zero-allocation reading for PRESTO .dat",
    "[timeseries][presto][lite]") {
    const auto dat_file = test_data_path("GBT_J1807-0847.dat");
    const auto inf_file = test_data_path("GBT_J1807-0847.inf");

    // Open via companion inference or explicit paths
    psrio::TimeSeriesReader reader(dat_file, inf_file);
    REQUIRE(reader.nsamples() == 131072U);
    REQUIRE(reader.header().backend == "VEGAS");

    // Read directly into caller's buffer
    std::vector<float> caller_buffer(reader.nsamples());
    const auto count = reader.read_data(caller_buffer);
    REQUIRE(count == 131072U);
    REQUIRE_THAT(caller_buffer[0], WithinAbs(444259.0F, 1e-2F));
    REQUIRE_THAT(caller_buffer[4], WithinAbs(446026.0F, 1e-2F));

    // Zero copy view of payload
    reader.rewind();
    const auto view = reader.view_bytes(sizeof(float) * 2);
    REQUIRE(view.size() == 8U);
}

TEST_CASE("TimeSeries to_tim and to_dat round-trips correctly",
          "[timeseries][roundtrip]") {
    psrio::Header hdr;
    hdr.source    = "J1807-0847";
    hdr.telescope = "GBT";
    hdr.backend   = "VEGAS";
    hdr.tstart    = 59313.309837;
    hdr.tsamp     = 0.00016384;
    hdr.dm        = 112.3802;
    hdr.fch1      = 920.0;
    hdr.foff      = 0.5;
    hdr.ra        = "18:07:38.0000";
    hdr.dec       = "-08:47:43.7500";
    hdr.raj       = 180738.0;
    hdr.dej       = -84743.75;

    const std::vector<float> sample_data{1.0F, 2.5F, -3.2F, 4.0F, 10.5F};
    const psrio::TimeSeries ts(sample_data, hdr);

    // 1. Round-trip .tim
    const auto tim_tmp =
        std::filesystem::temp_directory_path() / "psrio_test_roundtrip.tim";
    ts.to_tim(tim_tmp);
    REQUIRE(std::filesystem::exists(tim_tmp));

    const auto ts_tim_back = psrio::TimeSeries::from_tim(tim_tmp);
    REQUIRE(ts_tim_back.nsamples() == 5U);
    REQUIRE(ts_tim_back.header().source == "J1807-0847");
    REQUIRE_THAT(ts_tim_back.header().dm, WithinAbs(112.3802, 1e-4));
    for (std::size_t i = 0; i < sample_data.size(); ++i) {
        REQUIRE_THAT(ts_tim_back.data()[i], WithinAbs(sample_data[i], 1e-5F));
    }

    // 2. Round-trip PRESTO .dat / .inf
    const auto presto_base =
        std::filesystem::temp_directory_path() / "psrio_test_roundtrip_presto";
    ts.to_dat(presto_base);
    const auto dat_tmp = presto_base.string() + ".dat";
    const auto inf_tmp = presto_base.string() + ".inf";
    REQUIRE(std::filesystem::exists(dat_tmp));
    REQUIRE(std::filesystem::exists(inf_tmp));

    const auto ts_dat_back = psrio::TimeSeries::from_dat(dat_tmp, inf_tmp);
    REQUIRE(ts_dat_back.nsamples() == 5U);
    REQUIRE(ts_dat_back.header().source == "J1807-0847");
    REQUIRE_THAT(ts_dat_back.header().dm, WithinAbs(112.3802, 1e-4));
    for (std::size_t i = 0; i < sample_data.size(); ++i) {
        REQUIRE_THAT(ts_dat_back.data()[i], WithinAbs(sample_data[i], 1e-5F));
    }

    std::error_code ec;
    std::filesystem::remove(tim_tmp, ec);
    std::filesystem::remove(dat_tmp, ec);
    std::filesystem::remove(inf_tmp, ec);
}

TEST_CASE("TimeSeries and Header auto-detect format via from_file and open",
          "[timeseries][autodetect]") {
    const auto tim_file = test_data_path("GBT_J1807-0847.tim");
    const auto dat_file = test_data_path("GBT_J1807-0847.dat");
    const auto inf_file = test_data_path("GBT_J1807-0847.inf");

    // Header::from_file
    const auto hdr_tim = psrio::Header::from_file(tim_file);
    REQUIRE(hdr_tim.source == "J1807-0847");
    REQUIRE(hdr_tim.nsamples == 131130U);

    const auto hdr_dat = psrio::Header::from_file(dat_file);
    REQUIRE(hdr_dat.source == "J1807-0847");
    REQUIRE(hdr_dat.nsamples == 131072U);

    const auto hdr_inf = psrio::Header::from_file(inf_file);
    REQUIRE(hdr_inf.source == "J1807-0847");
    REQUIRE(hdr_inf.nsamples == 131072U);

    // TimeSeries::from_file
    const auto ts_tim = psrio::TimeSeries::from_file(tim_file);
    REQUIRE(ts_tim.nsamples() == 131130U);
    REQUIRE(ts_tim.header().source == "J1807-0847");

    const auto ts_dat = psrio::TimeSeries::from_file(dat_file);
    REQUIRE(ts_dat.nsamples() == 131072U);
    REQUIRE(ts_dat.header().source == "J1807-0847");

    // TimeSeries::open
    const auto reader = psrio::TimeSeries::open(tim_file);
    REQUIRE(reader.nsamples() == 131130U);
}

TEST_CASE("TimeSeriesReader error validation and boundary handling",
          "[timeseries][validation]") {
    const auto tim_file = test_data_path("GBT_J1807-0847.tim");
    psrio::TimeSeriesReader reader(tim_file);

    // Seeking past end throws
    REQUIRE_THROWS_AS(reader.seek(reader.nsamples() + 1),
                      psrio::ValidationError);

    // read_data with mismatched buffer size throws
    std::vector<float> small_buf(10);
    REQUIRE_THROWS_AS(reader.read_data(small_buf), psrio::ValidationError);

    // read_samples with mismatched destination count throws
    std::vector<float> sample_buf(5);
    REQUIRE_THROWS_AS(reader.read_samples(10, sample_buf),
                      psrio::ValidationError);

    // Missing companion inf file throws IoError
    const auto non_existent_inf =
        std::filesystem::temp_directory_path() / "non_existent.inf";
    const auto non_existent_dat =
        std::filesystem::temp_directory_path() / "non_existent.dat";
    REQUIRE_THROWS_AS(
        psrio::TimeSeriesReader(non_existent_dat, non_existent_inf),
        psrio::IoError);
}

TEST_CASE("Header::make_inf produces valid parseable PRESTO metadata",
          "[timeseries][make_inf]") {
    const auto inf_file = test_data_path("GBT_J1807-0847.inf");
    const auto orig_hdr = psrio::Header::from_inffile(inf_file);

    const auto tmp_inf =
        std::filesystem::temp_directory_path() / "test_make_inf.inf";
    orig_hdr.make_inf(tmp_inf);
    REQUIRE(std::filesystem::exists(tmp_inf));

    const auto parsed_hdr = psrio::Header::from_inffile(tmp_inf);
    REQUIRE(parsed_hdr.source == orig_hdr.source);
    REQUIRE(parsed_hdr.telescope == orig_hdr.telescope);
    REQUIRE(parsed_hdr.backend == orig_hdr.backend);
    REQUIRE(parsed_hdr.nsamples == orig_hdr.nsamples);
    REQUIRE_THAT(parsed_hdr.tsamp, WithinAbs(orig_hdr.tsamp, 1e-12));
    REQUIRE_THAT(parsed_hdr.dm, WithinAbs(orig_hdr.dm, 1e-5));
    REQUIRE_THAT(parsed_hdr.tstart, WithinAbs(orig_hdr.tstart, 1e-9));

    std::error_code ec;
    std::filesystem::remove(tmp_inf, ec);
}
