#include "psrio/detail/concepts.hpp"
#include "psrio/formats/presto.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <span>
#include <string_view>
#include <system_error>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;

namespace {

class TempPrestoFile {
public:
    TempPrestoFile(std::string_view inf_text,
                   std::span<const float> dat_samples) {
        m_base = std::filesystem::temp_directory_path() / "psrio_test_presto";
        m_dat_path = m_base.string() + ".dat";
        m_inf_path = m_base.string() + ".inf";

        {
            std::ofstream inf_out(m_inf_path);
            inf_out << inf_text;
        }
        {
            std::ofstream dat_out(m_dat_path, std::ios::binary);
            dat_out.write(reinterpret_cast<const char*>(dat_samples.data()),
                          static_cast<std::streamsize>(dat_samples.size() *
                                                       sizeof(float)));
        }
    }

    TempPrestoFile(const TempPrestoFile&)            = delete;
    TempPrestoFile& operator=(const TempPrestoFile&) = delete;
    TempPrestoFile(TempPrestoFile&&)                 = delete;
    TempPrestoFile& operator=(TempPrestoFile&&)      = delete;

    ~TempPrestoFile() {
        std::error_code ec;
        std::filesystem::remove(m_dat_path, ec);
        std::filesystem::remove(m_inf_path, ec);
    }

    [[nodiscard]] const std::filesystem::path& dat_path() const noexcept {
        return m_dat_path;
    }
    [[nodiscard]] const std::filesystem::path& inf_path() const noexcept {
        return m_inf_path;
    }
    [[nodiscard]] const std::filesystem::path& base_path() const noexcept {
        return m_base;
    }

private:
    std::filesystem::path m_base;
    std::filesystem::path m_dat_path;
    std::filesystem::path m_inf_path;
};

constexpr std::string_view kSampleInf = R"(
 Data file name without the suffix          =  GBT_Lband_PSR
 Telescope used                             =  GBT
 Instrument used                            =  SPIGOT
 Object being observed                      =  B1937+21
 J2000 Right Ascension (hh:mm:ss.ssss)      =  19:39:38.5602
 J2000 Declination     (dd:mm:ss.ssss)      =  21:34:59.1430
 Data observed by                           =  R. Ransom
 Epoch of observation (MJD)                 =  53234.123456789012
 Barycentered?           (1=yes, 0=no)      =  1
 Number of bins in the data                 =  1024
 Width of each time bin (s)                 =  0.000081920000
 Any breaking sub-band cuts? (1=yes, 0=no)  =  0
 Dispersion measure (cm-3 pc)               =  71.02500
 Central freq of low channel (Mhz)          =  1170.000000
 Total bandwidth (Mhz)                      =  64.000000
 Number of channels                         =  1024
 Channel bandwidth (Mhz)                    =  0.062500
 Data analyzed by                           =  PRESTO
 Any additional notes:
    Test observation note.
)";

} // namespace

TEST_CASE("parse PRESTO .inf ASCII metadata", "[presto][inf]") {
    const auto hdr =
        psrio::formats::presto::parse_inf(kSampleInf, "GBT_Lband_PSR");

    REQUIRE(hdr.filename == "GBT_Lband_PSR");
    REQUIRE(hdr.telescope == "GBT");
    REQUIRE(hdr.backend == "SPIGOT");
    REQUIRE(hdr.source == "B1937+21");
    REQUIRE(hdr.ra == "19:39:38.5602");
    REQUIRE(hdr.dec == "21:34:59.1430");
    REQUIRE_THAT(hdr.raj, WithinAbs(193938.5602, 1e-4));
    REQUIRE_THAT(hdr.dej, WithinAbs(213459.1430, 1e-4));
    REQUIRE_THAT(hdr.tstart, WithinAbs(53234.123456789012, 1e-12));
    REQUIRE(hdr.barycentric == 1);
    REQUIRE(hdr.nsamples == 1024U);
    REQUIRE_THAT(hdr.tsamp, WithinAbs(0.00008192, 1e-10));
    REQUIRE_THAT(hdr.dm, WithinAbs(71.025, 1e-4));
    REQUIRE_THAT(hdr.fch1, WithinAbs(1170.0, 1e-4));
    REQUIRE_THAT(hdr.foff, WithinAbs(0.0625, 1e-4));
    REQUIRE(hdr.nbits == 32);
    REQUIRE(hdr.nchans == 1);
}

TEST_CASE("presto TimeSeriesReader satisfies TimeSeriesReader concept",
          "[presto][concept]") {
    static_assert(psrio::concepts::TimeSeriesReader<
                  psrio::formats::presto::TimeSeriesReader>);
}

TEST_CASE("presto TimeSeriesReader streams floats and exposes zero-copy views",
          "[presto][reader]") {
    std::vector<float> test_data(1024);
    for (std::size_t i = 0; i < test_data.size(); ++i) {
        test_data[i] = static_cast<float>(i) * 1.5F;
    }

    const TempPrestoFile presto(kSampleInf, test_data);

    // Can open by .dat, .inf, or base path
    psrio::formats::presto::TimeSeriesReader reader(presto.dat_path());
    REQUIRE(reader.header().nsamples == 1024U);
    REQUIRE(reader.header().source == "B1937+21");
    REQUIRE(reader.tell() == 0U);

    // Read first chunk
    std::array<float, 4> first_block{};
    REQUIRE(reader.read_samples(4, first_block) == 4U);
    REQUIRE(first_block[0] == 0.0F);
    REQUIRE(first_block[1] == 1.5F);
    REQUIRE(first_block[2] == 3.0F);
    REQUIRE(first_block[3] == 4.5F);
    REQUIRE(reader.tell() == 4U);

    // Zero-copy view
    const auto viewed = reader.view_samples(4);
    REQUIRE(viewed.size() == 4U);
    REQUIRE(viewed[0] == 6.0F);
    REQUIRE(viewed[1] == 7.5F);
    REQUIRE(viewed[2] == 9.0F);
    REQUIRE(viewed[3] == 10.5F);
    REQUIRE(reader.tell() == 8U);

    // Seek and tell
    reader.seek(1020);
    REQUIRE(reader.tell() == 1020U);
    std::array<float, 8> tail{};
    tail.fill(-1.0F);
    // Asking for 8 when only 4 remain returns 4
    REQUIRE(reader.read(4, std::span<float>{tail}.first(4)) == 4U);
    REQUIRE(tail[0] == 1020.0F * 1.5F);
    REQUIRE(tail[3] == 1023.0F * 1.5F);
    REQUIRE(reader.tell() == 1024U);

    // At EOF, read returns 0
    std::array<float, 2> empty_read{};
    REQUIRE(reader.read(0, std::span<float>{empty_read}.first(0)) == 0U);

    // Rewind
    reader.rewind();
    REQUIRE(reader.tell() == 0U);

    // Allocating read
    const auto read_all = reader.read_samples(1024);
    REQUIRE(read_all.size() == 1024U);
    REQUIRE(read_all[500] == 500.0F * 1.5F);
}

TEST_CASE(
    "presto TimeSeriesReader warns on truncated file vs declared nsamples",
    "[presto][reader]") {
    std::vector<float> short_data(10, 42.0F);
    const TempPrestoFile presto(kSampleInf,
                                short_data); // Declares 1024, only has 10

    const psrio::formats::presto::TimeSeriesReader reader(presto.inf_path());
    REQUIRE(reader.header().nsamples == 10U);
    REQUIRE_FALSE(reader.warnings().empty());
    REQUIRE_THAT(reader.warnings().front(), ContainsSubstring("declared 1024"));
}
