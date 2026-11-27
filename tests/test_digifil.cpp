#include "psrio/block_source.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/digifil.hpp"
#include "psrio/formats/sigproc.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {

namespace digifil = psrio::digifil;
using psrio::formats::sigproc::FilterbankHeader;
using psrio::formats::sigproc::FilterbankReader;

const std::filesystem::path kDataDir{PSRIO_TEST_DATA_DIR};

class TempPath {
public:
    TempPath() {
        static std::atomic<int> counter{0};
        m_path = std::filesystem::temp_directory_path() /
                 ("psrio_digifil_" + std::to_string(::getpid()) + "_" +
                  std::to_string(counter.fetch_add(1)) + ".fil");
    }
    ~TempPath() { std::filesystem::remove(m_path); }

    TempPath(const TempPath&)            = delete;
    TempPath& operator=(const TempPath&) = delete;
    TempPath(TempPath&&)                 = delete;
    TempPath& operator=(TempPath&&)      = delete;

    [[nodiscard]] const std::filesystem::path& path() const noexcept {
        return m_path;
    }

private:
    std::filesystem::path m_path;
};

std::uint64_t fnv1a(std::span<const std::byte> bytes) {
    std::uint64_t hash = 0xcbf29ce484222325ULL;
    for (const auto byte : bytes) {
        hash ^= std::to_integer<std::uint64_t>(byte);
        hash *= 0x100000001b3ULL;
    }
    return hash;
}

std::vector<std::byte> payload(const std::filesystem::path& path) {
    const FilterbankReader reader(path);
    std::ifstream in(path, std::ios::binary);
    in.seekg(static_cast<std::streamoff>(reader.header().header_bytes));
    std::vector<char> raw((std::istreambuf_iterator<char>(in)),
                          std::istreambuf_iterator<char>());
    std::vector<std::byte> bytes(raw.size());
    for (std::size_t i = 0; i < raw.size(); ++i) {
        bytes[i] = static_cast<std::byte>(raw[i]);
    }
    return bytes;
}

digifil::InputLayout
layout(std::uint64_t nchans, std::uint64_t nifs = 1, double foff = -1.0) {
    return {
        .nchans = nchans,
        .nifs   = nifs,
        .nbits  = 8,
        .tsamp  = 1e-3,
        .foff   = foff,
    };
}

} // namespace

TEST_CASE("digifil quantise matches the dspsr digitiser",
          "[digifil][quantise]") {
    CHECK(digifil::quantise(0.0F, 127.5F, 1.0F, 255U) == 128U);
    CHECK(digifil::quantise(-0.6F, 127.5F, 1.0F, 255U) == 127U);
    CHECK(digifil::quantise(1000.0F, 127.5F, 1.0F, 255U) == 255U);
    CHECK(digifil::quantise(-1000.0F, 127.5F, 1.0F, 255U) == 0U);
    CHECK(digifil::quantise(std::numeric_limits<float>::quiet_NaN(), 127.5F,
                            1.0F, 255U) == 0U);
    CHECK(digifil::quantise(std::numeric_limits<float>::infinity(), 1.5F, 1.0F,
                            3U) == 3U);
    CHECK(digifil::quantise(0.4F, 0.0F, 1.0F, 3U) == 0U);
    CHECK(digifil::quantise(0.5F, 0.0F, 1.0F, 3U) == 1U);
}

TEST_CASE("digifil block_samples follows the dspsr RAM limit",
          "[digifil][gulp]") {
    // 64 channels at 2 bits with one float copy: 4.25 * 64 = 272 bytes.
    const auto expected =
        static_cast<std::uint64_t>(digifil::kDefaultMaxRamBytes / 272.0);
    const auto rounded = ((expected + 3U) / 4U) * 4U;
    CHECK(digifil::block_samples(2, 1, 64, 1, 4,
                                 digifil::kDefaultMaxRamBytes) == rounded);
    CHECK(digifil::block_samples(8, 1, 1U << 30U, 3, 4, 1.0) == 4U);
    CHECK_THROWS_AS(digifil::block_samples(8, 1, 64, 1, 4, 0.0),
                    psrio::ValidationError);
}

TEST_CASE("digifil Rescaler normalises each column per interval",
          "[digifil][rescale]") {
    digifil::Rescaler rescaler(2, 4, false);
    // Column 0 alternates 1/3 (mean 2, sigma 1); column 1 is constant 5.
    std::vector<float> data{1, 5, 3, 5, 1, 5, 3, 5};
    rescaler.process(data, 4);
    const std::vector<float> expected{-1, 0, 1, 0, -1, 0, 1, 0};
    CHECK(data == expected);
    CHECK(rescaler.scales()[1] == 1.0F);
}

TEST_CASE("digifil Rescaler updates within the first interval per block",
          "[digifil][rescale]") {
    digifil::Rescaler rescaler(1, 100, false);
    std::vector<float> first{2, 4};
    rescaler.process(first, 2);
    CHECK(first == std::vector<float>{-1, 1});
    std::vector<float> second{0, 6};
    rescaler.process(second, 2);
    // Statistics now cover {2, 4, 0, 6}: mean 3, variance 5.
    CHECK(rescaler.offsets()[0] == -3.0F);
    CHECK(second[1] == Catch::Approx(3.0 / std::sqrt(5.0)));
}

TEST_CASE("digifil Rescaler in constant mode freezes after one interval",
          "[digifil][rescale]") {
    digifil::Rescaler rescaler(1, 2, true);
    std::vector<float> data{1, 3, 100, 200, 7};
    rescaler.process(data, 5);
    CHECK(data == std::vector<float>{-1, 1, 98, 198, 5});
    CHECK_THROWS_AS(digifil::Rescaler(0, 2, false), psrio::ValidationError);
    CHECK_THROWS_AS(digifil::Rescaler(1, 0, false), psrio::ValidationError);
}

TEST_CASE("digifil Config validation", "[digifil][config]") {
    digifil::Config config;
    CHECK_NOTHROW(config.validate());
    for (const int nbits : {0, 3, 32, -8}) {
        config.nbits = nbits;
        CHECK_THROWS_AS(config.validate(), psrio::ValidationError);
    }
    config           = {};
    config.scale_fac = 0.0F;
    CHECK_THROWS_AS(config.validate(), psrio::ValidationError);
    config                 = {};
    config.rescale_seconds = -1.0;
    CHECK_THROWS_AS(config.validate(), psrio::ValidationError);
    config      = {};
    config.gulp = 0U;
    CHECK_THROWS_AS(config.validate(), psrio::ValidationError);
}

TEST_CASE("digifil Pipeline rejects layouts it cannot produce",
          "[digifil][pipeline]") {
    digifil::Config config;
    config.fscrunch = 3;
    CHECK_THROWS_AS(digifil::Pipeline(config, layout(8)),
                    psrio::ValidationError);
    config.fscrunch = 16;
    CHECK_THROWS_AS(digifil::Pipeline(config, layout(8)),
                    psrio::ValidationError);
    config             = {};
    config.poln_select = 2;
    CHECK_THROWS_AS(digifil::Pipeline(config, layout(8, 2)),
                    psrio::ValidationError);
    config       = {};
    config.nbits = 1;
    CHECK_THROWS_AS(digifil::Pipeline(config, layout(4)),
                    psrio::ValidationError);
    config                 = {};
    config.rescale_seconds = 1e-4;
    CHECK_THROWS_AS(digifil::Pipeline(config, layout(8)),
                    psrio::ValidationError);
}

TEST_CASE("digifil Pipeline without rescaling sums and flips channels",
          "[digifil][pipeline]") {
    digifil::Config config;
    config.rescale_seconds = 0.0;
    config.fscrunch        = 2;
    config.tscrunch        = 2;
    // Ascending input band: output channels are reversed.
    digifil::Pipeline pipeline(config, layout(4, 1, 1.0));
    REQUIRE(pipeline.nchans_out() == 2U);
    REQUIRE(pipeline.bytes_per_sample() == 2U);

    std::vector<float> block{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    // Rows 0-1 complete one output sample; row 2 is carried.
    const auto out = pipeline.push(block, 3);
    REQUIRE(out.size() == 2U);
    // Scale 1/4 undoes the sum of four values: (3+4+7+8)/4 = 5.5 -> 6.
    CHECK(std::to_integer<int>(out[0]) == 6);
    CHECK(std::to_integer<int>(out[1]) == 4);

    std::vector<float> tail{13, 14, 15, 16};
    const auto out2 = pipeline.push(tail, 1);
    REQUIRE(out2.size() == 2U);
    CHECK(std::to_integer<int>(out2[0]) == 14);
    CHECK(std::to_integer<int>(out2[1]) == 12);
}

TEST_CASE("digifil Pipeline selects one IF", "[digifil][pipeline]") {
    digifil::Config config;
    config.rescale_seconds = 0.0;
    config.poln_select     = 1;
    digifil::Pipeline pipeline(config, layout(2, 2));
    REQUIRE(pipeline.nifs_out() == 1U);
    std::vector<float> block{1, 2, 30, 40, 5, 6, 70, 80};
    const auto out = pipeline.push(block, 2);
    REQUIRE(out.size() == 4U);
    CHECK(std::to_integer<int>(out[0]) == 30);
    CHECK(std::to_integer<int>(out[1]) == 40);
    CHECK(std::to_integer<int>(out[2]) == 70);
    CHECK(std::to_integer<int>(out[3]) == 80);
}

TEST_CASE("digifil output_header keeps the band centre", "[digifil][header]") {
    FilterbankHeader in;
    in.nchans      = 64;
    in.nbits       = 2;
    in.tsamp       = 0.00032;
    in.tstart      = 50000.0;
    in.fch1        = 1510.0;
    in.foff        = -1.09;
    in.signed_data = 1;
    digifil::Config config;
    config.nbits    = 2;
    config.tscrunch = 3;
    config.fscrunch = 4;
    const digifil::Pipeline pipeline(config, {
                                                 .nchans = 64,
                                                 .nifs   = 1,
                                                 .nbits  = 2,
                                                 .tsamp  = 0.00032,
                                                 .foff   = -1.09,
                                             });
    const auto out = digifil::output_header(in, pipeline, 1000);
    CHECK(out.nchans == 16);
    CHECK(out.nbits == 2);
    CHECK(out.tsamp.value() == Catch::Approx(0.00096));
    CHECK(out.foff.value() == Catch::Approx(-4.36));
    CHECK(out.fch1.value() == Catch::Approx(1508.365));
    CHECK(out.tstart.value() ==
          Catch::Approx(50000.0 + (1000 * 0.00032 / 86400.0)));
    CHECK(out.signed_data == std::optional<std::int8_t>{0});
    CHECK_FALSE(out.declared_nsamples.has_value());
}

TEST_CASE("digifil run matches sigpyproc spp_digifil on tutorial.fil",
          "[digifil][run]") {
    struct Case {
        digifil::Config config;
        std::uint64_t nbytes;
        std::uint64_t hash;
    };
    // Payload sizes and FNV-1a hashes of spp_digifil (sigpyproc) output.
    const auto make = [](int nbits, std::uint64_t t, std::uint64_t f,
                         double seconds, bool constant, float scale,
                         std::optional<std::uint64_t> gulp) {
        digifil::Config config;
        config.nbits            = nbits;
        config.tscrunch         = t;
        config.fscrunch         = f;
        config.rescale_seconds  = seconds;
        config.rescale_constant = constant;
        config.scale_fac        = scale;
        config.gulp             = gulp;
        return config;
    };
    const auto test_case =
        GENERATE_COPY(Case{make(8, 1, 1, 10.0, false, 1.0F, std::nullopt),
                           12001280U, 0x727351dd14512dd4ULL},
                      Case{make(2, 3, 4, 0.5, false, 1.0F, std::nullopt),
                           250024U, 0x7cb769a67dda47b9ULL},
                      Case{make(4, 2, 1, 0.3, true, 1.0F, 1000U), 3000320U,
                           0x58e043579d639303ULL},
                      Case{make(16, 5, 1, 0.0, false, 1.0F, std::nullopt),
                           4800512U, 0x30de26581fadcfc6ULL},
                      Case{make(-32, 7, 2, 10.0, false, 1.0F, 333U), 3428864U,
                           0x0b7481d324e3fdf0ULL},
                      Case{make(1, 1, 1, 1.0, false, 2.0F, 777U), 1500160U,
                           0x09e4a26ff64f1c0aULL});

    const auto input = kDataDir / "tutorial.fil";
    FilterbankReader reader(input);
    const auto header = reader.header();
    psrio::BlockSource source(std::move(reader));
    const TempPath out;
    const auto summary =
        digifil::run(source, header, out.path(), test_case.config);
    CHECK(summary.nsamples_in == 187520U);

    const auto bytes = payload(out.path());
    CHECK(bytes.size() == test_case.nbytes);
    CHECK(fnv1a(bytes) == test_case.hash);
}

TEST_CASE("digifil run honours the sample range", "[digifil][run]") {
    const auto input = kDataDir / "tutorial.fil";
    FilterbankReader reader(input);
    const auto header = reader.header();
    psrio::BlockSource source(std::move(reader));
    digifil::Config config;
    config.rescale_seconds = 0.0;
    config.tscrunch        = 4;
    const TempPath out;
    const auto summary =
        digifil::run(source, header, out.path(), config, 100, 1001);
    CHECK(summary.nsamples_in == 1001U);
    CHECK(summary.nsamples_out == 250U);

    const FilterbankReader result(out.path());
    CHECK(result.nsamples() == 250U);
    CHECK(result.header().tstart.value() ==
          Catch::Approx(header.mjd_after_samples(100)));
    CHECK_THROWS_AS(digifil::run(source, header, out.path(), config, 187520),
                    psrio::ValidationError);
    CHECK_THROWS_AS(digifil::run(source, header, out.path(), config, 0, 187521),
                    psrio::ValidationError);
}

TEST_CASE("digifil Source streams the payload run() writes",
          "[digifil][source]") {
    static_assert(psrio::concepts::BlockReader<digifil::Source>);
    const auto make = [](int nbits, std::uint64_t t, std::uint64_t f,
                         double seconds, std::optional<std::uint64_t> gulp) {
        digifil::Config config;
        config.nbits           = nbits;
        config.tscrunch        = t;
        config.fscrunch        = f;
        config.rescale_seconds = seconds;
        config.gulp            = gulp;
        return config;
    };
    const auto config = GENERATE_COPY(make(8, 1, 1, 10.0, std::nullopt),
                                      make(2, 3, 4, 0.5, std::nullopt),
                                      make(-32, 7, 2, 10.0, 333U));
    const auto input  = kDataDir / "tutorial.fil";
    FilterbankReader reader(input);
    const auto header = reader.header();
    psrio::BlockSource source(std::move(reader));
    const TempPath out;
    (void)digifil::run(source, header, out.path(), config);
    const auto expected = payload(out.path());

    digifil::Source stream(psrio::BlockSource(FilterbankReader(input)), config);
    const auto stride = stream.bytes_per_sample();
    CHECK(stream.nsamples() * stride == expected.size());
    CHECK(stream.header().nchans ==
          FilterbankReader(out.path()).header().nchans);
    CHECK(stream.refdm() == header.refdm);

    // Uneven reads, through the type-erased handle.
    psrio::BlockSource wrapped(std::move(stream));
    std::vector<std::byte> got;
    std::vector<std::byte> piece(997U * stride);
    while (true) {
        const auto n = wrapped.read_block(997U, piece);
        got.insert(got.end(), piece.begin(),
                   piece.begin() + static_cast<std::ptrdiff_t>(n * stride));
        if (n < 997U) {
            break;
        }
    }
    CHECK(got == expected);
    CHECK(wrapped.tell() == wrapped.nsamples());

    // A rewind restarts the pipeline from the first sample.
    wrapped.rewind();
    const auto again = wrapped.read_block(wrapped.nsamples());
    CHECK(again == expected);
}

TEST_CASE("digifil Source unpacks 8-bit output and keeps refdm",
          "[digifil][source]") {
    auto header  = FilterbankReader(kDataDir / "tutorial.fil").header();
    header.refdm = 413.0;
    psrio::BlockSource source(FilterbankReader(kDataDir / "tutorial.fil"));
    const digifil::Config config;
    digifil::Source stream(source, header, config, 1000, 20000);
    CHECK(stream.refdm() == std::optional<double>(413.0));
    CHECK(stream.header().refdm == std::optional<double>(413.0));
    CHECK(stream.nsamples() == 20000U);
    const auto bytes = stream.read_block(100);
    stream.seek(0);
    std::vector<std::uint8_t> levels(100U * stream.nchans());
    CHECK(stream.read_samples(100, std::span<std::uint8_t>(levels)) == 100U);
    for (std::size_t i = 0; i < levels.size(); ++i) {
        REQUIRE(levels[i] == std::to_integer<std::uint8_t>(bytes[i]));
    }
    stream.set_fswap(true);
    stream.seek(0);
    const auto swapped = stream.read_block(1);
    CHECK(std::to_integer<int>(swapped.front()) ==
          std::to_integer<int>(bytes[stream.nchans() - 1U]));
}

TEST_CASE("digifil Rescaler clipping ignores impulsive outliers",
          "[digifil][rescale]") {
    // Gaussian noise of sigma 2 around 50, with every 20th value a spike of
    // +400 standing in for impulsive RFI.
    std::uint64_t state = 12345;
    const auto uniform  = [&state] {
        state = (state * 6364136223846793005ULL) + 1442695040888963407ULL;
        return (static_cast<double>(state >> 11U) + 0.5) / 9007199254740992.0;
    };
    std::vector<float> values(6000);
    for (std::size_t i = 0; i < values.size(); ++i) {
        const double gauss = std::sqrt(-2.0 * std::log(uniform())) *
                             std::cos(6.283185307179586 * uniform());
        values[i]          = static_cast<float>(50.0 + (2.0 * gauss) +
                                                (i % 20U == 0U ? 400.0 : 0.0));
    }
    digifil::Rescaler plain(1, 2000, false);
    digifil::Rescaler robust(1, 2000, false, 3.0);
    for (std::size_t start = 0; start < values.size(); start += 1000U) {
        std::vector<float> block(
            values.begin() + static_cast<std::ptrdiff_t>(start),
            values.begin() + static_cast<std::ptrdiff_t>(start + 1000U));
        auto copy = block;
        plain.process(block, 1000);
        robust.process(copy, 1000);
    }
    // Plain statistics see the spikes: sigma near 87. Clipped ones see the
    // noise: sigma near 2.
    CHECK(1.0F / plain.scales()[0] > 50.0F);
    CHECK(1.0F / robust.scales()[0] == Catch::Approx(2.0).epsilon(0.1));
    CHECK(-robust.offsets()[0] == Catch::Approx(50.0).epsilon(0.01));

    // A clip of zero is dspsr's rescaler.
    digifil::Rescaler zero(1, 2000, false, 0.0);
    digifil::Rescaler fresh(1, 2000, false);
    auto a = values;
    auto b = values;
    zero.process(a, 6000);
    fresh.process(b, 6000);
    CHECK(a == b);

    digifil::Config config;
    config.rescale_clip = -1.0;
    CHECK_THROWS_AS(config.validate(), psrio::ValidationError);
    config.rescale_clip = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(config.validate(), psrio::ValidationError);

    // NaN values do not poison clipped rescaler sums.
    digifil::Rescaler nan_safe(1, 4, false, 3.0);
    std::vector<float> with_nan{10.0F, 12.0F,
                                std::numeric_limits<float>::quiet_NaN(), 14.0F};
    nan_safe.process(with_nan, 4);
    CHECK(std::isfinite(nan_safe.offsets()[0]));
    CHECK(std::isfinite(nan_safe.scales()[0]));
    CHECK(-nan_safe.offsets()[0] == Catch::Approx(12.0F));
}

TEST_CASE("digifil Pipeline pscrunch and run progress callback",
          "[digifil][pipeline][progress]") {
    digifil::Config config;
    config.rescale_seconds = 0.0;
    config.pscrunch        = true;
    digifil::Pipeline pipeline(config, layout(2, 2));
    REQUIRE(pipeline.nifs_out() == 1U);
    std::vector<float> block{1, 2, 30, 40, 5, 6, 70, 80};
    const auto out = pipeline.push(block, 2);
    REQUIRE(out.size() == 4U);
    CHECK(std::to_integer<int>(out[0]) == 31);
    CHECK(std::to_integer<int>(out[1]) == 42);
    CHECK(std::to_integer<int>(out[2]) == 75);
    CHECK(std::to_integer<int>(out[3]) == 86);

    FilterbankReader reader(kDataDir / "tutorial.fil");
    const auto header = reader.header();
    psrio::BlockSource source(std::move(reader));
    digifil::Config run_cfg;
    run_cfg.gulp          = 20000U;
    std::uint64_t prev_in = 0;
    std::size_t callbacks = 0;
    const TempPath out_path;
    const auto summary =
        digifil::run(source, header, out_path.path(), run_cfg, 0, std::nullopt,
                     [&](const digifil::Progress& prog) {
                         CHECK(prog.samples_in_done >= prev_in);
                         CHECK(prog.samples_in_done <= prog.samples_in_total);
                         prev_in = prog.samples_in_done;
                         ++callbacks;
                     });
    CHECK(callbacks > 1U);
    CHECK(prev_in == summary.nsamples_in);
}
