#include "psrio/block_source.hpp"
#include "psrio/clean.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/formats/memory.hpp"
#include "psrio/stats.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <optional>
#include <random>
#include <span>
#include <utility>
#include <vector>

namespace {

// 64 values with outliers; the expected masks come from sigpyproc 2.0.0
// (double_mad_mask and iqrm_mask) on the same values.
const std::vector<double> kValues{
    10.001, 10.299, 9.726,  9.109,  9.545,  25.0,   10.06,  11.34,
    9.508,  9.38,   10.49,  10.357, 10.105, 9.07,   9.971,  10.695,
    8.656,  -3.0,   8.099,  8.71,   8.158,  9.765,  8.733,  10.271,
    10.157, 9.813,  7.483,  9.461,  9.951,  10.113, 8.47,   9.522,
    9.021,  9.191,  11.061, 9.192,  9.967,  10.884, 9.416,  9.888,
    18.0,   10.064, 8.775,  10.076, 11.359, 8.453,  10.859, 10.119,
    9.359,  12.0,   14.5,   14.6,   14.7,   10.577, 9.811,  10.683,
    9.933,  10.667, 11.439, 9.324,  10.203, 9.537,  10.127, 8.813};

std::vector<std::size_t> flagged(const std::vector<std::uint8_t>& mask) {
    std::vector<std::size_t> out;
    for (std::size_t index = 0; index < mask.size(); ++index) {
        if (mask[index] != 0U) {
            out.push_back(index);
        }
    }
    return out;
}

psrio::BlockSource memory_source(const std::vector<float>& data,
                                 std::uint64_t nchans) {
    psrio::MemoryInfo info;
    info.nchans = nchans;
    info.nbits  = 32;
    info.tsamp  = 1e-3;
    info.fch1   = 1500.0;
    info.foff   = -1.0;
    std::vector<std::byte> bytes(data.size() * sizeof(float));
    std::memcpy(bytes.data(), data.data(), bytes.size());
    return psrio::BlockSource(psrio::MemoryBlock(info, std::move(bytes)));
}

} // namespace

static_assert(psrio::concepts::BlockReader<psrio::clean::Source>);

TEST_CASE("clean masks match sigpyproc's double MAD and IQRM", "[clean]") {
    using Index = std::vector<std::size_t>;
    CHECK(flagged(psrio::clean::mask_doublemad(kValues, 3.0)) ==
          Index{5, 17, 40, 50, 51, 52});
    CHECK(flagged(psrio::clean::mask_doublemad(kValues, 2.0)) ==
          Index{5, 17, 26, 40, 49, 50, 51, 52});
    CHECK(flagged(psrio::clean::mask_iqrm(kValues, 3.0)) ==
          Index{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 12, 13, 14,
                15, 16, 17, 18, 19, 20, 21, 22, 35, 36, 37, 38, 39, 40,
                41, 42, 43, 44, 45, 47, 48, 50, 51, 52, 53, 54, 55, 56});
    CHECK_THROWS_AS(psrio::clean::mask_doublemad(kValues, 0.0),
                    psrio::ValidationError);
    CHECK_THROWS_AS(psrio::clean::mask_iqrm(kValues, 3.0, 0),
                    psrio::ValidationError);
}

TEST_CASE("channel_mask flags a channel with an outlying variance", "[clean]") {
    constexpr std::size_t kRows = 2000;
    constexpr std::size_t kCols = 32;
    std::mt19937 rng(3);
    std::normal_distribution<float> gauss(0.0F, 1.0F);
    std::vector<float> data(kRows * kCols);
    for (std::size_t row = 0; row < kRows; ++row) {
        for (std::size_t col = 0; col < kCols; ++col) {
            const float width         = col == 9 ? 8.0F : 1.0F;
            data[(row * kCols) + col] = 100.0F + (width * gauss(rng));
        }
    }
    psrio::stats::ChannelStats stats(kCols);
    stats.push(data, kRows);
    const auto mask =
        psrio::clean::channel_mask(stats, psrio::clean::Method::kMad, 3.0);
    CHECK(mask[9] == 1U);
    const auto fill = psrio::clean::fill_value(stats, mask);
    CHECK(fill == Catch::Approx(100.0).margin(0.1));
}

TEST_CASE("ZeroDm removes the bandpass-weighted channel sum", "[clean]") {
    const std::vector<double> bandpass{1.0, 2.0, 3.0, 4.0};
    const psrio::clean::ZeroDm zerodm(bandpass);
    std::vector<float> block{2.0F, 3.0F, 4.0F, 5.0F, 1.0F, 2.0F, 3.0F, 4.0F};
    zerodm.apply(block, 2);
    // Row 0 sums to 14: x - 14 * bp / 10 + bp.
    CHECK(block[0] == Catch::Approx(2.0 - 1.4 + 1.0));
    CHECK(block[3] == Catch::Approx(5.0 - 5.6 + 4.0));
    // Row 1 equals the bandpass, so it is unchanged.
    CHECK(block[5] == Catch::Approx(2.0));
}

TEST_CASE("zap_rows replaces a sample range", "[clean]") {
    std::mt19937_64 rng(11);
    std::vector<float> block(30);
    for (std::size_t index = 0; index < block.size(); ++index) {
        block[index] = static_cast<float>(index);
    }
    const auto original = block;
    SECTION("zero") {
        psrio::clean::zap_rows(block, 10, 3, 2, 4, psrio::clean::ZapMode::kZero,
                               rng);
        CHECK(block[6] == 0.0F);
        CHECK(block[11] == 0.0F);
        CHECK(block[12] == original[12]);
    }
    SECTION("samples copies other rows") {
        psrio::clean::zap_rows(block, 10, 3, 2, 4,
                               psrio::clean::ZapMode::kSamples, rng);
        for (std::size_t row = 2; row < 4; ++row) {
            const auto first = static_cast<std::size_t>(block[row * 3]);
            CHECK(first % 3U == 0U);
            const auto source = first / 3U;
            CHECK((source < 2U || source >= 4U));
            CHECK(block[(row * 3) + 1] == original[(source * 3) + 1]);
        }
    }
    SECTION("gaussian draws near the other rows") {
        psrio::clean::zap_rows(block, 10, 3, 2, 4,
                               psrio::clean::ZapMode::kGaussian, rng);
        CHECK(block[6] != original[6]);
        CHECK(block[6] > -40.0F);
        CHECK(block[6] < 60.0F);
    }
}

TEST_CASE("clean::Source masks, removes zero-DM, and reverses channels",
          "[clean]") {
    constexpr std::size_t kRows = 50;
    constexpr std::size_t kCols = 4;
    std::vector<float> data(kRows * kCols);
    for (std::size_t index = 0; index < data.size(); ++index) {
        data[index] = static_cast<float>(index % 13);
    }
    auto input = memory_source(data, kCols);
    psrio::clean::Source cleaned(input, psrio::clean::Steps{
                                            .mask            = {0, 1, 0, 0},
                                            .fill            = 7.5F,
                                            .zerodm_bandpass = std::nullopt,
                                        });
    CHECK(cleaned.nbits() == 32);
    CHECK(cleaned.sigproc_header().nbits == 32);
    auto got = std::vector<float>(kRows * kCols);
    REQUIRE(cleaned.read_samples(kRows, got) == kRows);
    for (std::size_t row = 0; row < kRows; ++row) {
        CHECK(got[(row * kCols) + 1] == 7.5F);
        CHECK(got[(row * kCols) + 2] == data[(row * kCols) + 2]);
    }

    cleaned.rewind();
    cleaned.set_fswap(true);
    psrio::BlockSource erased(std::move(cleaned));
    const auto reversed = erased.read_samples<float>(1);
    CHECK(reversed[2] == 7.5F);
    CHECK(reversed[3] == data[0]);

    // The mask indexes columns in file order even when the input was set to
    // reverse its channels.
    auto swapped = memory_source(data, kCols);
    swapped.set_fswap(true);
    psrio::clean::Source file_order(swapped,
                                    psrio::clean::Steps{
                                        .mask            = {0, 1, 0, 0},
                                        .fill            = 7.5F,
                                        .zerodm_bandpass = std::nullopt,
                                    });
    CHECK_FALSE(swapped.fswap_enabled());
    auto first = std::vector<float>(kCols);
    REQUIRE(file_order.read_samples(1, first) == 1U);
    CHECK(first[0] == data[0]);
    CHECK(first[1] == 7.5F);

    CHECK_THROWS_AS(psrio::clean::Source(input,
                                         psrio::clean::Steps{
                                             .mask            = {1, 0},
                                             .fill            = 0.0F,
                                             .zerodm_bandpass = std::nullopt,
                                         }),
                    psrio::ValidationError);
}
