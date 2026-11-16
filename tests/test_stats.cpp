#include "psrio/block_source.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/formats/memory.hpp"
#include "psrio/stats.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <random>
#include <span>
#include <vector>

namespace {

struct Reference {
    std::vector<double> mean, var, skew, kurt;
    std::vector<float> lo, hi;
};

/// Two-pass moments of a [nrows, ncols] array, as sigpyproc defines them.
Reference
two_pass(const std::vector<float>& data, std::size_t nrows, std::size_t ncols) {
    Reference ref;
    for (std::size_t col = 0; col < ncols; ++col) {
        double sum = 0.0;
        float lo   = data[col];
        float hi   = data[col];
        for (std::size_t row = 0; row < nrows; ++row) {
            const float value = data[(row * ncols) + col];
            sum += value;
            lo = std::min(lo, value);
            hi = std::max(hi, value);
        }
        const double mean = sum / static_cast<double>(nrows);
        double m2         = 0.0;
        double m3         = 0.0;
        double m4         = 0.0;
        for (std::size_t row = 0; row < nrows; ++row) {
            const double d = data[(row * ncols) + col] - mean;
            m2 += d * d;
            m3 += d * d * d;
            m4 += d * d * d * d;
        }
        const auto n = static_cast<double>(nrows);
        ref.mean.push_back(mean);
        ref.var.push_back(m2 / n);
        ref.skew.push_back(m2 > 0.0 ? std::sqrt(n) * m3 / std::pow(m2, 1.5)
                                    : 0.0);
        ref.kurt.push_back((m2 > 0.0 ? n * m4 / (m2 * m2) : 0.0) - 3.0);
        ref.lo.push_back(lo);
        ref.hi.push_back(hi);
    }
    return ref;
}

/// Gaussian columns with column-dependent offsets and widths, one column
/// constant and one heavy-tailed, so every moment is exercised.
std::vector<float> make_data(std::size_t nrows, std::size_t ncols) {
    std::mt19937 rng(42);
    std::normal_distribution<float> gauss(0.0F, 1.0F);
    std::exponential_distribution<float> expo(1.0F);
    std::vector<float> data(nrows * ncols);
    for (std::size_t row = 0; row < nrows; ++row) {
        for (std::size_t col = 0; col < ncols; ++col) {
            float value = 1000.0F + (10.0F * static_cast<float>(col)) +
                          ((1.0F + static_cast<float>(col % 7)) * gauss(rng));
            if (col == 3) {
                value = 77.0F;
            }
            if (col == 5) {
                value = 500.0F + (20.0F * expo(rng));
            }
            data[(row * ncols) + col] = value;
        }
    }
    return data;
}

void check_against(const psrio::stats::ChannelStats& stats,
                   const Reference& ref) {
    const auto var  = stats.var();
    const auto skew = stats.skew();
    const auto kurt = stats.kurtosis();
    for (std::size_t col = 0; col < ref.mean.size(); ++col) {
        INFO("column " << col);
        CHECK(stats.mean()[col] == Catch::Approx(ref.mean[col]).epsilon(1e-12));
        CHECK(var[col] ==
              Catch::Approx(ref.var[col]).epsilon(1e-9).margin(1e-12));
        CHECK(skew[col] ==
              Catch::Approx(ref.skew[col]).epsilon(1e-7).margin(1e-9));
        CHECK(kurt[col] ==
              Catch::Approx(ref.kurt[col]).epsilon(1e-7).margin(1e-9));
        CHECK(stats.minima()[col] == ref.lo[col]);
        CHECK(stats.maxima()[col] == ref.hi[col]);
    }
}

} // namespace

TEST_CASE("ChannelStats matches two-pass moments across blocks", "[stats]") {
    constexpr std::size_t kRows = 5000;
    constexpr std::size_t kCols = 130; // more than two tiles, one partial
    const auto data             = make_data(kRows, kCols);
    const auto ref              = two_pass(data, kRows, kCols);

    SECTION("one block") {
        psrio::stats::ChannelStats stats(kCols);
        stats.push(data, kRows);
        CHECK(stats.count() == kRows);
        check_against(stats, ref);
    }
    SECTION("uneven blocks") {
        psrio::stats::ChannelStats stats(kCols);
        std::size_t row = 0;
        for (const std::size_t size : {1U, 2U, 997U, 3000U, 1000U}) {
            stats.push(
                std::span<const float>(data).subspan(row * kCols, size * kCols),
                size);
            row += size;
        }
        REQUIRE(row == kRows);
        check_against(stats, ref);
    }
    SECTION("merged halves") {
        psrio::stats::ChannelStats first(kCols);
        psrio::stats::ChannelStats second(kCols);
        first.push(std::span<const float>(data).first(1234 * kCols), 1234);
        second.push(std::span<const float>(data).subspan(1234 * kCols),
                    kRows - 1234);
        first.merge(second);
        check_against(first, ref);
    }
    SECTION("compute() over a BlockSource") {
        psrio::MemoryInfo info;
        info.nchans = kCols;
        info.nbits  = 32;
        info.tsamp  = 1e-3;
        info.fch1   = 1500.0;
        info.foff   = -1.0;
        std::vector<std::byte> bytes(data.size() * sizeof(float));
        std::memcpy(bytes.data(), data.data(), bytes.size());
        psrio::BlockSource source(psrio::MemoryBlock(info, std::move(bytes)));
        const auto stats = psrio::stats::compute(
            source, psrio::stats::ComputeOptions{.gulp = 777});
        check_against(stats, ref);

        const auto window = psrio::stats::compute(
            source, psrio::stats::ComputeOptions{
                        .start  = 100,
                        .nsamps = 200,
                        .gulp   = 64,
                        .order  = psrio::stats::Moments::kMean});
        CHECK(window.count() == 200U);
        double sum = 0.0;
        for (std::size_t row = 100; row < 300; ++row) {
            sum += data[row * kCols];
        }
        CHECK(window.sum()[0] == Catch::Approx(sum).epsilon(1e-12));
        CHECK_THROWS_AS(window.var(), psrio::ValidationError);
    }
}
