#include "psrio/detail/exceptions.hpp"
#include "psrio/scrunch.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <span>
#include <vector>

namespace {

std::vector<float> iota_rows(std::size_t nrows, std::size_t row_len) {
    std::vector<float> data(nrows * row_len);
    std::iota(data.begin(), data.end(), 0.0F);
    return data;
}

} // namespace

TEST_CASE("fscrunch_inplace sums channel groups per row",
          "[scrunch][fscrunch]") {
    auto data      = iota_rows(3, 6);
    const auto out = psrio::fscrunch_inplace(std::span<float>(data), 3, 6, 3);
    REQUIRE(out.size() == 6U);
    CHECK(out.data() == data.data());
    const std::vector<float> expected{3, 12, 21, 30, 39, 48};
    CHECK(std::vector<float>(out.begin(), out.end()) == expected);
}

TEST_CASE("fscrunch_inplace with factor 1 is a no-op", "[scrunch][fscrunch]") {
    auto data           = iota_rows(2, 4);
    const auto original = data;
    const auto out = psrio::fscrunch_inplace(std::span<float>(data), 2, 4, 1);
    CHECK(out.size() == 8U);
    CHECK(data == original);
}

TEST_CASE("fscrunch_inplace rejects bad factors and short buffers",
          "[scrunch][fscrunch]") {
    auto data = iota_rows(2, 6);
    CHECK_THROWS_AS(psrio::fscrunch_inplace(data, 2, 6, 0),
                    psrio::ValidationError);
    CHECK_THROWS_AS(psrio::fscrunch_inplace(data, 2, 6, 4),
                    psrio::ValidationError);
    CHECK_THROWS_AS(psrio::fscrunch_inplace(data, 3, 6, 2),
                    psrio::ValidationError);
}

TEST_CASE("TimeScruncher output does not depend on block boundaries",
          "[scrunch][tscrunch]") {
    constexpr std::size_t kRowLen = 3;
    constexpr std::size_t kRows   = 23;
    const auto factor             = GENERATE(1U, 2U, 3U, 5U, 7U);
    const auto gulp               = GENERATE(1U, 2U, 4U, 6U, 23U);
    auto data                     = iota_rows(kRows, kRowLen);

    std::vector<float> expected;
    for (std::size_t g = 0; g + factor <= kRows; g += factor) {
        for (std::size_t c = 0; c < kRowLen; ++c) {
            float acc = data[(g * kRowLen) + c];
            for (std::size_t r = 1; r < factor; ++r) {
                acc += data[((g + r) * kRowLen) + c];
            }
            expected.push_back(acc);
        }
    }

    psrio::TimeScruncher scruncher(kRowLen, factor);
    std::vector<float> got;
    auto input = data;
    for (std::size_t start = 0; start < kRows; start += gulp) {
        const auto nrows = std::min<std::size_t>(gulp, kRows - start);
        const auto out   = scruncher.push(
            std::span<float>(input).subspan(start * kRowLen), nrows);
        got.insert(got.end(), out.begin(), out.end());
    }
    CHECK(got == expected);
    CHECK(input == data);
    CHECK(scruncher.pending() == kRows % factor);
    scruncher.reset();
    CHECK(scruncher.pending() == 0U);
}

TEST_CASE("TimeScruncher rejects invalid arguments", "[scrunch][tscrunch]") {
    CHECK_THROWS_AS(psrio::TimeScruncher(0, 2), psrio::ValidationError);
    CHECK_THROWS_AS(psrio::TimeScruncher(2, 0), psrio::ValidationError);
    psrio::TimeScruncher scruncher(4, 2);
    std::vector<float> data(7);
    CHECK_THROWS_AS(scruncher.push(data, 2), psrio::ValidationError);
}
