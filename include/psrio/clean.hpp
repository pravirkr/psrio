#pragma once

#include "psrio/astro.hpp"
#include "psrio/block_source.hpp"
#include "psrio/common/types.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/parallel.hpp"
#include "psrio/packed.hpp"
#include "psrio/stats.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <fmt/format.h>
#include <memory>
#include <numbers>
#include <optional>
#include <random>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

namespace psrio::clean {

/// Outlier test applied across channels to each channel statistic.
enum class Method : std::uint8_t {
    /// Double median absolute deviation (separate left and right scales).
    kMad,
    /// Inter-quartile range mitigation over lagged differences.
    kIqrm,
};

/// Normal-distribution constants used by sigpyproc's robust scales.
inline constexpr double kNormMad = 0.6744897501960817; // Phi^-1(0.75)
inline constexpr double kNormIqr = 1.3489795003921634; // IQR of N(0, 1)

namespace clean_detail {

/// numpy.median: the mean of the two middle values for an even count.
[[nodiscard]] inline double median(std::vector<double> values) {
    if (values.empty()) {
        throw ValidationError("psrio: median of no values");
    }
    const auto mid = values.size() / 2U;
    std::ranges::nth_element(values,
                             values.begin() + static_cast<std::ptrdiff_t>(mid));
    const double upper = values[mid];
    if (values.size() % 2U == 1U) {
        return upper;
    }
    const double lower = *std::max_element(
        values.begin(), values.begin() + static_cast<std::ptrdiff_t>(mid));
    return 0.5 * (lower + upper);
}

/// numpy.percentile with linear interpolation, @p q in [0, 100].
[[nodiscard]] inline double percentile(std::vector<double> values, double q) {
    if (values.empty()) {
        throw ValidationError("psrio: percentile of no values");
    }
    std::ranges::sort(values);
    const double index = (q / 100.0) * static_cast<double>(values.size() - 1U);
    const auto below   = static_cast<std::size_t>(std::floor(index));
    const auto above   = std::min(below + 1U, values.size() - 1U);
    const double frac  = index - static_cast<double>(below);
    return values[below] + ((values[above] - values[below]) * frac);
}

/// numpy.isclose(value, 0) with the default tolerances.
[[nodiscard]] inline bool near_zero(double value) noexcept {
    return std::abs(value) <= 1.0e-8;
}

/// |z| > threshold for z = (x - loc) / scale, computed in float32 as
/// sigpyproc's estimate_zscore. A zero scale counts as 1.
inline void flag_outliers(std::span<const float> data,
                          double loc,
                          std::span<const double> scales,
                          double threshold,
                          std::vector<std::uint8_t>& mask) {
    for (std::size_t index = 0; index < data.size(); ++index) {
        const double scale = near_zero(scales[index]) ? 1.0 : scales[index];
        const auto z       = static_cast<float>(
            static_cast<float>(static_cast<double>(data[index]) - loc) /
            static_cast<float>(scale));
        if (std::abs(z) > threshold) {
            mask[index] = 1;
        }
    }
}

[[nodiscard]] inline std::vector<float> as_float(std::span<const double> x) {
    std::vector<float> out(x.size());
    std::ranges::transform(
        x, out.begin(), [](double value) { return static_cast<float>(value); });
    return out;
}

inline void check_threshold(double threshold) {
    if (!std::isfinite(threshold) || threshold <= 0.0) {
        throw ValidationError(fmt::format(
            "psrio: clean threshold must be positive, got {}", threshold));
    }
}

} // namespace clean_detail

/// Outliers of @p values by the double MAD: |x - median| over the MAD of
/// the values on the same side of the median, as sigpyproc's
/// `double_mad_mask`. Returns 1 for an outlier.
/// @throws ValidationError if @p threshold is not positive.
[[nodiscard]] inline std::vector<std::uint8_t>
mask_doublemad(std::span<const double> values, double threshold) {
    clean_detail::check_threshold(threshold);
    std::vector<std::uint8_t> mask(values.size(), 0);
    if (values.empty()) {
        return mask;
    }
    const auto data = clean_detail::as_float(values);
    std::vector<double> wide(data.begin(), data.end());
    const double loc = clean_detail::median(wide);
    std::vector<double> left;
    std::vector<double> right;
    for (const double value : wide) {
        if (value <= loc) {
            left.push_back(std::abs(value - loc));
        }
        if (value >= loc) {
            right.push_back(std::abs(value - loc));
        }
    }
    const auto side_scale = [](const std::vector<double>& side) {
        const double mad = clean_detail::median(side) / kNormMad;
        if (!clean_detail::near_zero(mad)) {
            return mad;
        }
        double sum = 0.0;
        for (const double value : side) {
            sum += value;
        }
        return (sum / static_cast<double>(side.size())) /
               std::sqrt(2.0 / std::numbers::pi);
    };
    const double mad_left  = side_scale(left);
    const double mad_right = side_scale(right);
    std::vector<double> scales(wide.size());
    for (std::size_t index = 0; index < wide.size(); ++index) {
        scales[index] = wide[index] < loc ? mad_left : mad_right;
    }
    clean_detail::flag_outliers(data, loc, scales, threshold, mask);
    return mask;
}

/// Outliers of @p values by IQRM: for every lag within @p radius, the
/// differences x[i] - x[i + lag] (edge-padded) are scored against their
/// median and IQR, and any score above @p threshold flags x[i]. As
/// sigpyproc's `iqrm_mask`. Returns 1 for an outlier.
/// @throws ValidationError if @p threshold or @p radius is not positive.
[[nodiscard]] inline std::vector<std::uint8_t>
mask_iqrm(std::span<const double> values, double threshold, int radius = 5) {
    clean_detail::check_threshold(threshold);
    if (radius < 1) {
        throw ValidationError(
            fmt::format("psrio: IQRM radius must be >= 1, got {}", radius));
    }
    std::vector<std::uint8_t> mask(values.size(), 0);
    const auto count = static_cast<std::ptrdiff_t>(values.size());
    if (count == 0) {
        return mask;
    }
    std::vector<double> diffs(values.size());
    for (int lag = -radius; lag <= radius; ++lag) {
        if (lag == 0) {
            continue;
        }
        for (std::ptrdiff_t index = 0; index < count; ++index) {
            const auto other =
                std::clamp<std::ptrdiff_t>(index + lag, 0, count - 1);
            diffs[static_cast<std::size_t>(index)] =
                values[static_cast<std::size_t>(index)] -
                values[static_cast<std::size_t>(other)];
        }
        const auto data = clean_detail::as_float(diffs);
        const std::vector<double> wide(data.begin(), data.end());
        const double loc = clean_detail::median(wide);
        const double iqr = (clean_detail::percentile(wide, 75.0) -
                            clean_detail::percentile(wide, 25.0)) /
                           kNormIqr;
        const std::vector<double> scales(wide.size(), iqr);
        clean_detail::flag_outliers(data, loc, scales, threshold, mask);
    }
    return mask;
}

/// Channels whose variance, skewness, or kurtosis is an outlier across
/// channels, as sigpyproc's `RFIMask.apply_method`. Returns 1 for a channel
/// to mask.
[[nodiscard]] inline std::vector<std::uint8_t> channel_mask(
    const stats::ChannelStats& stats, Method method, double threshold) {
    const auto test = [&](const std::vector<double>& values) {
        return method == Method::kMad ? mask_doublemad(values, threshold)
                                      : mask_iqrm(values, threshold);
    };
    auto mask       = test(stats.var());
    const auto skew = test(stats.skew());
    const auto kurt = test(stats.kurtosis());
    for (std::size_t index = 0; index < mask.size(); ++index) {
        mask[index] =
            static_cast<std::uint8_t>(mask[index] | skew[index] | kurt[index]);
    }
    return mask;
}

/// Value for masked channels: the median of the unmasked channel means.
/// @throws ValidationError if every channel is masked.
[[nodiscard]] inline double fill_value(const stats::ChannelStats& stats,
                                       std::span<const std::uint8_t> mask) {
    std::vector<double> kept;
    const auto mean = stats.mean();
    for (std::size_t index = 0; index < mean.size(); ++index) {
        if (index >= mask.size() || mask[index] == 0U) {
            kept.push_back(mean[index]);
        }
    }
    if (kept.empty()) {
        throw ValidationError("psrio: every channel is masked");
    }
    return clean_detail::median(std::move(kept));
}

/// Set the masked columns of @p nrows rows to @p value.
inline void apply_mask(std::span<float> block,
                       std::uint64_t nrows,
                       std::span<const std::uint8_t> mask,
                       float value) {
    std::vector<std::size_t> masked;
    for (std::size_t index = 0; index < mask.size(); ++index) {
        if (mask[index] != 0U) {
            masked.push_back(index);
        }
    }
    if (masked.empty()) {
        return;
    }
    const auto ncols = mask.size();
    [[maybe_unused]] const bool wide =
        ::psrio::detail::parallel_worthwhile(nrows * masked.size());
    PSRIO_OMP(omp parallel for schedule(static) if (wide))
    for (std::uint64_t row = 0; row < nrows; ++row) {
        float* values = block.data() + (static_cast<std::size_t>(row) * ncols);
        for (const auto col : masked) {
            values[col] = value;
        }
    }
}

/// Zero-DM removal of Eatough, Keane & Lyne (2009), as PRESTO and
/// sigpyproc: each sample's channel sum, shared out by the bandpass
/// weights, is subtracted and the bandpass added back.
class ZeroDm {
public:
    /// @p bandpass is the mean of each column.
    explicit ZeroDm(std::span<const double> bandpass)
        : m_bandpass(bandpass.begin(), bandpass.end()),
          m_weights(bandpass.size()) {
        double total = 0.0;
        for (const double value : bandpass) {
            total += value;
        }
        for (std::size_t col = 0; col < m_weights.size(); ++col) {
            m_weights[col] = std::abs(total) < 1.0e-8
                                 ? 1.0 / static_cast<double>(m_weights.size())
                                 : bandpass[col] / total;
        }
    }

    [[nodiscard]] std::uint64_t ncols() const noexcept {
        return m_bandpass.size();
    }

    void apply(std::span<float> block, std::uint64_t nrows) const {
        const auto ncols = m_bandpass.size();
        if (block.size() != static_cast<std::size_t>(nrows) * ncols) {
            throw ValidationError("psrio: zero-DM block has the wrong size");
        }
        const double* weights  = m_weights.data();
        const double* bandpass = m_bandpass.data();
        [[maybe_unused]] const bool wide =
            ::psrio::detail::parallel_worthwhile(nrows * ncols);
        PSRIO_OMP(omp parallel for schedule(static) if (wide))
        for (std::uint64_t row = 0; row < nrows; ++row) {
            float* values =
                block.data() + (static_cast<std::size_t>(row) * ncols);
            double zerodm = 0.0;
            PSRIO_OMP(omp simd reduction(+ : zerodm))
            for (std::size_t col = 0; col < ncols; ++col) {
                zerodm += values[col];
            }
            PSRIO_OMP(omp simd)
            for (std::size_t col = 0; col < ncols; ++col) {
                values[col] =
                    static_cast<float>(static_cast<double>(values[col]) -
                                       (zerodm * weights[col]) + bandpass[col]);
            }
        }
    }

private:
    std::vector<double> m_bandpass;
    std::vector<double> m_weights;
};

/// How zap_rows() replaces samples.
enum class ZapMode : std::uint8_t {
    /// Zeros.
    kZero,
    /// Gaussian draws with each channel's mean and standard deviation over
    /// the block's other samples.
    kGaussian,
    /// Copies of random other samples of the block.
    kSamples,
};

/// Replace rows [@p lo, @p hi) of a [nrows, ncols] block, as sigpyproc's
/// `zap_block`. With no other rows in the block, the rows become zeros.
inline void zap_rows(std::span<float> block,
                     std::uint64_t nrows,
                     std::uint64_t ncols,
                     std::uint64_t lo,
                     std::uint64_t hi,
                     ZapMode mode,
                     std::mt19937_64& rng) {
    hi = std::min(hi, nrows);
    if (lo >= hi) {
        return;
    }
    const auto row_at = [&](std::uint64_t row) {
        return block.subspan(static_cast<std::size_t>(row * ncols),
                             static_cast<std::size_t>(ncols));
    };
    const auto kept = nrows - (hi - lo);
    if (mode == ZapMode::kZero || kept == 0U) {
        std::fill(block.begin() + static_cast<std::ptrdiff_t>(lo * ncols),
                  block.begin() + static_cast<std::ptrdiff_t>(hi * ncols),
                  0.0F);
        return;
    }
    const auto other = [&](std::uint64_t index) {
        return index < lo ? index : index + (hi - lo);
    };
    if (mode == ZapMode::kSamples) {
        std::uniform_int_distribution<std::uint64_t> pick(0, kept - 1U);
        for (auto row = lo; row < hi; ++row) {
            const auto source = row_at(other(pick(rng)));
            std::ranges::copy(source, row_at(row).begin());
        }
        return;
    }
    std::vector<double> mean(static_cast<std::size_t>(ncols), 0.0);
    std::vector<double> var(static_cast<std::size_t>(ncols), 0.0);
    for (std::uint64_t index = 0; index < kept; ++index) {
        const auto values = row_at(other(index));
        for (std::size_t col = 0; col < values.size(); ++col) {
            mean[col] += values[col];
        }
    }
    for (auto& value : mean) {
        value /= static_cast<double>(kept);
    }
    for (std::uint64_t index = 0; index < kept; ++index) {
        const auto values = row_at(other(index));
        for (std::size_t col = 0; col < values.size(); ++col) {
            const double d = values[col] - mean[col];
            var[col] += d * d;
        }
    }
    std::normal_distribution<double> gauss(0.0, 1.0);
    for (auto row = lo; row < hi; ++row) {
        auto values = row_at(row);
        for (std::size_t col = 0; col < values.size(); ++col) {
            double sd = std::sqrt(var[col] / static_cast<double>(kept));
            if (sd == 0.0) {
                sd = 1.0;
            }
            values[col] = static_cast<float>(mean[col] + (sd * gauss(rng)));
        }
    }
}

/// Cleaning applied by Source to each block, in this order.
struct Steps {
    /// 1 for each column (IF-major channel) to replace with fill.
    std::vector<std::uint8_t> mask;
    float fill{0.0F};
    /// Zero-DM removal with this bandpass (column means).
    std::optional<std::vector<double>> zerodm_bandpass;
};

/**
 * @brief A BlockSource read through channel masking and zero-DM removal.
 *
 * Models concepts::BlockReader as a float32 stream with the input's
 * channels, so it can feed digifil::Source or a writer without an
 * intermediate file. Masked columns are set to Steps::fill, then zero-DM
 * removal (when asked for) runs on the masked samples.
 */
class Source {
public:
    /// Read @p source, which this object then owns.
    /// @throws ValidationError if @p steps does not fit the input.
    Source(BlockSource&& source, Steps steps)
        : m_owned(std::make_unique<BlockSource>(std::move(source))),
          m_input(m_owned.get()),
          m_steps(std::move(steps)) {
        init();
    }

    /// Read @p source, which must outlive this object.
    Source(BlockSource& source, Steps steps)
        : m_input(&source),
          m_steps(std::move(steps)) {
        init();
    }

    [[nodiscard]] const Steps& steps() const noexcept { return m_steps; }

    [[nodiscard]] std::uint64_t nchans() const { return m_input->nchans(); }
    [[nodiscard]] std::uint64_t nifs() const { return m_input->nifs(); }
    [[nodiscard]] static int nbits() noexcept { return 32; }
    [[nodiscard]] static SampleType sample_type() noexcept {
        return SampleType::kFloat32;
    }
    [[nodiscard]] std::uint64_t bytes_per_sample() const {
        return m_row * sizeof(float);
    }
    [[nodiscard]] std::uint64_t nsamples() const { return m_input->nsamples(); }
    [[nodiscard]] bool has_nsamples() const { return m_input->has_nsamples(); }
    [[nodiscard]] double tsamp() const { return m_input->tsamp(); }
    [[nodiscard]] double tstart() const { return m_input->tstart(); }
    [[nodiscard]] std::string_view source_name() const {
        return m_input->source_name();
    }
    [[nodiscard]] std::string_view telescope() const {
        return m_input->telescope();
    }
    [[nodiscard]] double raj() const { return m_input->raj(); }
    [[nodiscard]] double dej() const { return m_input->dej(); }
    [[nodiscard]] double fch1() const { return m_input->fch1(); }
    [[nodiscard]] double foff() const { return m_input->foff(); }
    [[nodiscard]] double bandwidth() const { return m_input->bandwidth(); }
    [[nodiscard]] double center_frequency() const {
        return m_input->center_frequency();
    }
    [[nodiscard]] int beam() const { return m_input->beam(); }
    [[nodiscard]] double spectra_rate() const {
        return m_input->spectra_rate();
    }
    [[nodiscard]] std::time_t utc_start() const { return m_input->utc_start(); }
    [[nodiscard]] std::optional<double> refdm() const {
        return m_input->refdm();
    }
    /// The input's SIGPROC header with float32 samples.
    [[nodiscard]] formats::sigproc::FilterbankHeader sigproc_header() const {
        auto header        = sigproc_header_of(*m_input);
        header.nbits       = 32;
        header.signed_data = std::nullopt;
        return header;
    }
    [[nodiscard]] std::uint64_t tell() const { return m_input->tell(); }

    /// Reverse the channels of each IF on read.
    void set_fswap(bool enable) noexcept { m_fswap = enable; }
    [[nodiscard]] bool fswap_enabled() const noexcept { return m_fswap; }

    void seek(std::uint64_t sample) { m_input->seek(sample); }
    void rewind() { m_input->rewind(); }
    void skip(std::int64_t delta) { m_input->skip(delta); }

    std::uint64_t read_samples(std::uint64_t count, std::span<float> dest) {
        if (dest.size() != static_cast<std::size_t>(count * m_row)) {
            throw ValidationError(fmt::format(
                "psrio: destination has {} values but the request needs {}",
                dest.size(), count * m_row));
        }
        const auto got = m_input->read_samples(count, dest);
        const auto out = dest.first(static_cast<std::size_t>(got * m_row));
        if (!m_steps.mask.empty()) {
            apply_mask(out, got, m_steps.mask, m_steps.fill);
        }
        if (m_zerodm.has_value()) {
            m_zerodm->apply(out, got);
        }
        if (m_fswap && got > 0U) {
            reverse_unpacked_channels(out, got, nchans(), nifs());
        }
        return got;
    }

    std::uint64_t read_samples(std::uint64_t /*count*/,
                               std::span<std::uint8_t> /*dest*/) {
        throw ValidationError("psrio: cleaned samples unpack to float");
    }
    std::uint64_t read_samples(std::uint64_t /*count*/,
                               std::span<std::uint16_t> /*dest*/) {
        throw ValidationError("psrio: cleaned samples unpack to float");
    }

    /// Little-endian float32 samples.
    std::uint64_t read_block(std::uint64_t count, std::span<std::byte> dest) {
        if (dest.size() !=
            static_cast<std::size_t>(count * bytes_per_sample())) {
            throw ValidationError(
                "psrio: destination does not hold the requested samples");
        }
        m_floats.resize(static_cast<std::size_t>(count * m_row));
        const auto got = read_samples(count, m_floats);
        std::memcpy(dest.data(), m_floats.data(),
                    static_cast<std::size_t>(got * bytes_per_sample()));
        return got;
    }

    std::uint64_t read_bytes(std::uint64_t nbytes, std::span<std::byte> dest) {
        const auto stride = bytes_per_sample();
        if (nbytes == 0U || nbytes % stride != 0U || dest.size() != nbytes) {
            throw ValidationError("psrio: byte request must be a positive "
                                  "multiple of the sample stride");
        }
        return read_block(nbytes / stride, dest) * stride;
    }

private:
    void init() {
        m_row = m_input->nifs() * m_input->nchans();
        if (!m_steps.mask.empty() && m_steps.mask.size() != m_row) {
            throw ValidationError(fmt::format(
                "psrio: channel mask has {} entries for {} channels",
                m_steps.mask.size(), m_row));
        }
        if (m_steps.zerodm_bandpass.has_value()) {
            if (m_steps.zerodm_bandpass->size() != m_row) {
                throw ValidationError(fmt::format(
                    "psrio: zero-DM bandpass has {} entries for {} channels",
                    m_steps.zerodm_bandpass->size(), m_row));
            }
            m_zerodm.emplace(*m_steps.zerodm_bandpass);
        }
    }

    std::unique_ptr<BlockSource> m_owned;
    BlockSource* m_input{nullptr};
    Steps m_steps;
    std::optional<ZeroDm> m_zerodm;
    std::uint64_t m_row{0};
    bool m_fswap{false};
    std::vector<float> m_floats;
};

} // namespace psrio::clean
