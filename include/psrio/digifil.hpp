#pragma once

#include "psrio/astro.hpp"
#include "psrio/block_source.hpp"
#include "psrio/common/types.hpp"
#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/packed_bits.hpp"
#include "psrio/detail/parallel.hpp"
#include "psrio/detail/skip.hpp"
#include "psrio/formats/psrfits/reader.hpp"
#include "psrio/formats/sigproc/header.hpp"
#include "psrio/formats/sigproc/writer.hpp"
#include "psrio/packed.hpp"
#include "psrio/scrunch.hpp"

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fmt/format.h>
#include <future>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

/// Filterbank redigitisation following dspsr's `digifil`.
///
/// The pipeline selects an IF, sums channel groups and time samples, rescales
/// each channel to zero mean and unit variance over an interval, and
/// digitises to 1, 2, 4, 8, or 16 bits, or writes float32. Output channels
/// always descend in frequency (negative foff), as dspsr writes them.
namespace psrio::digifil {

/// dspsr SigProcDigitizer mean for @p nbits output levels.
/// @throws ValidationError for an unsupported bit depth.
[[nodiscard]] inline float digitiser_mean(int nbits) {
    switch (nbits) {
    case 1:
        return 0.5F;
    case 2:
        return 1.5F;
    case 4:
        return 7.5F;
    case 8:
        return 127.5F;
    case 16:
        return 32768.0F;
    default:
        throw ValidationError(
            fmt::format("psrio: digifil cannot digitise to {} bits", nbits));
    }
}

/// Sigma range spanned by the 4, 8, and 16-bit digitiser.
inline constexpr float kDigitiserSigma = 6.0F;

/// dspsr's default RAM limit for sizing input blocks, in bytes.
inline constexpr double kDefaultMaxRamBytes = 256.0 * 1024.0 * 1024.0;

/// Input block resolution of SIGPROC data in dspsr, in samples.
inline constexpr std::uint64_t kSigprocResolution = 4;

/// Pipeline settings, named after the dspsr `digifil` options.
struct Config {
    /// Output bits per sample: 1, 2, 4, 8, 16, or -32 for float32 (-b).
    int nbits{8};
    /// Rescale interval in seconds; 0 disables rescaling (-I).
    double rescale_seconds{10.0};
    /// Freeze the offset and scale after the first interval (-c).
    bool rescale_constant{false};
    /// Robust rescaling: once a channel has statistics, values more than
    /// this many standard deviations from its mean are left out of the next
    /// statistics, so bursts of RFI do not inflate the scale. 0 keeps
    /// dspsr's plain mean and variance (--rescale-clip).
    double rescale_clip{0.0};
    /// Time decimation factor; 0 means 1 (-t).
    std::uint64_t tscrunch{1};
    /// Frequency decimation factor; 0 means 1 (-f).
    std::uint64_t fscrunch{1};
    /// Extra scale applied before digitising (-s).
    float scale_fac{1.0F};
    /// Keep only this input IF, or all IFs when negative (-P).
    int poln_select{-1};
    /// Input samples per block. Rescaling statistics depend on the block
    /// boundaries during the first interval. Defaults to dspsr's block for
    /// max_ram_mb (default_gulp()).
    std::optional<std::uint64_t> gulp;
    /// RAM limit in MiB that sets the default block, as dspsr `-U`.
    double max_ram_mb{256.0};
    /// Search-mode PSRFITS: apply `DAT_SCL`, `DAT_OFFS`, and `ZERO_OFF`
    /// (dspsr `--scloffs`). Used by decode_options() when the caller opens
    /// the file; the pipeline itself does not read it.
    bool scloffs{false};

    [[nodiscard]] double max_ram_bytes() const noexcept {
        return max_ram_mb * 1024.0 * 1024.0;
    }

    /// @throws ValidationError if a setting is out of range.
    void validate() const {
        if (nbits != -32) {
            (void)digitiser_mean(nbits);
        }
        if (!std::isfinite(rescale_seconds) || rescale_seconds < 0.0) {
            throw ValidationError(fmt::format(
                "psrio: digifil rescale interval must be finite and >= 0, "
                "got {}",
                rescale_seconds));
        }
        if (!std::isfinite(rescale_clip) || rescale_clip < 0.0) {
            throw ValidationError(fmt::format(
                "psrio: digifil rescale clip must be finite and >= 0, got {}",
                rescale_clip));
        }
        if (!std::isfinite(scale_fac) || scale_fac == 0.0F) {
            throw ValidationError(fmt::format(
                "psrio: digifil scale factor must be finite and non-zero, "
                "got {}",
                scale_fac));
        }
        if (gulp.has_value() && *gulp == 0U) {
            throw ValidationError("psrio: digifil gulp must be positive");
        }
        if (!std::isfinite(max_ram_mb) || max_ram_mb <= 0.0) {
            throw ValidationError(fmt::format(
                "psrio: digifil RAM limit must be positive, got {} MB",
                max_ram_mb));
        }
    }
};

/// How dspsr's digifil decodes search-mode PSRFITS for @p config: DAT_SCL
/// and DAT_OFFS only with Config::scloffs, no DAT_WTS, and `-P` picking a
/// stored product (AA, BB, ...) instead of an IF of the decoded stream.
/// Open the file with these options and run the pipeline with
/// poln_select -1.
[[nodiscard]] inline formats::psrfits::DecodeOptions
decode_options(const Config& config) noexcept {
    return formats::psrfits::DecodeOptions{
        .scloffs = config.scloffs,
        .weights = false,
        .product = config.poln_select >= 0 ? config.poln_select : -1,
    };
}

/// Input samples per block that dspsr reads under a RAM limit.
///
/// dspsr holds floor(max_ram / nbyte) samples, where nbyte is the input
/// sample size plus four bytes per float buffer, rounded up to a multiple of
/// @p resolution.
/// @throws ValidationError if @p max_ram_bytes is not positive.
[[nodiscard]] inline std::uint64_t block_samples(int nbits_in,
                                                 std::uint64_t npol,
                                                 std::uint64_t nchans,
                                                 int copies,
                                                 std::uint64_t resolution,
                                                 double max_ram_bytes) {
    if (!std::isfinite(max_ram_bytes) || max_ram_bytes <= 0.0) {
        throw ValidationError("psrio: digifil RAM limit must be positive");
    }
    const double nbyte = ((static_cast<double>(std::abs(nbits_in)) / 8.0) +
                          (4.0 * static_cast<double>(copies))) *
                         static_cast<double>(npol) *
                         static_cast<double>(nchans);
    const auto nblock =
        static_cast<std::uint64_t>(std::floor(max_ram_bytes / nbyte));
    resolution         = std::max<std::uint64_t>(resolution, 1U);
    const auto rounded = ((nblock + resolution - 1U) / resolution) * resolution;
    return std::max(resolution, rounded);
}

/// Block length `spp_digifil` uses for one search-mode PSRFITS file.
///
/// The stored bit depth, polarisation count, and `NSBLK` size the block.
/// Summing the two feeds counts as one extra float buffer, which is how
/// sigpyproc sizes a coherence or `PPQQ` file. The calibrated stream is one
/// float IF, so sizing from that stream picks a different block.
[[nodiscard]] inline std::uint64_t
psrfits_block_samples(int stored_nbits,
                      std::uint64_t npol,
                      std::uint64_t nchans,
                      std::uint64_t subint_samples,
                      bool sum_feeds,
                      std::uint64_t tscrunch,
                      std::uint64_t fscrunch,
                      double max_ram_bytes = kDefaultMaxRamBytes) {
    const auto time_factor = std::max<std::uint64_t>(tscrunch, 1U);
    const auto chan_factor = std::max<std::uint64_t>(fscrunch, 1U);
    const int copies       = 1 + (chan_factor > 1U ? 1 : 0) +
                             (time_factor > 1U ? 1 : 0) + (sum_feeds ? 1 : 0);
    return block_samples(stored_nbits, npol, nchans, copies, subint_samples,
                         max_ram_bytes);
}

/// Per-column rescaling to zero mean and unit variance, as dspsr `Rescale`.
///
/// Sums and sums of squares accumulate in double over intervals of
/// nsample rows, which may span blocks. During the first interval the
/// statistics are recomputed at the end of every block, so block boundaries
/// matter there. Each row is rescaled as (x + offset) * scale with the latest
/// statistics.
///
/// With a clip of k > 0, a column that already has statistics adds only the
/// values within k standard deviations of its mean to the next statistics,
/// which then average over the values kept. The first block is taken whole.
class Rescaler {
public:
    /// @throws ValidationError if @p ncols or @p nsample is 0, or @p clip is
    ///         negative.
    Rescaler(std::uint64_t ncols,
             std::uint64_t nsample,
             bool constant,
             double clip = 0.0)
        : m_nsample(nsample),
          m_constant(constant),
          m_clip(clip),
          m_sums(static_cast<std::size_t>(ncols), 0.0),
          m_sumsq(static_cast<std::size_t>(ncols), 0.0),
          m_offset(static_cast<std::size_t>(ncols), 0.0F),
          m_scale(static_cast<std::size_t>(ncols), 1.0F) {
        if (ncols == 0U || nsample == 0U) {
            throw ValidationError(
                "psrio: rescaler needs at least one column and sample");
        }
        if (!(clip >= 0.0) || !std::isfinite(clip)) {
            throw ValidationError("psrio: rescaler clip must be >= 0");
        }
        if (clip > 0.0) {
            m_count.assign(static_cast<std::size_t>(ncols), 0U);
            m_lower.assign(static_cast<std::size_t>(ncols), 0.0);
            m_upper.assign(static_cast<std::size_t>(ncols), 0.0);
        }
    }

    [[nodiscard]] std::span<const float> offsets() const noexcept {
        return m_offset;
    }

    [[nodiscard]] std::span<const float> scales() const noexcept {
        return m_scale;
    }

    /// Rescale @p nrows rows of ncols values in place.
    void process(std::span<float> data, std::uint64_t nrows) {
        const auto ncols = m_sums.size();
        const auto rows  = static_cast<std::size_t>(nrows);
        if (data.size() / ncols < rows) {
            throw ValidationError("psrio: rescaler input is smaller than rows");
        }
        if (m_constant && !m_first) {
            apply(data.first(rows * ncols));
            return;
        }
        std::size_t start = 0;
        while (start < rows) {
            const auto step = std::min<std::size_t>(
                rows - start, static_cast<std::size_t>(m_nsample - m_isample));
            auto end = start + step;
            if (m_first || !m_constant) {
                accumulate(data.subspan(start * ncols, step * ncols));
                m_isample += step;
            }
            if (m_isample == m_nsample || (m_first && end == rows)) {
                compute();
                if (m_isample == m_nsample) {
                    m_first   = false;
                    m_isample = 0;
                    std::ranges::fill(m_sums, 0.0);
                    std::ranges::fill(m_sumsq, 0.0);
                    std::ranges::fill(m_count, 0U);
                }
            }
            if (m_constant && !m_first) {
                end = rows;
            }
            apply(data.subspan(start * ncols, (end - start) * ncols));
            start = end;
        }
    }

private:
    void accumulate(std::span<const float> rows) {
        // Threads own tiles of columns, and each column still sums its rows
        // in order, so the sums match a serial pass bit for bit.
        constexpr std::size_t kTile = 64;
        const auto ncols            = m_sums.size();
        const auto nrows            = rows.size() / ncols;
        const auto ntiles           = (ncols + kTile - 1U) / kTile;
        [[maybe_unused]] const bool wide =
            ::psrio::detail::parallel_worthwhile(rows.size());
        const float* data = rows.data();
        double* sums      = m_sums.data();
        double* sumsq     = m_sumsq.data();
        if (m_clip > 0.0) {
            const bool clip     = m_have_stats;
            std::uint64_t* kept = m_count.data();
            const double* lower = m_lower.data();
            const double* upper = m_upper.data();
            PSRIO_OMP(omp parallel for schedule(static) if (wide))
            for (std::size_t tile = 0; tile < ntiles; ++tile) {
                const auto first = tile * kTile;
                const auto last  = std::min(first + kTile, ncols);
                for (std::size_t row = 0; row < nrows; ++row) {
                    const float* values = data + (row * ncols);
                    for (std::size_t col = first; col < last; ++col) {
                        const auto value = static_cast<double>(values[col]);
                        if (clip &&
                            (value < lower[col] || value > upper[col])) {
                            continue;
                        }
                        sums[col] += value;
                        sumsq[col] += value * value;
                        kept[col] += 1U;
                    }
                }
            }
            return;
        }
        PSRIO_OMP(omp parallel for schedule(static) if (wide))
        for (std::size_t tile = 0; tile < ntiles; ++tile) {
            const auto first = tile * kTile;
            const auto last  = std::min(first + kTile, ncols);
            for (std::size_t row = 0; row < nrows; ++row) {
                const float* values = data + (row * ncols);
                for (std::size_t col = first; col < last; ++col) {
                    const auto value = static_cast<double>(values[col]);
                    sums[col] += value;
                    sumsq[col] += value * value;
                }
            }
        }
    }

    void compute() {
        if (m_clip > 0.0) {
            compute_clipped();
            return;
        }
        const double recip = 1.0 / static_cast<double>(m_isample);
        for (std::size_t col = 0; col < m_sums.size(); ++col) {
            const double mean     = m_sums[col] * recip;
            const double variance = (m_sumsq[col] * recip) - (mean * mean);
            m_offset[col]         = static_cast<float>(-mean);
            m_scale[col] = variance > 0.0
                               ? static_cast<float>(1.0 / std::sqrt(variance))
                               : 1.0F;
        }
    }

    void compute_clipped() {
        for (std::size_t col = 0; col < m_sums.size(); ++col) {
            if (m_count[col] == 0U) {
                continue; // every value was clipped: keep the last statistics
            }
            const double recip    = 1.0 / static_cast<double>(m_count[col]);
            const double mean     = m_sums[col] * recip;
            const double variance = (m_sumsq[col] * recip) - (mean * mean);
            const double sigma    = variance > 0.0 ? std::sqrt(variance) : 0.0;
            m_offset[col]         = static_cast<float>(-mean);
            m_scale[col] = sigma > 0.0 ? static_cast<float>(1.0 / sigma) : 1.0F;
            // A constant channel has no spread to clip against.
            const double half = sigma > 0.0
                                    ? m_clip * sigma
                                    : std::numeric_limits<double>::infinity();
            m_lower[col]      = mean - half;
            m_upper[col]      = mean + half;
        }
        m_have_stats = true;
    }

    void apply(std::span<float> rows) const {
        const auto ncols = m_offset.size();
        const auto nrows = rows.size() / ncols;
        [[maybe_unused]] const bool wide =
            ::psrio::detail::parallel_worthwhile(rows.size());
        float* data         = rows.data();
        const float* offset = m_offset.data();
        const float* scale  = m_scale.data();
        PSRIO_OMP(omp parallel for schedule(static) if (wide))
        for (std::size_t row = 0; row < nrows; ++row) {
            float* values = data + (row * ncols);
            for (std::size_t col = 0; col < ncols; ++col) {
                const float shifted = values[col] + offset[col];
                values[col]         = shifted * scale[col];
            }
        }
    }

    std::uint64_t m_nsample;
    bool m_constant;
    double m_clip{0.0};
    bool m_have_stats{false};
    bool m_first{true};
    std::uint64_t m_isample{0};
    std::vector<double> m_sums;
    std::vector<double> m_sumsq;
    std::vector<float> m_offset;
    std::vector<float> m_scale;
    std::vector<std::uint64_t> m_count;
    std::vector<double> m_lower;
    std::vector<double> m_upper;
};

/// Quantise one value as dspsr's SigProcDigitizer: round(value * scale +
/// mean) clamped to [0, @p digi_max]. NaN maps to 0.
[[nodiscard]] inline std::uint32_t quantise(float value,
                                            float mean,
                                            float scale,
                                            std::uint32_t digi_max) noexcept {
    constexpr double kLimit = 2147483647.0;
    const float scaled      = value * scale;
    const float shifted     = scaled + mean;
    double rounded          = static_cast<double>(shifted) + 0.5;
    if (!(rounded >= -kLimit)) {
        rounded = -kLimit;
    } else if (rounded > kLimit) {
        rounded = kLimit;
    }
    const auto level = static_cast<std::int64_t>(rounded);
    if (level < 0) {
        return 0;
    }
    return std::min(static_cast<std::uint32_t>(level), digi_max);
}

/// Geometry of the data entering the pipeline.
struct InputLayout {
    std::uint64_t nchans{0};
    std::uint64_t nifs{1};
    int nbits{8};
    double tsamp{0.0};
    double foff{0.0};
};

/// Stateful digifil pipeline acting on successive input blocks.
///
/// push() takes a block of [nsamp, nifs, nchans] floats and returns the
/// encoded SIGPROC payload for the output samples it completes.
class Pipeline {
public:
    /// @throws ValidationError if @p config does not fit @p input.
    Pipeline(const Config& config, const InputLayout& input)
        : m_config(config),
          m_in(input) {
        config.validate();
        if (input.nchans == 0U || input.nifs == 0U) {
            throw ValidationError("psrio: digifil input has no channels");
        }
        if (!(input.tsamp > 0.0)) {
            throw ValidationError(
                "psrio: digifil input tsamp must be positive");
        }
        const auto ffac = std::max<std::uint64_t>(config.fscrunch, 1U);
        if (ffac > input.nchans || input.nchans % ffac != 0U) {
            throw ValidationError(
                fmt::format("psrio: fscrunch factor {} must divide nchans {}",
                            ffac, input.nchans));
        }
        if (config.poln_select >= 0 &&
            std::cmp_greater_equal(config.poln_select, input.nifs)) {
            throw ValidationError(
                fmt::format("psrio: IF {} is out of range for {}-IF input",
                            config.poln_select, input.nifs));
        }
        m_ffac     = ffac;
        m_tfac     = std::max<std::uint64_t>(config.tscrunch, 1U);
        m_nifs_out = config.poln_select >= 0 ? 1U : input.nifs;
        m_nchans   = input.nchans / ffac;
        m_flip     = input.foff > 0.0;
        m_is_float = config.nbits == -32;
        m_nbits    = m_is_float ? 32 : config.nbits;

        float scale_total = 1.0F;
        if (m_ffac > 1U) {
            scale_total *= static_cast<float>(m_ffac);
        }
        if (m_tfac > 1U) {
            scale_total *= static_cast<float>(m_tfac);
        }
        const bool rescale_on = config.rescale_seconds != 0.0;
        m_scale_in            = rescale_on ? 1.0F : scale_total;
        if (rescale_on) {
            const double rate_out =
                1.0 / (input.tsamp * static_cast<double>(m_tfac));
            const auto nsample =
                static_cast<std::uint64_t>(config.rescale_seconds * rate_out);
            if (nsample == 0U) {
                throw ValidationError(fmt::format(
                    "psrio: rescale interval {} s is shorter than one output "
                    "sample",
                    config.rescale_seconds));
            }
            m_rescaler.emplace(m_nifs_out * m_nchans, nsample,
                               config.rescale_constant, config.rescale_clip);
        }
        if (!m_is_float) {
            const float mean = digitiser_mean(m_nbits);
            const float scale =
                (rescale_on && m_nbits >= 4) ? mean / kDigitiserSigma : 1.0F;
            m_digi_scale = scale / (m_scale_in * config.scale_fac);
            m_digi_max =
                (std::uint32_t{1} << static_cast<unsigned>(m_nbits)) - 1U;
            m_pol_mean.assign(static_cast<std::size_t>(m_nifs_out), 0.0F);
            for (std::size_t pol = 0; pol < m_pol_mean.size(); ++pol) {
                if (rescale_on || pol >= 2U) {
                    m_pol_mean[pol] = mean;
                }
            }
            if (m_nbits < 8 &&
                (m_nchans % (8U / static_cast<unsigned>(m_nbits))) != 0U) {
                throw ValidationError(fmt::format(
                    "psrio: {} output channels do not pack into whole bytes "
                    "at {} bits",
                    m_nchans, m_nbits));
            }
        }
        m_tscrunch = TimeScruncher(m_nifs_out * m_nchans, m_tfac);
    }

    [[nodiscard]] std::uint64_t nchans_out() const noexcept { return m_nchans; }
    [[nodiscard]] std::uint64_t nifs_out() const noexcept { return m_nifs_out; }
    [[nodiscard]] int nbits_out() const noexcept { return m_nbits; }
    [[nodiscard]] std::uint64_t tscrunch() const noexcept { return m_tfac; }
    [[nodiscard]] std::uint64_t fscrunch() const noexcept { return m_ffac; }

    /// Bytes in one output sample.
    [[nodiscard]] std::uint64_t bytes_per_sample() const noexcept {
        return (m_nifs_out * m_nchans * static_cast<std::uint64_t>(m_nbits)) /
               8U;
    }

    /// Float buffers dspsr allocates for this configuration, for
    /// block_samples().
    [[nodiscard]] int float_copies() const noexcept {
        return 1 + (m_ffac > 1U ? 1 : 0) + (m_tfac > 1U ? 1 : 0);
    }

    /// Process @p nsamp input samples from @p block, which is modified.
    /// @return Encoded output samples; valid until the next push().
    /// @throws ValidationError if @p block is smaller than nsamp samples.
    std::span<const std::byte> push(std::span<float> block,
                                    std::uint64_t nsamp) {
        const auto row_in = m_in.nifs * m_in.nchans;
        if (block.size() / row_in < nsamp) {
            throw ValidationError(
                "psrio: digifil block is smaller than its sample count");
        }
        auto data = block.first(static_cast<std::size_t>(nsamp * row_in));
        if (m_config.poln_select >= 0) {
            const auto pol   = static_cast<std::size_t>(m_config.poln_select);
            const auto width = static_cast<std::size_t>(m_in.nchans);
            for (std::size_t s = 0; s < nsamp; ++s) {
                std::memmove(data.data() + (s * width),
                             data.data() + ((s * row_in) + (pol * width)),
                             width * sizeof(float));
            }
            data = data.first(static_cast<std::size_t>(nsamp) * width);
        }
        data = fscrunch_inplace(data, nsamp * m_nifs_out, m_in.nchans, m_ffac);

        const auto row_out = static_cast<std::size_t>(m_nifs_out * m_nchans);
        const auto work    = m_tscrunch.push(data, nsamp);
        const auto nout    = work.size() / row_out;
        if (nout == 0U) {
            return {};
        }
        if (m_rescaler.has_value()) {
            m_rescaler->process(work, nout);
        }
        return encode(work, nout);
    }

private:
    std::span<const std::byte> encode(std::span<const float> data,
                                      std::size_t nout) {
        const auto nrows  = nout * static_cast<std::size_t>(m_nifs_out);
        const auto nchans = static_cast<std::size_t>(m_nchans);
        const auto nbytes = static_cast<std::size_t>(nout * bytes_per_sample());
        m_bytes.resize(nbytes);
        const auto channel = [this, nchans](std::size_t ochan) {
            return m_flip ? nchans - 1U - ochan : ochan;
        };

        [[maybe_unused]] const bool wide =
            ::psrio::detail::parallel_worthwhile(nrows * nchans);
        if (m_is_float) {
            PSRIO_OMP(omp parallel for schedule(static) if (wide))
            for (std::size_t row = 0; row < nrows; ++row) {
                for (std::size_t ochan = 0; ochan < nchans; ++ochan) {
                    const float value = ::psrio::detail::to_little_endian(
                        data[(row * nchans) + channel(ochan)] / m_scale_in);
                    std::memcpy(m_bytes.data() +
                                    (((row * nchans) + ochan) * sizeof(float)),
                                &value, sizeof(float));
                }
            }
            return m_bytes;
        }

        const auto npol = static_cast<std::size_t>(m_nifs_out);
        if (m_nbits == 16) {
            PSRIO_OMP(omp parallel for schedule(static) if (wide))
            for (std::size_t row = 0; row < nrows; ++row) {
                const float mean = m_pol_mean[row % npol];
                for (std::size_t ochan = 0; ochan < nchans; ++ochan) {
                    const auto level = ::psrio::detail::to_little_endian(
                        static_cast<std::uint16_t>(
                            quantise(data[(row * nchans) + channel(ochan)],
                                     mean, m_digi_scale, m_digi_max)));
                    std::memcpy(m_bytes.data() +
                                    (((row * nchans) + ochan) * 2U),
                                &level, 2U);
                }
            }
            return m_bytes;
        }

        const bool packed = m_nbits < 8;
        if (packed) {
            m_levels.resize(nrows * nchans);
        }
        PSRIO_OMP(omp parallel for schedule(static) if (wide))
        for (std::size_t row = 0; row < nrows; ++row) {
            const float mean = m_pol_mean[row % npol];
            for (std::size_t ochan = 0; ochan < nchans; ++ochan) {
                const auto level = static_cast<std::uint8_t>(
                    quantise(data[(row * nchans) + channel(ochan)], mean,
                             m_digi_scale, m_digi_max));
                if (packed) {
                    m_levels[(row * nchans) + ochan] = level;
                } else {
                    m_bytes[(row * nchans) + ochan] = std::byte{level};
                }
            }
        }
        if (packed) {
            ::psrio::detail::pack_sub_byte(m_levels, m_bytes, m_nbits,
                                           BitOrder::kLsbFirst);
        }
        return m_bytes;
    }

    Config m_config;
    InputLayout m_in;
    std::uint64_t m_ffac{1};
    std::uint64_t m_tfac{1};
    std::uint64_t m_nifs_out{1};
    std::uint64_t m_nchans{0};
    bool m_flip{false};
    bool m_is_float{false};
    int m_nbits{8};
    float m_scale_in{1.0F};
    float m_digi_scale{1.0F};
    std::uint32_t m_digi_max{0};
    std::vector<float> m_pol_mean;
    std::optional<Rescaler> m_rescaler;
    TimeScruncher m_tscrunch{1, 1};
    std::vector<std::uint8_t> m_levels;
    std::vector<std::byte> m_bytes;
};

/// SIGPROC header of the digifil output for @p input.
///
/// The band keeps its centre and width; channels descend in frequency. tstart
/// moves to input sample @p start. A signed flag is cleared, since the output
/// levels are unsigned. A frequency table is dropped.
/// @throws ValidationError if @p input has no fch1 and foff.
[[nodiscard]] inline formats::sigproc::FilterbankHeader
output_header(const formats::sigproc::FilterbankHeader& input,
              const Pipeline& pipeline,
              std::uint64_t start = 0) {
    if (!input.fch1.has_value() || !input.foff.has_value() ||
        !input.tsamp.has_value()) {
        throw ValidationError(
            "psrio: digifil needs fch1, foff, and tsamp in the input header");
    }
    auto out            = input;
    const auto nin      = static_cast<double>(input.nchans);
    const auto nout     = static_cast<double>(pipeline.nchans_out());
    const double bw     = std::abs(*input.foff) * nin;
    const double centre = *input.fch1 + (0.5 * (*input.foff) * (nin - 1.0));
    out.frequency_table.clear();
    out.nchans = static_cast<std::int32_t>(pipeline.nchans_out());
    out.nifs   = static_cast<std::int32_t>(pipeline.nifs_out());
    out.nbits  = pipeline.nbits_out();
    out.foff   = -bw / nout;
    out.fch1   = centre + (0.5 * bw) - (0.5 * bw / nout);
    out.tsamp  = *input.tsamp * static_cast<double>(pipeline.tscrunch());
    if (input.tstart.has_value()) {
        out.tstart = input.mjd_after_samples(start);
    }
    out.declared_nsamples.reset();
    if (input.samples_are_signed()) {
        out.signed_data = 0;
    }
    return out;
}

/// Default input block for @p source: dspsr's block for a RAM limit
/// (256 MB unless given).
[[nodiscard]] inline std::uint64_t
default_gulp(const BlockSource& source,
             const Pipeline& pipeline,
             double max_ram_bytes = kDefaultMaxRamBytes) {
    return block_samples(source.nbits(), source.nifs(), source.nchans(),
                         pipeline.float_copies(), kSigprocResolution,
                         max_ram_bytes);
}

/// SIGPROC header that describes @p source, from the reader's own header when
/// it has one, else from its metadata getters.
[[nodiscard]] inline formats::sigproc::FilterbankHeader
header_of(const BlockSource& source) {
    if (auto header = source.sigproc_header()) {
        return *std::move(header);
    }
    formats::sigproc::FilterbankHeader header;
    header.nchans    = static_cast<std::int32_t>(source.nchans());
    header.nifs      = static_cast<std::int32_t>(source.nifs());
    header.nbits     = source.nbits();
    header.tsamp     = source.tsamp();
    header.tstart    = source.tstart();
    header.fch1      = source.fch1();
    header.foff      = source.foff();
    header.data_type = 1;
    header.src_raj   = source.raj();
    header.src_dej   = source.dej();
    header.ibeam     = source.beam();
    header.refdm     = source.refdm();
    if (const auto name = source.source_name(); !name.empty()) {
        header.source_name = std::string(name);
    }
    header.telescope_id =
        formats::sigproc::telescope_id_from_name(source.telescope());
    return header;
}

/**
 * @brief digifil output streamed from a BlockSource, without a file.
 *
 * Wraps an input BlockSource (any format, including a multi-file PSRFITS
 * stream) and runs the digifil Pipeline as samples are read. It models
 * `concepts::BlockReader`, so it can itself be wrapped in a BlockSource.
 * The stream is the SIGPROC payload digifil writes: nifs_out IFs of
 * nchans_out channels at nbits_out bits, channels descending in frequency,
 * sub-byte samples packed least significant bit first. header() is the
 * SIGPROC header digifil writes, with refdm kept.
 *
 * Reads run forward. While one input block is processed, the next is read on
 * another thread, so file reads overlap the rescaling. seek(), skip(), and
 * rewind() restart the pipeline at the matching input sample, so rescaling
 * statistics restart there as they would for a fresh run over that range.
 * Output samples are input samples / tscrunch.
 *
 * The wrapped source must not be used by anyone else while this object
 * reads it.
 */
class Source {
public:
    /// Read @p source, which this object then owns.
    /// @param start First input sample.
    /// @param nsamps Input samples to read; all after @p start by default.
    /// @throws ValidationError if the configuration does not fit the input.
    Source(BlockSource source,
           const Config& config,
           std::uint64_t start                 = 0,
           std::optional<std::uint64_t> nsamps = std::nullopt)
        : m_owned(std::make_unique<BlockSource>(std::move(source))),
          m_input(m_owned.get()),
          m_config(config) {
        init(header_of(*m_input), start, nsamps);
    }

    /// Read @p source, which must outlive this object. @p header is the full
    /// input header (FilterbankReader::header(), or one built for another
    /// format).
    Source(BlockSource& source,
           const formats::sigproc::FilterbankHeader& header,
           const Config& config,
           std::uint64_t start                 = 0,
           std::optional<std::uint64_t> nsamps = std::nullopt)
        : m_input(&source),
          m_config(config) {
        init(header, start, nsamps);
    }

    Source(const Source&)            = delete;
    Source& operator=(const Source&) = delete;
    Source(Source&&) noexcept        = default;
    // A read may be in flight into this object's buffers.
    Source& operator=(Source&&) = delete;
    ~Source() { wait(); }

    /// SIGPROC header of the output stream.
    [[nodiscard]] const formats::sigproc::FilterbankHeader&
    header() const noexcept {
        return m_header;
    }
    [[nodiscard]] const Pipeline& pipeline() const noexcept {
        return *m_pipeline;
    }
    /// Input samples per block pushed through the pipeline.
    [[nodiscard]] std::uint64_t gulp() const noexcept { return m_gulp; }
    /// Input samples read so far.
    [[nodiscard]] std::uint64_t input_samples_read() const noexcept {
        return m_consumed;
    }

    [[nodiscard]] std::optional<double> refdm() const noexcept {
        return m_header.refdm;
    }
    [[nodiscard]] std::uint64_t nchans() const noexcept {
        return m_pipeline->nchans_out();
    }
    [[nodiscard]] std::uint64_t nifs() const noexcept {
        return m_pipeline->nifs_out();
    }
    [[nodiscard]] int nbits() const noexcept { return m_pipeline->nbits_out(); }
    [[nodiscard]] SampleType sample_type() const {
        return sample_type_from_nbits(nbits(), false);
    }
    [[nodiscard]] std::uint64_t bytes_per_sample() const noexcept {
        return m_pipeline->bytes_per_sample();
    }
    [[nodiscard]] std::uint64_t nsamples() const noexcept { return m_nsamples; }
    [[nodiscard]] static bool has_nsamples() noexcept { return true; }
    [[nodiscard]] double tsamp() const { return m_header.tsamp.value_or(0.0); }
    [[nodiscard]] double tstart() const {
        return m_header.tstart.value_or(0.0);
    }
    [[nodiscard]] std::string_view source_name() const {
        return m_input->source_name();
    }
    [[nodiscard]] std::string_view telescope() const {
        return m_input->telescope();
    }
    [[nodiscard]] double raj() const { return m_input->raj(); }
    [[nodiscard]] double dej() const { return m_input->dej(); }
    [[nodiscard]] double fch1() const { return m_header.fch1.value_or(0.0); }
    [[nodiscard]] double foff() const { return m_header.foff.value_or(0.0); }
    [[nodiscard]] double bandwidth() const {
        return std::abs(foff()) * static_cast<double>(nchans());
    }
    [[nodiscard]] double center_frequency() const {
        return fch1() + (foff() * (static_cast<double>(nchans()) - 1.0) / 2.0);
    }
    [[nodiscard]] int beam() const { return m_input->beam(); }
    [[nodiscard]] double spectra_rate() const {
        const double interval = tsamp();
        return interval > 0.0 ? 1.0 / interval : 0.0;
    }
    [[nodiscard]] std::time_t utc_start() const {
        return astro::mjd_to_time(tstart());
    }
    [[nodiscard]] std::uint64_t tell() const noexcept { return m_sample; }

    /// Reverse the output channels (to ascending frequency) on read_block()
    /// and read_samples().
    void set_fswap(bool enable) noexcept { m_fswap = enable; }
    [[nodiscard]] bool fswap_enabled() const noexcept { return m_fswap; }

    /// Restart the pipeline so the next output sample is @p sample.
    /// @throws ValidationError if @p sample is past nsamples().
    void seek(std::uint64_t sample) {
        if (sample > m_nsamples) {
            throw ValidationError("psrio: seek is past the readable samples");
        }
        wait();
        const auto tfac = m_pipeline->tscrunch();
        const auto skip = sample * tfac;
        restart(m_start + skip, m_count - skip);
        m_sample = sample;
    }
    void rewind() { seek(0); }
    void skip(std::int64_t delta) {
        seek(::psrio::detail::apply_skip(m_sample, delta, m_nsamples));
    }

    /// The next encoded output, as the pipeline produces it, valid until the
    /// next call. Empty at the end. Does not honour fswap.
    std::span<const std::byte> next() {
        if (m_ready_at < m_ready.size()) {
            const auto rest =
                std::span<const std::byte>(m_ready).subspan(m_ready_at);
            m_ready_at = m_ready.size();
            m_sample += rest.size() / bytes_per_sample();
            return rest;
        }
        while (true) {
            const auto got = take_block();
            if (got == 0U) {
                return {};
            }
            const auto bytes = m_pipeline->push(m_current, got);
            if (!bytes.empty()) {
                m_sample += bytes.size() / bytes_per_sample();
                return bytes;
            }
        }
    }

    /// Copy up to @p count output samples into @p dest.
    /// @return Samples copied; short only at the end.
    std::uint64_t read_block(std::uint64_t count, std::span<std::byte> dest) {
        const auto stride = bytes_per_sample();
        if (dest.size() / stride < count) {
            throw ValidationError(fmt::format(
                "psrio: destination has {} bytes but the block needs {}",
                dest.size(), count * stride));
        }
        std::uint64_t done = 0;
        while (done < count) {
            if (m_ready_at == m_ready.size() && !refill()) {
                break;
            }
            const auto have = (m_ready.size() - m_ready_at) / stride;
            const auto take = std::min<std::uint64_t>(have, count - done);
            std::memcpy(dest.data() + (done * stride),
                        m_ready.data() + m_ready_at,
                        static_cast<std::size_t>(take * stride));
            m_ready_at += static_cast<std::size_t>(take * stride);
            done += take;
        }
        m_sample += done;
        if (m_fswap && done > 0U) {
            reverse_channels(
                dest.first(static_cast<std::size_t>(done * stride)),
                done * nifs(), nchans(), nbits());
        }
        return done;
    }

    /// Convenience allocating read of up to @p count output samples.
    [[nodiscard]] std::vector<std::byte> read_block(std::uint64_t count) {
        std::vector<std::byte> out(
            static_cast<std::size_t>(count * bytes_per_sample()));
        const auto got = read_block(count, out);
        out.resize(static_cast<std::size_t>(got * bytes_per_sample()));
        return out;
    }

    /// Unpack up to @p count output samples into @p dest, which holds
    /// count * nifs * nchans values. T is float for any nbits, uint8_t for 8
    /// bits or fewer, uint16_t for 16 bits.
    template <typename T>
        requires std::same_as<T, float> || std::same_as<T, std::uint8_t> ||
                 std::same_as<T, std::uint16_t>
    std::uint64_t read_samples(std::uint64_t count, std::span<T> dest) {
        const auto per_sample = nifs() * nchans();
        if (dest.size() != count * per_sample) {
            throw ValidationError(fmt::format(
                "psrio: destination has {} values but the request needs {}",
                dest.size(), count * per_sample));
        }
        const int bits = nbits();
        if constexpr (std::same_as<T, std::uint8_t>) {
            if (bits > 8) {
                throw ValidationError(
                    "psrio: uint8_t output requires nbits of 8 or less");
            }
        } else if constexpr (std::same_as<T, std::uint16_t>) {
            if (bits != 16) {
                throw ValidationError(
                    "psrio: uint16_t output requires 16-bit samples");
            }
        }
        m_unpacked.resize(static_cast<std::size_t>(count * bytes_per_sample()));
        const bool fswap = m_fswap;
        m_fswap          = false;
        const auto got   = read_block(count, m_unpacked);
        m_fswap          = fswap;
        const auto raw =
            std::span<const std::byte>(m_unpacked)
                .first(static_cast<std::size_t>(got * bytes_per_sample()));
        const auto out = dest.first(static_cast<std::size_t>(got * per_sample));
        if (bits < 8) {
            if constexpr (!std::same_as<T, std::uint16_t>) {
                ::psrio::detail::unpack_sub_byte(raw, out, bits,
                                                 BitOrder::kLsbFirst);
            }
        } else if (bits == 8) {
            if constexpr (!std::same_as<T, std::uint16_t>) {
                ::psrio::detail::unpack_8bit(raw, out, false);
            }
        } else if (bits == 16) {
            if constexpr (!std::same_as<T, std::uint8_t>) {
                ::psrio::detail::unpack_16le(raw, out);
            }
        } else if constexpr (std::same_as<T, float>) {
            ::psrio::detail::unpack_32le(raw, out);
        } else {
            throw ValidationError("psrio: 32-bit samples unpack to float");
        }
        if (m_fswap && got > 0U) {
            reverse_unpacked_channels(out, got, nchans(), nifs());
        }
        return got;
    }

    /// Copy the next @p nbytes output bytes, a whole number of samples.
    std::uint64_t read_bytes(std::uint64_t nbytes, std::span<std::byte> dest) {
        const auto stride = bytes_per_sample();
        if (nbytes == 0U || nbytes % stride != 0U || dest.size() != nbytes) {
            throw ValidationError("psrio: byte request must be a positive "
                                  "multiple of the sample stride");
        }
        if (nbytes / stride > m_nsamples - m_sample) {
            throw ValidationError(
                "psrio: byte request extends past the readable samples");
        }
        const bool fswap = m_fswap;
        m_fswap          = false;
        const auto got   = read_block(nbytes / stride, dest);
        m_fswap          = fswap;
        return got * stride;
    }

private:
    void init(const formats::sigproc::FilterbankHeader& header,
              std::uint64_t start,
              std::optional<std::uint64_t> nsamps) {
        const auto total = m_input->nsamples();
        if (start >= total) {
            throw ValidationError(fmt::format(
                "psrio: start sample {} is outside the input ({} samples)",
                start, total));
        }
        const auto count = nsamps.value_or(total - start);
        if (count == 0U || count > total - start) {
            throw ValidationError(fmt::format(
                "psrio: {} samples from {} exceed the input ({} samples)",
                count, start, total));
        }
        m_start    = start;
        m_count    = count;
        m_pipeline = std::make_unique<Pipeline>(m_config, layout());
        m_gulp     = m_config.gulp.value_or(
            default_gulp(*m_input, *m_pipeline, m_config.max_ram_bytes()));
        m_header   = output_header(header, *m_pipeline, start);
        m_nsamples = count / m_pipeline->tscrunch();
        m_input->set_fswap(false);
        restart(start, count);
    }

    [[nodiscard]] InputLayout layout() const {
        return InputLayout{
            .nchans = m_input->nchans(),
            .nifs   = m_input->nifs(),
            .nbits  = m_input->nbits(),
            .tsamp  = m_input->tsamp(),
            .foff   = m_input->foff(),
        };
    }

    /// Start reading @p count input samples at @p start with a new pipeline.
    void restart(std::uint64_t start, std::uint64_t count) {
        wait();
        if (m_consumed != 0U || m_requested != 0U) {
            m_pipeline = std::make_unique<Pipeline>(m_config, layout());
        }
        m_consumed  = 0;
        m_requested = 0;
        m_ready.clear();
        m_ready_at = 0;
        m_input->seek(start);
        m_count_left = count;
        launch();
    }

    /// Read the next input block on another thread.
    void launch() {
        const auto row  = m_input->nifs() * m_input->nchans();
        const auto want = std::min(m_gulp, m_count_left - m_requested);
        if (want == 0U) {
            return;
        }
        m_next.resize(static_cast<std::size_t>(want * row));
        m_requested += want;
        m_pending =
            std::async(std::launch::async, [input = m_input, want,
                                            block = std::span<float>(m_next)] {
                return input->read_samples(want, block);
            });
    }

    /// Wait for the block in flight, start the one after, and return how
    /// many samples m_current now holds.
    std::uint64_t take_block() {
        if (!m_pending.valid()) {
            return 0;
        }
        const auto got = m_pending.get();
        std::swap(m_current, m_next);
        m_consumed += got;
        if (got > 0U) {
            launch();
        }
        return got;
    }

    bool refill() {
        const auto bytes = next_from_pipeline();
        if (bytes.empty()) {
            return false;
        }
        m_ready.assign(bytes.begin(), bytes.end());
        m_ready_at = 0;
        return true;
    }

    std::span<const std::byte> next_from_pipeline() {
        while (true) {
            const auto got = take_block();
            if (got == 0U) {
                return {};
            }
            const auto bytes = m_pipeline->push(m_current, got);
            if (!bytes.empty()) {
                return bytes;
            }
        }
    }

    void wait() noexcept {
        if (m_pending.valid()) {
            m_pending.wait();
        }
    }

    std::unique_ptr<BlockSource> m_owned;
    BlockSource* m_input{nullptr};
    Config m_config;
    formats::sigproc::FilterbankHeader m_header;
    std::unique_ptr<Pipeline> m_pipeline;
    std::uint64_t m_start{0};
    std::uint64_t m_count{0};
    std::uint64_t m_count_left{0};
    std::uint64_t m_gulp{0};
    std::uint64_t m_consumed{0};
    std::uint64_t m_requested{0};
    std::uint64_t m_sample{0};
    std::uint64_t m_nsamples{0};
    bool m_fswap{false};
    std::vector<std::byte> m_ready;
    std::size_t m_ready_at{0};
    std::vector<float> m_current;
    std::vector<float> m_next;
    std::vector<std::byte> m_unpacked;
    // Last, so it is destroyed (and joined) before the buffers it fills.
    std::future<std::uint64_t> m_pending;
};

/// Counts reported by run().
struct Summary {
    std::uint64_t nsamples_in{0};
    std::uint64_t nsamples_out{0};
    std::uint64_t gulp{0};
};

/// Stream samples [start, start + nsamps) of @p source through digifil into
/// a new SIGPROC file at @p output.
///
/// @p header is the full input header (FilterbankReader::header(), or one
/// built for another format); @p source supplies the samples as float in
/// file channel order. Channel reversal on @p source is turned off.
/// @throws ValidationError if the configuration or the sample range does not
///         fit the input.
/// @throws IoError if the output cannot be written.
inline Summary run(BlockSource& source,
                   const formats::sigproc::FilterbankHeader& header,
                   const std::filesystem::path& output,
                   const Config& config,
                   std::uint64_t start                 = 0,
                   std::optional<std::uint64_t> nsamps = std::nullopt) {
    Source stream(source, header, config, start, nsamps);
    formats::sigproc::FilterbankWriter writer(output, stream.header());
    while (true) {
        const auto bytes = stream.next();
        if (bytes.empty()) {
            break;
        }
        writer.write_bytes(bytes);
    }
    writer.close();
    return Summary{
        .nsamples_in  = stream.input_samples_read(),
        .nsamples_out = writer.nsamples_written(),
        .gulp         = stream.gulp(),
    };
}

} // namespace psrio::digifil
