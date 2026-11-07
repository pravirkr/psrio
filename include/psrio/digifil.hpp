#pragma once

#include "psrio/block_source.hpp"
#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/packed_bits.hpp"
#include "psrio/formats/sigproc/header.hpp"
#include "psrio/formats/sigproc/writer.hpp"
#include "psrio/scrunch.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fmt/format.h>
#include <limits>
#include <optional>
#include <span>
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
    /// Time decimation factor; 0 means 1 (-t).
    std::uint64_t tscrunch{1};
    /// Frequency decimation factor; 0 means 1 (-f).
    std::uint64_t fscrunch{1};
    /// Extra scale applied before digitising (-s).
    float scale_fac{1.0F};
    /// Keep only this input IF, or all IFs when negative (-P).
    int poln_select{-1};
    /// Input samples per block. Rescaling statistics depend on the block
    /// boundaries during the first interval. Defaults to dspsr's block for a
    /// 256 MB RAM limit (default_gulp()).
    std::optional<std::uint64_t> gulp;

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
        if (!std::isfinite(scale_fac) || scale_fac == 0.0F) {
            throw ValidationError(fmt::format(
                "psrio: digifil scale factor must be finite and non-zero, "
                "got {}",
                scale_fac));
        }
        if (gulp.has_value() && *gulp == 0U) {
            throw ValidationError("psrio: digifil gulp must be positive");
        }
    }
};

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

/// Per-column rescaling to zero mean and unit variance, as dspsr `Rescale`.
///
/// Sums and sums of squares accumulate in double over intervals of
/// nsample rows, which may span blocks. During the first interval the
/// statistics are recomputed at the end of every block, so block boundaries
/// matter there. Each row is rescaled as (x + offset) * scale with the latest
/// statistics.
class Rescaler {
public:
    /// @throws ValidationError if @p ncols or @p nsample is 0.
    Rescaler(std::uint64_t ncols, std::uint64_t nsample, bool constant)
        : m_nsample(nsample),
          m_constant(constant),
          m_sums(static_cast<std::size_t>(ncols), 0.0),
          m_sumsq(static_cast<std::size_t>(ncols), 0.0),
          m_offset(static_cast<std::size_t>(ncols), 0.0F),
          m_scale(static_cast<std::size_t>(ncols), 1.0F) {
        if (ncols == 0U || nsample == 0U) {
            throw ValidationError(
                "psrio: rescaler needs at least one column and sample");
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
        const auto ncols = m_sums.size();
        for (std::size_t base = 0; base < rows.size(); base += ncols) {
            for (std::size_t col = 0; col < ncols; ++col) {
                const auto value = static_cast<double>(rows[base + col]);
                m_sums[col] += value;
                m_sumsq[col] += value * value;
            }
        }
    }

    void compute() {
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

    void apply(std::span<float> rows) const {
        const auto ncols = m_offset.size();
        for (std::size_t base = 0; base < rows.size(); base += ncols) {
            for (std::size_t col = 0; col < ncols; ++col) {
                const float shifted = rows[base + col] + m_offset[col];
                rows[base + col]    = shifted * m_scale[col];
            }
        }
    }

    std::uint64_t m_nsample;
    bool m_constant;
    bool m_first{true};
    std::uint64_t m_isample{0};
    std::vector<double> m_sums;
    std::vector<double> m_sumsq;
    std::vector<float> m_offset;
    std::vector<float> m_scale;
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
                               config.rescale_constant);
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

        if (m_is_float) {
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
    const auto total = source.nsamples();
    if (start >= total) {
        throw ValidationError(fmt::format(
            "psrio: start sample {} is outside the input ({} samples)", start,
            total));
    }
    const auto count = nsamps.value_or(total - start);
    if (count == 0U || count > total - start) {
        throw ValidationError(fmt::format(
            "psrio: {} samples from {} exceed the input ({} samples)", count,
            start, total));
    }

    Pipeline pipeline(config, InputLayout{
                                  .nchans = source.nchans(),
                                  .nifs   = source.nifs(),
                                  .nbits  = source.nbits(),
                                  .tsamp  = source.tsamp(),
                                  .foff   = source.foff(),
                              });
    const auto gulp = config.gulp.value_or(block_samples(
        source.nbits(), source.nifs(), source.nchans(), pipeline.float_copies(),
        kSigprocResolution, kDefaultMaxRamBytes));

    formats::sigproc::FilterbankWriter writer(
        output, output_header(header, pipeline, start));

    source.set_fswap(false);
    source.seek(start);
    const auto row = source.nifs() * source.nchans();
    std::vector<float> block(
        static_cast<std::size_t>(std::min(gulp, count) * row));
    Summary summary{.nsamples_in = 0, .nsamples_out = 0, .gulp = gulp};
    while (summary.nsamples_in < count) {
        const auto want = std::min(gulp, count - summary.nsamples_in);
        const auto got  = source.read_samples(
            want, std::span<float>(block).first(
                      static_cast<std::size_t>(want * row)));
        if (got == 0U) {
            break;
        }
        const auto bytes = pipeline.push(block, got);
        writer.write_bytes(bytes);
        summary.nsamples_in += got;
    }
    writer.close();
    summary.nsamples_out = writer.nsamples_written();
    return summary;
}

} // namespace psrio::digifil
