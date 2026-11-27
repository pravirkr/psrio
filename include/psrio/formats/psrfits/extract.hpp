#pragma once

#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/mjd.hpp"
#include "psrio/detail/parallel.hpp"
#include "psrio/formats/psrfits/archive.hpp"
#include "psrio/formats/psrfits/stream.hpp"
#include "psrio/formats/psrfits/table.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace psrio::formats::psrfits {

/// dspsr's and psrchive's dispersion constant: a delay is
/// `DM / 2.41e-4 * (1/f1^2 - 1/f2^2)` seconds with frequencies in MHz.
inline constexpr double kDspsrDispersion = 2.41e-4;

/**
 * @brief Band and per-channel dispersion delays, computed as dspsr does.
 *
 * Channel centres come from OBSFREQ, OBSBW, and NCHAN
 * (`dsp::Observation::get_centre_frequency`); DAT_FREQ is not used.
 * Delays follow `dsp::Dedispersion::SampleDelay` and `dsp::SampleDelay`:
 * `D_c = floor(DM/2.41e-4 * (1/f0^2 - 1/fc^2) * rate + 0.5)` relative to the
 * centre frequency f0, a span centre delay from the mean of the first and
 * last channel frequencies, and `applied_c = span + zero - D_c >= 0`, the
 * number of samples channel c is advanced by.
 */
struct DspsrBand {
    double centre_frequency{0.0};
    double bandwidth{0.0};
    double rate{0.0};
    double dm{0.0};
    std::vector<double> frequencies;
    /// Samples each channel is advanced by (all zero without dedispersion).
    std::vector<std::int64_t> applied;
    /// dspsr's `zero_delay`: output time stamps move by zero/rate.
    std::int64_t zero_delay{0};
    /// Delay of the span centre relative to f0 (normally 0).
    std::int64_t span_centre{0};
    /// Largest |D_c| over the band, in samples.
    std::int64_t max_abs_delay{0};
    bool dedisperse{false};

    [[nodiscard]] std::int64_t max_applied() const {
        return applied.empty() ? 0 : *std::ranges::max_element(applied);
    }
    [[nodiscard]] std::int64_t min_applied() const {
        return applied.empty() ? 0 : *std::ranges::min_element(applied);
    }
};

/// dspsr channel centre frequencies (MHz) from the band centre and width.
[[nodiscard]] inline std::vector<double> dspsr_channel_frequencies(
    double centre, double bandwidth, std::uint64_t nchan) {
    std::vector<double> freqs(static_cast<std::size_t>(nchan));
    const double base = centre - (0.5 * bandwidth) +
                        (0.5 * bandwidth / static_cast<double>(nchan));
    for (std::uint64_t chan = 0; chan < nchan; ++chan) {
        const auto channel = static_cast<double>(chan);
        freqs[static_cast<std::size_t>(chan)] =
            base + (channel * bandwidth / static_cast<double>(nchan));
    }
    return freqs;
}

/// Seconds a signal at @p from_mhz arrives before one at @p to_mhz
/// (positive when @p from_mhz is the higher frequency), dspsr's constant.
[[nodiscard]] inline double
dispersion_offset(double dm, double from_mhz, double to_mhz) {
    return (dm / kDspsrDispersion) *
           ((1.0 / (to_mhz * to_mhz)) - (1.0 / (from_mhz * from_mhz)));
}

/// dispersion_offset() from @p freq to @p centre (MHz) in samples, rounded
/// to nearest with halves up: the per-channel shift dspsr applies.
[[nodiscard]] inline std::int64_t
dspsr_sample_delay(double dm, double centre, double freq, double rate) {
    const double seconds = dispersion_offset(dm, freq, centre);
    return static_cast<std::int64_t>(std::floor((seconds * rate) + 0.5));
}

[[nodiscard]] inline DspsrBand make_dspsr_band(double centre,
                                               double bandwidth,
                                               std::uint64_t nchan,
                                               double rate,
                                               double dm,
                                               bool dedisperse) {
    if (nchan == 0U || !(rate > 0.0)) {
        throw ValidationError("psrio: a band needs channels and a sample rate");
    }
    DspsrBand band;
    band.centre_frequency = centre;
    band.bandwidth        = bandwidth;
    band.rate             = rate;
    band.dm               = dm;
    band.dedisperse       = dedisperse;
    band.frequencies      = dspsr_channel_frequencies(centre, bandwidth, nchan);
    std::vector<std::int64_t> delays(band.frequencies.size());
    for (std::size_t chan = 0; chan < delays.size(); ++chan) {
        delays[chan] =
            dspsr_sample_delay(dm, centre, band.frequencies[chan], rate);
        band.max_abs_delay =
            std::max(band.max_abs_delay,
                     delays[chan] < 0 ? -delays[chan] : delays[chan]);
    }
    band.applied.assign(delays.size(), 0);
    if (!dedisperse) {
        return band;
    }
    const double mid = (band.frequencies.front() + band.frequencies.back()) / 2;
    const double span_delay = dispersion_offset(dm, mid, centre);
    band.span_centre =
        static_cast<std::int64_t>(std::floor((span_delay * rate) + 0.5));
    std::int64_t zero = 0;
    for (const auto delay : delays) {
        zero = std::max(zero, delay - band.span_centre);
    }
    band.zero_delay = zero;
    for (std::size_t chan = 0; chan < delays.size(); ++chan) {
        band.applied[chan] = band.span_centre + zero - delays[chan];
    }
    return band;
}

/**
 * @brief One folded sub-integration taken from a search-mode stream.
 *
 * Output sample `o` (counted from the start of the stream) of channel c is
 * input sample `o + applied_c`. The window is output samples
 * `[first, first + nsamples)`, folded into `nbin` bins of
 * `nsamples / nbin` consecutive samples, so bin 0 starts the window.
 */
struct ExtractWindow {
    std::uint64_t first{0};
    std::uint64_t nsamples{0};
    std::uint32_t nbin{0};

    [[nodiscard]] std::uint64_t samples_per_bin() const {
        return nbin == 0U ? 0U : nsamples / nbin;
    }
};

/// Input samples [first, last) a window reads.
struct InputSpan {
    std::uint64_t first{0};
    std::uint64_t last{0};
};

[[nodiscard]] inline InputSpan input_span(const ExtractWindow& window,
                                          const DspsrBand& band) {
    return InputSpan{
        .first = window.first + static_cast<std::uint64_t>(band.min_applied()),
        .last  = window.first + window.nsamples +
                 static_cast<std::uint64_t>(band.max_applied()),
    };
}

inline void check_window(const ExtractWindow& window,
                         const DspsrBand& band,
                         std::uint64_t stream_samples) {
    if (window.nbin == 0U || window.nsamples == 0U) {
        throw ValidationError("psrio: the window needs samples and bins");
    }
    if (window.nsamples % window.nbin != 0U) {
        throw ValidationError(
            fmt::format("psrio: {} samples do not divide into {} bins",
                        window.nsamples, window.nbin));
    }
    const auto span = input_span(window, band);
    if (span.last > stream_samples) {
        throw ValidationError(fmt::format(
            "psrio: the window needs input samples {}..{} but the data end at "
            "{}",
            span.first, span.last, stream_samples));
    }
}

/// Parse a decimal MJD string without losing the fraction to a double.
[[nodiscard]] inline ::psrio::detail::Mjd parse_mjd(std::string_view text) {
    text             = detail::trim(text);
    const auto dot   = text.find('.');
    const auto ipart = text.substr(0, dot);
    int days         = 0;
    const auto [ptr, ec] =
        std::from_chars(ipart.data(), ipart.data() + ipart.size(), days);
    if (ec != std::errc{} || ptr != ipart.data() + ipart.size() || days < 0) {
        throw ValidationError(fmt::format("psrio: '{}' is not an MJD", text));
    }
    double frac = 0.0;
    if (dot != std::string_view::npos) {
        const auto digits = std::string{"0"} + std::string{text.substr(dot)};
        const auto* begin = digits.data();
        const auto* end   = digits.data() + digits.size();
        const auto [fptr, fec] = std::from_chars(begin, end, frac);
        if (fec != std::errc{} || fptr != end) {
            throw ValidationError(
                fmt::format("psrio: '{}' is not an MJD", text));
        }
    }
    return ::psrio::detail::Mjd(days, 0, 0.0) + (frac * 86400.0);
}

/// MJD as text with 15 decimal places in the day fraction.
[[nodiscard]] inline std::string format_mjd(const ::psrio::detail::Mjd& mjd) {
    const double frac =
        (static_cast<double>(mjd.second()) + mjd.fraction()) / 86400.0;
    const auto text = fmt::format("{:.15f}", frac);
    if (text.starts_with("1.")) {
        return fmt::format("{}.000000000000000", mjd.day() + 1);
    }
    return fmt::format("{}{}", mjd.day(), text.substr(1));
}

/**
 * @brief Times stored for a window, following dspsr's arithmetic.
 *
 * dspsr stamps the block that holds the window as
 * `STT + load/rate (+ zero_delay/rate)`, with `load` the start of the
 * SUBINT row at or before the window, and the sub-integration as that plus
 * `(first - load)/rate`. The epoch is the mid-time pulled back to a multiple
 * of the period counted from MJD 0 (`dsp::PhaseSeries::get_mid_time`).
 */
struct WindowTiming {
    ::psrio::detail::Mjd start;
    double period{0.0};
    double tsubint{0.0};
    double offs_sub{0.0};
};

[[nodiscard]] inline WindowTiming
window_timing(const ::psrio::detail::Mjd& stream_start,
              std::uint64_t nsblk,
              const ExtractWindow& window,
              const DspsrBand& band) {
    const double rate = band.rate;
    const auto load   = (window.first / nsblk) * nsblk;
    auto block        = stream_start + (static_cast<double>(load) / rate);
    if (band.dedisperse) {
        block = block + (static_cast<double>(band.zero_delay) / rate);
    }
    WindowTiming out;
    out.start      = block + (static_cast<double>(window.first - load) / rate);
    const auto end = out.start + (static_cast<double>(window.nsamples) / rate);
    out.period     = static_cast<double>(window.nsamples) / rate;
    out.tsubint    = static_cast<double>(window.nsamples) / rate;
    auto mid       = ::psrio::detail::Mjd::midpoint(out.start, end);
    const double phase = std::fmod(mid.in_seconds(), out.period) / out.period;
    mid -= phase * out.period;
    out.offs_sub = mid.seconds_since(out.start);
    return out;
}

/// Time stamp of output sample @p sample (start of the sample).
[[nodiscard]] inline ::psrio::detail::Mjd
output_time(const ::psrio::detail::Mjd& stream_start,
            std::uint64_t sample,
            const DspsrBand& band) {
    return stream_start +
           (static_cast<double>(static_cast<std::int64_t>(sample) +
                                (band.dedisperse ? band.zero_delay : 0)) /
            band.rate);
}

/**
 * @brief Window of @p nsamples output samples centred on a burst.
 *
 * @p toa is the arrival time at @p toa_freq (MHz). It is moved to the band
 * centre with dspsr's constant, and the window is placed so that the
 * output sample holding that time is sample `nsamples / 2`, the first sample
 * of bin `nbin / 2`.
 * @throws ValidationError if the window would start before the data.
 */
[[nodiscard]] inline ExtractWindow
centred_window(const ::psrio::detail::Mjd& stream_start,
               const ::psrio::detail::Mjd& toa,
               double toa_freq,
               std::uint64_t nsamples,
               std::uint32_t nbin,
               const DspsrBand& band) {
    const double shift =
        dispersion_offset(band.dm, toa_freq, band.centre_frequency);
    const auto at_centre = toa + shift;
    const double seconds = at_centre.seconds_since(stream_start);
    const auto sample =
        static_cast<std::int64_t>(std::llround(seconds * band.rate));
    // Channel f0 of output o reads input o + span + zero.
    const std::int64_t output =
        sample - (band.dedisperse ? band.span_centre + band.zero_delay : 0);
    const std::int64_t first =
        output - static_cast<std::int64_t>(nsamples / 2U);
    if (first < 0 || first + band.min_applied() < 0) {
        throw ValidationError(fmt::format(
            "psrio: the window would start {:.6f} s before the data",
            -static_cast<double>(first) / band.rate));
    }
    return ExtractWindow{
        .first    = static_cast<std::uint64_t>(first),
        .nsamples = nsamples,
        .nbin     = nbin,
    };
}

/**
 * @brief Fold a window of every stored product into @p out.
 *
 * Reads the input once, in chunks of at most @p chunk_bytes. Each sample is
 * unpacked as dspsr does (SearchStream::read_products), summed in float in
 * time order per bin, and divided by the bin's sample count in double, as
 * `dsp::Archiver` does.
 *
 * @param out `[pol][chan][bin]` amplitudes, NPOL * NCHAN * nbin floats.
 */
inline void fold_window(SearchStream& stream,
                        const ExtractWindow& window,
                        const DspsrBand& band,
                        std::span<float> out,
                        std::uint64_t chunk_bytes = std::uint64_t{64} << 20U) {
    check_window(window, band, stream.nsamples());
    const auto nchan = static_cast<std::size_t>(stream.nchans());
    const auto npol  = static_cast<std::size_t>(stream.input_npol());
    const auto nbin  = static_cast<std::size_t>(window.nbin);
    if (band.applied.size() != nchan) {
        throw ValidationError("psrio: band and stream disagree on NCHAN");
    }
    if (out.size() != npol * nchan * nbin) {
        throw ValidationError(
            fmt::format("psrio: fold output has {} values, needs {}",
                        out.size(), npol * nchan * nbin));
    }
    std::ranges::fill(out, 0.0F);
    const auto per_sample = npol * nchan;
    const auto span       = input_span(window, band);
    const auto spb        = static_cast<std::int64_t>(window.samples_per_bin());
    const auto nsamples   = static_cast<std::int64_t>(window.nsamples);
    const auto chunk =
        std::max<std::uint64_t>(1U, chunk_bytes / (per_sample * sizeof(float)));
    std::vector<float> buffer;
    stream.seek(span.first);
    for (std::uint64_t at = span.first; at < span.last;) {
        const auto count = std::min(chunk, span.last - at);
        buffer.resize(static_cast<std::size_t>(count) * per_sample);
        const auto got = stream.read_products(count, buffer);
        if (got != count) {
            throw FormatError("psrio: the data ended inside the window");
        }
        const auto base   = static_cast<std::int64_t>(at) -
                            static_cast<std::int64_t>(window.first);
        const float* data = buffer.data();
        const auto ncount = static_cast<std::int64_t>(count);
        [[maybe_unused]] const bool wide =
            ::psrio::detail::parallel_worthwhile(count * per_sample);
        PSRIO_OMP(omp parallel for default(none) schedule(static) if (wide)
                      shared(nchan, npol, nbin, per_sample, band, spb, nsamples,
                             ncount, base, data, out))
        for (std::size_t chan = 0; chan < nchan; ++chan) {
            const auto shift = band.applied[chan];
            // Sample t of the chunk is output index base + t - shift.
            const auto lo = std::max<std::int64_t>(0, shift - base);
            const auto hi =
                std::min<std::int64_t>(ncount, nsamples + shift - base);
            for (std::size_t pol = 0; pol < npol; ++pol) {
                float* prof = out.data() + (((pol * nchan) + chan) * nbin);
                const float* input = data + (pol * nchan) + chan;
                for (std::int64_t t = lo; t < hi; ++t) {
                    const auto index = base + t - shift;
                    prof[static_cast<std::size_t>(index / spb)] +=
                        input[static_cast<std::size_t>(t) * per_sample];
                }
            }
        }
        at += count;
    }
    const auto hits = static_cast<double>(spb);
    for (auto& value : out) {
        value = static_cast<float>(static_cast<double>(value) / (1.0 * hits));
    }
}

/// psrchive state string for the archive's POL_TYPE.
[[nodiscard]] inline std::string_view archive_pol_type(std::string_view input,
                                                       std::uint64_t npol) {
    if (npol == 1U) {
        return "INTEN";
    }
    if (npol == 2U) {
        return "AABB";
    }
    if (input == "IQUV" || input == "STOKES" || input == "STOKE") {
        return "IQUV";
    }
    return "AABBCRCI";
}

} // namespace psrio::formats::psrfits
