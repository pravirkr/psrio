#include "common.hpp"

#include <CLI/CLI.hpp>
#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fmt/format.h>
#include <limits>
#include <memory>
#include <optional>
#include <psrio/detail/exceptions.hpp>
#include <psrio/detail/mjd.hpp>
#include <psrio/formats/psrfits.hpp>
#include <psrio/formats/psrfits/archive.hpp>
#include <psrio/formats/psrfits/extract.hpp>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace psrio::cli {

namespace {

namespace psrfits = ::psrio::formats::psrfits;

struct CutoutOptions {
    std::vector<std::filesystem::path> inputs;
    std::filesystem::path output;
    std::string extension{"ar"};
    double dm{0.0};
    bool dedisperse{false};
    std::uint32_t nbin{0};
    bool scloffs{false};
    std::optional<std::string> toa;
    std::string toa_freq{"top"};
    std::optional<double> start_seconds;
    std::optional<double> length_seconds;
    double margin{0.1};
};

/// Frequency named by --toa-freq: top, bottom, centre, inf, or MHz.
double toa_frequency(std::string_view text, const psrfits::DspsrBand& band) {
    const auto [low, high] = std::ranges::minmax(band.frequencies);
    if (text == "top") {
        return high;
    }
    if (text == "bottom") {
        return low;
    }
    if (text == "centre" || text == "center") {
        return band.centre_frequency;
    }
    if (text == "inf") {
        return std::numeric_limits<double>::infinity();
    }
    double mhz = 0.0;
    const auto [ptr, ec] =
        std::from_chars(text.data(), text.data() + text.size(), mhz);
    if (ec != std::errc{} || ptr != text.data() + text.size() || !(mhz > 0.0)) {
        throw ValidationError(fmt::format(
            "psrio: --toa-freq '{}' is not top, bottom, centre, inf, or a "
            "frequency in MHz",
            text));
    }
    return mhz;
}

std::string command_line(const CutoutOptions& opts) {
    std::string cmd =
        fmt::format("psrio cutout -D {} -b {}", opts.dm, opts.nbin);
    if (opts.dedisperse) {
        cmd += " -K";
    }
    if (opts.scloffs) {
        cmd += " --scloffs";
    }
    if (opts.toa.has_value()) {
        cmd += fmt::format(" --toa {} --toa-freq {}", *opts.toa, opts.toa_freq);
    }
    if (opts.start_seconds.has_value()) {
        cmd += fmt::format(" -S {}", *opts.start_seconds);
    }
    if (opts.length_seconds.has_value()) {
        cmd += fmt::format(" -T {}", *opts.length_seconds);
    }
    for (const auto& path : opts.inputs) {
        cmd += ' ';
        cmd += path.filename().string();
    }
    return cmd;
}

void run_cutout(const CutoutOptions& opts, const Context& ctx) {
    if (opts.toa.has_value() == opts.start_seconds.has_value()) {
        throw ValidationError("psrio: give either --toa or -S, not both");
    }
    if (opts.nbin == 0U) {
        throw ValidationError("psrio: -b must be at least 1");
    }
    for (const auto& path : opts.inputs) {
        if (!psrfits::is_search_path(path)) {
            throw ValidationError(fmt::format(
                "psrio: cutout reads search-mode PSRFITS (.sf, .fits), not "
                "'{}'",
                path.string()));
        }
    }
    psrfits::SearchStream stream(opts.inputs, psrfits::DecodeOptions{
                                                  .scloffs = opts.scloffs,
                                                  .weights = false,
                                                  .product = -1,
                                              });
    const auto& first   = stream.file(0);
    const auto& primary = first.table().primary;
    const auto backend  = std::string(first.backend());
    if (opts.scloffs && (backend == "GUPPI" || backend == "PUPPI" ||
                         backend == "VEGAS" || backend == "COBALT")) {
        ctx.warn("psrio cutout: dspsr ignores -scloffs for {} data; this "
                 "output applies DAT_SCL and DAT_OFFS",
                 backend);
    }
    const ::psrio::detail::Mjd stream_start(primary.integer("STT_IMJD"),
                                            primary.integer("STT_SMJD"),
                                            primary.number("STT_OFFS"));
    const double rate = 1.0 / stream.tsamp();
    const auto band   = psrfits::make_dspsr_band(
        primary.number("OBSFREQ"), primary.number("OBSBW"), stream.nchans(),
        rate, opts.dm, opts.dedisperse);

    // Window length: whole samples per bin.
    std::uint64_t wanted = 0;
    if (opts.length_seconds.has_value()) {
        if (!(*opts.length_seconds > 0.0)) {
            throw ValidationError("psrio: -T must be positive");
        }
        wanted = static_cast<std::uint64_t>(
            std::llround(*opts.length_seconds * rate));
    } else if (opts.dedisperse) {
        throw ValidationError("psrio: -K needs the window length -T");
    } else {
        // The sweep on either side of the band centre, plus a margin.
        wanted = static_cast<std::uint64_t>(
            std::ceil(2.0 * static_cast<double>(band.max_abs_delay) *
                      (1.0 + opts.margin)));
    }
    wanted              = std::max<std::uint64_t>(wanted, opts.nbin);
    const auto per_bin  = (wanted + opts.nbin - 1U) / opts.nbin;
    const auto nsamples = per_bin * opts.nbin;

    psrfits::ExtractWindow window{};
    if (opts.toa.has_value()) {
        const auto toa  = psrfits::parse_mjd(*opts.toa);
        const auto freq = toa_frequency(opts.toa_freq, band);
        window = psrfits::centred_window(stream_start, toa, freq, nsamples,
                                         opts.nbin, band);
    } else {
        // -S is the band-centre arrival time of the window start. With -K,
        // output sample o holds that time for input o + span + zero.
        const auto at_centre = std::llround(*opts.start_seconds * rate);
        const auto first =
            at_centre -
            (opts.dedisperse ? band.span_centre + band.zero_delay : 0);
        if (*opts.start_seconds < 0.0 || first < 0) {
            throw ValidationError(fmt::format(
                "psrio: -S {} s starts the window before the data (the "
                "earliest with -K is {:.6f} s)",
                *opts.start_seconds,
                static_cast<double>(band.span_centre + band.zero_delay) /
                    rate));
        }
        window = psrfits::ExtractWindow{
            .first    = static_cast<std::uint64_t>(first),
            .nsamples = nsamples,
            .nbin     = opts.nbin,
        };
    }
    psrfits::check_window(window, band, stream.nsamples());

    const auto timing = psrfits::window_timing(
        stream_start, stream.subint_samples(), window, band);
    const auto span = psrfits::input_span(window, band);
    ctx.info("psrio cutout: window starts at MJD {} ({} samples, {} per bin, "
             "P = {:.9g} s)",
             psrfits::format_mjd(timing.start), nsamples, per_bin,
             timing.period);
    ctx.debug("psrio cutout: output samples {}..{}, input samples {}..{}, "
              "zero delay {}, largest shift {}",
              window.first, window.first + nsamples, span.first, span.last,
              band.zero_delay, band.max_applied());

    const auto npol  = stream.input_npol();
    const auto nchan = stream.nchans();
    std::vector<float> amps(static_cast<std::size_t>(npol * nchan * opts.nbin));
    psrfits::fold_window(stream, window, band, amps);

    psrfits::FoldArchiveSpec spec;
    spec.input_primary    = primary;
    spec.input_subint     = first.table().subint;
    spec.start            = timing.start;
    spec.period           = timing.period;
    spec.tsubint          = timing.tsubint;
    spec.offs_sub         = timing.offs_sub;
    spec.nbin             = opts.nbin;
    spec.nchan            = static_cast<std::uint32_t>(nchan);
    spec.npol             = static_cast<std::uint32_t>(npol);
    spec.pol_type         = psrfits::archive_pol_type(stream.pol_type(), npol);
    spec.dm               = opts.dm;
    spec.dedispersed      = opts.dedisperse;
    spec.centre_frequency = band.centre_frequency;
    spec.bandwidth        = band.bandwidth;
    spec.frequencies      = band.frequencies;
    spec.command          = command_line(opts);
    spec.created          = std::time(nullptr);

    std::filesystem::path output = opts.output;
    if (output.empty()) {
        auto mjd = psrfits::format_mjd(timing.start);
        mjd.resize(mjd.size() - 2U);
        output = fmt::format("{}_{}.{}", opts.inputs.front().stem().string(),
                             mjd, opts.extension);
    }
    require_distinct(opts.inputs, output);
    psrfits::write_fold_archive(output, spec, amps);
    ctx.info("psrio cutout: wrote {} ({} pol, {} chan, {} bin)",
             output.string(), npol, nchan, opts.nbin);
}

} // namespace

void add_cutout(CLI::App& app, Context& ctx) {
    const auto opts = std::make_shared<CutoutOptions>();
    auto* sub       = app.add_subcommand(
        "cutout", "Cut a burst out of search-mode PSRFITS into a PSRCHIVE "
                  "archive, as dspsr -s would");
    sub->add_option("inputs", opts->inputs,
                    "Search-mode PSRFITS files (.sf, .fits) in time order")
        ->required()
        ->expected(1, -1)
        ->check(CLI::ExistingFile);
    sub->add_option("-D,--dm", opts->dm, "Dispersion measure (pc cm^-3)")
        ->required();
    sub->add_flag("-K,--dedisperse", opts->dedisperse,
                  "Remove inter-channel dispersion delays (dspsr -K)");
    sub->add_option("-b,--nbin", opts->nbin, "Phase bins in the profile")
        ->required();
    sub->add_flag("--scloffs", opts->scloffs,
                  "Apply DAT_SCL, DAT_OFFS and ZERO_OFF (dspsr -scloffs)");
    sub->add_option("--toa", opts->toa,
                    "Burst arrival time (MJD) at --toa-freq; the window is "
                    "centred on it");
    sub->add_option("--toa-freq", opts->toa_freq,
                    "Frequency of --toa: top, bottom, centre, inf, or MHz")
        ->capture_default_str();
    sub->add_option("-S,--start", opts->start_seconds,
                    "Window start, seconds from the start of the data "
                    "(band-centre time with -K)");
    sub->add_option("-T,--length", opts->length_seconds,
                    "Window length in seconds, rounded up to whole samples "
                    "per bin (default without -K: the band's sweep, twice, "
                    "plus --margin)");
    sub->add_option("--margin", opts->margin,
                    "Fraction added to the default window length")
        ->capture_default_str()
        ->check(CLI::NonNegativeNumber);
    sub->add_option("-O,--output", opts->output,
                    "Output archive (default: <first input stem>_<MJD>.<ext>)");
    sub->add_option("-e,--ext", opts->extension,
                    "Extension of the default name")
        ->capture_default_str();
    add_threads_option(*sub);
    sub->callback([opts, &ctx] { run_cutout(*opts, ctx); });
}

} // namespace psrio::cli
