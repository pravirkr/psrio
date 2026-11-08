#include "common.hpp"

#include <CLI/CLI.hpp>
#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <dmt/algorithms/ddmt.hpp>
#include <dmt/algorithms/fdmt.hpp>
#include <dmt/algorithms/sdmt.hpp>
#include <dmt/common/backend.hpp>
#include <dmt/common/types.hpp>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <iterator>
#include <memory>
#include <optional>
#include <psrio/detail/exceptions.hpp>
#include <psrio/header.hpp>
#include <psrio/timeseries.hpp>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>
#include <vector>

namespace psrio::cli {

namespace {

struct DedisperseOptions {
    std::vector<std::filesystem::path> inputs;
    std::filesystem::path output;
    std::vector<double> dms;
    std::optional<double> dm_start;
    std::optional<double> dm_end;
    std::optional<double> dm_step;
    std::string algo{"auto"};
    std::string fref{"top"};
    std::string ignorechans;
    std::uint64_t gulp{kDefaultGulp};
    int threads{0};
    std::string format{"auto"};
};

/// Largest DM count for which `auto` picks DDMT over SDMT.
constexpr std::size_t kDirectMaxTrials = 8;

std::vector<float> dm_trials(const DedisperseOptions& opts) {
    std::vector<float> dms;
    const bool grid = opts.dm_start || opts.dm_end || opts.dm_step;
    if (grid == !opts.dms.empty()) {
        throw ValidationError(
            "psrio: give either --dm or --dm-start, --dm-end, and --dm-step");
    }
    if (!grid) {
        std::ranges::transform(
            opts.dms, std::back_inserter(dms),
            [](double dm) { return static_cast<float>(dm); });
    } else {
        if (!opts.dm_start || !opts.dm_end || !opts.dm_step ||
            !(*opts.dm_step > 0.0) || *opts.dm_end < *opts.dm_start) {
            throw ValidationError(
                "psrio: a DM grid needs --dm-start <= --dm-end and a positive "
                "--dm-step");
        }
        const auto ntrials =
            static_cast<std::size_t>(std::floor(
                ((*opts.dm_end - *opts.dm_start) / *opts.dm_step) + 1e-9)) +
            1U;
        for (std::size_t i = 0; i < ntrials; ++i) {
            dms.push_back(static_cast<float>(
                *opts.dm_start + (static_cast<double>(i) * *opts.dm_step)));
        }
    }
    for (const float dm : dms) {
        if (!std::isfinite(dm) || dm < 0.0F) {
            throw ValidationError(fmt::format(
                "psrio: DM {} must be finite and non-negative", dm));
        }
    }
    return dms;
}

std::uint64_t parse_channel(std::string_view text, std::uint64_t nchans) {
    std::uint64_t value  = 0;
    const auto* end      = text.data() + text.size();
    const auto [ptr, ec] = std::from_chars(text.data(), end, value);
    if (ec != std::errc{} || ptr != end) {
        throw ValidationError(
            fmt::format("psrio: '{}' is not a channel index", text));
    }
    if (value >= nchans) {
        throw ValidationError(fmt::format(
            "psrio: channel {} is outside the {} channels", value, nchans));
    }
    return value;
}

/// Parse "0,1,100-120" (or a file of such tokens) into file channel indices.
std::vector<std::uint64_t> parse_channels(std::string spec,
                                          std::uint64_t nchans) {
    if (std::filesystem::is_regular_file(spec)) {
        std::ifstream in(spec);
        spec.assign(std::istreambuf_iterator<char>(in),
                    std::istreambuf_iterator<char>());
    }
    std::ranges::replace_if(
        spec,
        [](char c) { return c == '\n' || c == '\r' || c == ' ' || c == '\t'; },
        ',');
    std::vector<std::uint64_t> channels;
    std::string_view rest = spec;
    while (!rest.empty()) {
        const auto comma = rest.find(',');
        const auto token = rest.substr(0, comma);
        rest = comma == std::string_view::npos ? std::string_view{}
                                               : rest.substr(comma + 1);
        if (token.empty()) {
            continue;
        }
        const auto dash = token.find('-');
        if (dash == std::string_view::npos) {
            channels.push_back(parse_channel(token, nchans));
            continue;
        }
        const auto lo = parse_channel(token.substr(0, dash), nchans);
        const auto hi = parse_channel(token.substr(dash + 1), nchans);
        if (hi < lo) {
            throw ValidationError(
                fmt::format("psrio: channel range '{}' is reversed", token));
        }
        for (auto chan = lo; chan <= hi; ++chan) {
            channels.push_back(chan);
        }
    }
    return channels;
}

/// Transpose @p nrows rows of @p ncols values into column-major @p out.
void transpose(std::span<const float> in,
               std::span<float> out,
               std::size_t nrows,
               std::size_t ncols) {
    constexpr std::size_t kTile = 64;
    for (std::size_t r0 = 0; r0 < nrows; r0 += kTile) {
        const auto r1 = std::min(nrows, r0 + kTile);
        for (std::size_t c0 = 0; c0 < ncols; c0 += kTile) {
            const auto c1 = std::min(ncols, c0 + kTile);
            for (std::size_t r = r0; r < r1; ++r) {
                for (std::size_t c = c0; c < c1; ++c) {
                    out[(c * nrows) + r] = in[(r * ncols) + c];
                }
            }
        }
    }
}

/// Output path for each trial and the time series format.
struct Destinations {
    std::vector<std::filesystem::path> paths;
    TimeSeriesFormat format{TimeSeriesFormat::kSigprocTim};
};

Destinations destinations(const DedisperseOptions& opts,
                          std::span<const float> dms) {
    Destinations out;
    auto prefix      = opts.output;
    const bool named = is_tim_path(prefix) || is_presto_path(prefix);
    if (opts.format == "dat" ||
        (opts.format == "auto" && is_presto_path(prefix))) {
        out.format = TimeSeriesFormat::kPrestoDat;
    }
    const auto* const suffix =
        out.format == TimeSeriesFormat::kPrestoDat ? ".dat" : ".tim";
    if (dms.size() == 1U && named) {
        auto path = prefix;
        path.replace_extension(suffix);
        out.paths.push_back(path);
        return out;
    }
    if (std::filesystem::is_directory(prefix)) {
        prefix /= opts.inputs.front().stem();
    } else if (named) {
        prefix.replace_extension();
    }
    for (const float dm : dms) {
        out.paths.emplace_back(
            fmt::format("{}_DM{:.2f}{}", prefix.string(), dm, suffix));
    }
    return out;
}

/// Time series metadata for one trial.
Header series_header(const Input& input, double dm, double offset_samples) {
    const auto& raw = input.header();
    auto header     = Header::from_sigproc(raw, input.path().string());
    const double df = std::abs(*raw.foff);
    const double top =
        std::max(*raw.fch1, *raw.fch1 + (*raw.foff * (raw.nchans - 1)));
    header.data_type = "time series";
    header.nbits     = 32;
    header.nifs      = 1;
    header.fch1      = top;
    header.foff      = -df;
    header.dm        = dm;
    header.nsamples  = 0;
    header.tstart += offset_samples * input.tsamp() / 86400.0;
    return header;
}

class Writers {
public:
    Writers(const Input& input,
            const Destinations& dest,
            std::span<const float> dms,
            double offset_samples) {
        for (std::size_t i = 0; i < dms.size(); ++i) {
            require_distinct(input.path(), dest.paths[i]);
            m_writers.emplace_back(dest.paths[i],
                                   series_header(input,
                                                 static_cast<double>(dms[i]),
                                                 offset_samples),
                                   dest.format);
        }
    }

    TimeSeriesWriter& operator[](std::size_t i) { return m_writers[i]; }

    void close() {
        for (auto& writer : m_writers) {
            writer.close();
        }
    }

    [[nodiscard]] std::uint64_t nsamples() const {
        return m_writers.empty() ? 0U : m_writers.front().nsamples_written();
    }

private:
    std::vector<TimeSeriesWriter> m_writers;
};

/// Reads time-major blocks from @p input and hands channel-major ones, with
/// channels ascending, to @p sink.
template <typename Sink>
void stream_channel_major(Input& input, std::uint64_t gulp, const Sink& sink) {
    const auto nchans = static_cast<std::size_t>(input.nchans());
    std::vector<float> time_major(static_cast<std::size_t>(gulp) * nchans);
    std::vector<float> chan_major(time_major.size());
    while (const auto got = input.read_floats(gulp, time_major)) {
        const auto n = static_cast<std::size_t>(got);
        transpose(std::span<const float>(time_major).first(n * nchans),
                  std::span<float>(chan_major).first(n * nchans), n, nchans);
        sink(std::span<const float>(chan_major).first(n * nchans), got);
    }
}

template <typename Engine>
void run_direct(Engine& engine,
                Input& input,
                Writers& writers,
                std::size_t ndm,
                std::uint64_t gulp) {
    std::vector<float> dmt;
    stream_channel_major(
        input, gulp,
        [&](std::span<const float> waterfall, std::uint64_t nsamp) {
            const auto nout = engine.get_output_nsamps(nsamp);
            dmt.resize(ndm * nout);
            engine.execute(waterfall, dmt);
            for (std::size_t d = 0; d < ndm; ++d) {
                writers[d].write_samples(
                    std::span<const float>(dmt).subspan(d * nout, nout));
            }
        });
}

/// Band edges and channel mask in dmt's ascending channel order.
struct Band {
    float f_min;
    float f_max;
    float tsamp;
    std::uint64_t nchans;
    std::span<const std::uint8_t> kill_mask;
};

void run_fdmt(const DedisperseOptions& opts,
              const Context& ctx,
              Input& input,
              const Band& band,
              std::span<const float> requested,
              const dmt::Exec& exec) {
    if (opts.fref == "center") {
        throw ValidationError("psrio: fdmt supports --fref top or bottom only");
    }
    const auto nchans = static_cast<std::size_t>(band.nchans);
    const auto gulp   = opts.gulp;
    dmt::algorithms::FDMT engine(
        band.f_min, band.f_max, static_cast<dmt::SizeType>(nchans),
        static_cast<dmt::SizeType>(gulp), band.tsamp,
        std::vector<float>(requested.begin(), requested.end()),
        /*use_box_smearing=*/false, "valid", exec);
    const auto& plan = engine.get_plan();
    const auto dms   = plan.get_dm_grid_final();
    if (dms.size() != requested.size()) {
        ctx.warn("fdmt merged {} DM trials into {} delay trials",
                 requested.size(), dms.size());
    }
    // Output t of trial d reads channel c at t + trace[c] (<= 0), so the
    // top-of-band sample of t is t - max_c(-trace[c]).
    std::vector<std::uint64_t> lookback(dms.size(), 0U);
    for (std::size_t d = 0; d < dms.size(); ++d) {
        for (const auto shift : plan.trace_dm(d)) {
            lookback[d] =
                std::max(lookback[d], static_cast<std::uint64_t>(
                                          std::max<dmt::IndexType>(-shift, 0)));
        }
    }
    const auto dt_max   = std::ranges::max(lookback);
    const auto nsamples = input.nsamples();
    if (nsamples <= dt_max) {
        throw ValidationError(fmt::format(
            "psrio: {} samples are too few for a {}-sample dispersion sweep",
            nsamples, dt_max));
    }
    const bool top = opts.fref == "top";
    Writers writers(input, destinations(opts, dms), dms,
                    top ? 0.0 : static_cast<double>(dt_max));
    const auto length = nsamples - dt_max;
    const auto ndm    = dms.size();
    std::vector<float> padded(nchans * static_cast<std::size_t>(gulp));
    std::vector<float> dmt_out(ndm * static_cast<std::size_t>(gulp));
    std::uint64_t base = 0;
    const auto sink    = [&](std::span<const float> waterfall,
                             std::uint64_t nsamp) {
        std::ranges::fill(padded, 0.0F);
        for (std::size_t c = 0; c < nchans; ++c) {
            if (band.kill_mask[c] == 0U) {
                continue;
            }
            std::ranges::copy(
                waterfall.subspan(c * nsamp, static_cast<std::size_t>(nsamp)),
                padded.begin() + static_cast<std::ptrdiff_t>(c * gulp));
        }
        engine.execute(padded, dmt_out);
        for (std::size_t d = 0; d < ndm; ++d) {
            const auto skip = top ? lookback[d] : dt_max;
            const auto lo   = std::max(base, skip);
            const auto hi   = std::min(base + gulp, skip + length);
            if (lo < hi) {
                writers[d].write_samples(
                    std::span<const float>(dmt_out).subspan(
                        (d * gulp) + (lo - base), hi - lo));
            }
        }
        base += gulp;
    };
    stream_channel_major(input, gulp, sink);
    writers.close();
    ctx.info("psrio dedisperse: fdmt wrote {} DM trials of {} samples", ndm,
             writers.nsamples());
}

void run_dedisperse(const DedisperseOptions& opts, const Context& ctx) {
    Input input(opts.inputs);
    auto& source       = input.block();
    const auto& header = input.header();
    if (input.nifs() != 1U) {
        throw ValidationError(fmt::format(
            "psrio: dedisperse needs one IF, not {}; select one with "
            "'psrio digifil -P'",
            input.nifs()));
    }
    if (!header.frequency_table.empty() || !header.fch1 || !header.foff ||
        input.nchans() < 2U) {
        throw ValidationError(
            "psrio: dedisperse needs at least two uniformly spaced channels "
            "(fch1 and foff)");
    }

    const auto dms     = dm_trials(opts);
    const auto nchans  = input.nchans();
    const double foff  = *header.foff;
    const double df    = std::abs(foff);
    const double f_ch0 = *header.fch1;
    const double f_chn = f_ch0 + (foff * static_cast<double>(nchans - 1U));
    const auto f_min = static_cast<float>(std::min(f_ch0, f_chn) - (0.5 * df));
    const auto f_max = static_cast<float>(std::max(f_ch0, f_chn) + (0.5 * df));
    const auto tsamp = static_cast<float>(input.tsamp());
    if (foff < 0.0) {
        source.set_fswap(true);
    }

    std::vector<std::uint8_t> kill_mask(static_cast<std::size_t>(nchans), 1U);
    if (!opts.ignorechans.empty()) {
        for (const auto chan : parse_channels(opts.ignorechans, nchans)) {
            const auto ascending = foff < 0.0 ? nchans - 1U - chan : chan;
            kill_mask[static_cast<std::size_t>(ascending)] = 0U;
        }
    }

    std::optional<float> f_ref;
    if (opts.fref == "top") {
        f_ref = f_max;
    } else if (opts.fref == "bottom") {
        f_ref = f_min;
    } else if (opts.fref == "center") {
        f_ref = 0.5F * (f_min + f_max);
    }

    const int threads = opts.threads > 0
                            ? opts.threads
                            : static_cast<int>(std::max(
                                  1U, std::thread::hardware_concurrency()));
    const auto exec   = dmt::Exec::cpu(threads);
    auto algo         = opts.algo;
    if (algo == "auto") {
        algo = dms.size() <= kDirectMaxTrials ? "ddmt" : "sdmt";
    }
    const auto nsamples = input.nsamples();

    if (algo == "fdmt") {
        run_fdmt(opts, ctx, input,
                 Band{
                     .f_min     = f_min,
                     .f_max     = f_max,
                     .tsamp     = tsamp,
                     .nchans    = nchans,
                     .kill_mask = kill_mask,
                 },
                 dms, exec);
        return;
    }

    const auto dest       = destinations(opts, dms);
    const auto run_engine = [&](auto& engine) {
        const auto& plan = engine.get_plan();
        const auto nout  = engine.get_output_nsamps(nsamples);
        if (nout == 0U) {
            throw ValidationError(fmt::format(
                "psrio: {} samples are shorter than the dispersion sweep",
                nsamples));
        }
        const auto max_delay = nsamples - nout;
        const auto gulp = std::max<std::uint64_t>(opts.gulp, 2U * max_delay);
        if (gulp != opts.gulp) {
            ctx.debug("psrio dedisperse: gulp raised to {} samples to cover "
                      "the {}-sample sweep",
                      gulp, max_delay);
        }
        Writers writers(input, dest, dms,
                        static_cast<double>(plan.get_time_offset()));
        run_direct(engine, input, writers, dms.size(), gulp);
        writers.close();
        ctx.info("psrio dedisperse: {} wrote {} DM trials of {} samples", algo,
                 dms.size(), writers.nsamples());
    };
    if (algo == "ddmt") {
        dmt::algorithms::DDMT engine(f_min, f_max,
                                     static_cast<dmt::SizeType>(nchans), tsamp,
                                     dms, exec, 32, kill_mask, 1, f_ref);
        run_engine(engine);
    } else {
        dmt::algorithms::SDMT engine(f_min, f_max,
                                     static_cast<dmt::SizeType>(nchans), tsamp,
                                     dms, exec, 32, kill_mask, 1, f_ref);
        run_engine(engine);
    }
}

} // namespace

void add_dedisperse(CLI::App& app, Context& ctx) {
    const auto opts = std::make_shared<DedisperseOptions>();
    auto* sub       = app.add_subcommand(
        "dedisperse",
        "Dedisperse a filterbank into time series with dmt (SIGPROC "
        "'dedisperse', PRESTO 'prepdata')");
    sub->add_option("inputs", opts->inputs,
                    "Input files in time order: .fil or .sf/.fits, or one .h5")
        ->required()
        ->expected(1, -1)
        ->check(CLI::ExistingFile);
    sub->add_option("-o,--output", opts->output,
                    "Output .tim/.dat for one DM; otherwise a prefix or "
                    "directory for <prefix>_DM<dm>.tim")
        ->required();
    auto* dm = sub->add_option("-d,--dm", opts->dms,
                               "DM trials in pc/cm^3, comma separated")
                   ->delimiter(',');
    auto* start =
        sub->add_option("--dm-start", opts->dm_start, "First DM of a grid");
    auto* end = sub->add_option("--dm-end", opts->dm_end, "Last DM of a grid");
    auto* step =
        sub->add_option("--dm-step", opts->dm_step, "Spacing of a DM grid");
    dm->excludes(start)->excludes(end)->excludes(step);
    sub->add_option("--algo", opts->algo,
                    "dmt engine: ddmt, sdmt, fdmt, or auto (ddmt for up to 8 "
                    "DMs, else sdmt)")
        ->capture_default_str()
        ->check(CLI::IsMember({"auto", "ddmt", "sdmt", "fdmt"}));
    sub->add_option("--fref", opts->fref,
                    "Frequency whose arrival time sample 0 marks: top, "
                    "bottom, or center of the band")
        ->capture_default_str()
        ->check(CLI::IsMember({"top", "bottom", "center"}));
    sub->add_option("--ignorechans,--zap-chans", opts->ignorechans,
                    "Channels to drop, in file order: '0,1,100-120' or a file "
                    "of them");
    sub->add_option("--gulp", opts->gulp,
                    "Time samples per block (raised to twice the sweep)")
        ->capture_default_str()
        ->check(CLI::PositiveNumber);
    sub->add_option("-j,--threads", opts->threads,
                    "Threads (default: all cores)")
        ->check(CLI::NonNegativeNumber);
    sub->add_option("--format", opts->format,
                    "Output format: tim, dat, or auto (from -o, else tim)")
        ->capture_default_str()
        ->check(CLI::IsMember({"auto", "tim", "dat"}));
    sub->callback([opts, &ctx] { run_dedisperse(*opts, ctx); });
}

} // namespace psrio::cli
