#include "common.hpp"
#include "progress.hpp"

#include <CLI/CLI.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <psrio/detail/exceptions.hpp>
#include <psrio/detail/readahead.hpp>
#include <psrio/digifil.hpp>
#include <string>
#include <system_error>
#include <vector>

namespace psrio::cli {

namespace {

struct DigifilOptions {
    std::vector<std::string> inputs;
    std::optional<std::filesystem::path> file_list;
    std::filesystem::path output;
    digifil::Config config;
    std::optional<std::uint64_t> gulp;
    std::optional<double> max_ram_mb;
    std::optional<double> block_size;
    double io_buffer_mb{128.0};
    std::string progress{"auto"};
    bool no_atomic{false};
    int npol{1};
};

/// dspsr's default output name: the first input's stem with `.fil`, in the
/// working directory.
std::filesystem::path default_output(const std::filesystem::path& input) {
    auto name = input.filename();
    name.replace_extension(".fil");
    return name;
}

std::filesystem::path staging_path_for(const std::filesystem::path& output) {
    auto staged = output;
    staged += ".psrio.part";
    return staged;
}

void run_digifil(const DigifilOptions& opts, const Context& ctx) {
    if (opts.npol != 1) {
        throw ValidationError(fmt::format(
            "psrio: digifil -d {} is not supported; output is Stokes I or "
            "one polarisation (-P)",
            opts.npol));
    }
    if (opts.gulp.has_value() && opts.max_ram_mb.has_value()) {
        throw ValidationError("psrio: give either --gulp or -U, not both");
    }
    if (!std::isfinite(opts.io_buffer_mb) || opts.io_buffer_mb <= 0.0) {
        throw ValidationError(
            fmt::format("psrio: --io-buffer-mb must be positive, got {}",
                        opts.io_buffer_mb));
    }
    const auto progress_mode = parse_progress_mode(opts.progress);
    const auto input_paths   = expand_input_paths(opts.inputs, opts.file_list);

    auto config = opts.config;
    config.gulp = opts.gulp;
    if (opts.max_ram_mb.has_value()) {
        config.max_ram_mb = *opts.max_ram_mb;
    }
    config.validate();

    const auto block_mb =
        static_cast<double>(::psrio::detail::ReadAhead::kBlockBytes) /
        (1024.0 * 1024.0);
    const auto read_ahead_blocks = std::max<std::size_t>(
        1U,
        static_cast<std::size_t>(std::llround(opts.io_buffer_mb / block_mb)));

    config.drop_cache = ctx.drop_cache;

    ProgressReporter reporter(ctx, progress_mode, "psrio digifil");
    ctx.debug("psrio digifil: opening {} input file{}...", input_paths.size(),
              input_paths.size() == 1U ? "" : "s");
    // PSRFITS is decoded as dspsr does: scales and offsets only with
    // --scloffs, no weights, and -P picking a stored product.
    Input input(
        input_paths, digifil::decode_options(config),
        [&reporter](std::size_t done, std::size_t total) {
            reporter.update_scan(done, total);
        },
        read_ahead_blocks, ctx.drop_cache, !ctx.force_contiguity);
    reporter.finish_scan();

    const auto output =
        opts.output.empty() ? default_output(input.path()) : opts.output;
    require_distinct(input.paths(), output);

    if (input.is_psrfits()) {
        if (!config.gulp.has_value()) {
            config.gulp = digifil::psrfits_block_samples(
                input.stored_nbits(), input.stored_npol(), input.nchans(),
                input.subint_samples(), input.sums_polarisations(),
                config.tscrunch, config.fscrunch, config.max_ram_bytes());
        }
        config.poln_select = -1;
    } else {
        if (config.scloffs) {
            ctx.warn("psrio digifil: --scloffs only applies to PSRFITS input");
        }
        if (opts.npol == 1 && config.poln_select < 0 && input.nifs() > 1U) {
            config.pscrunch = true;
        }
    }

    const auto nsamples_in  = input.nsamples();
    const auto tsamp_in     = input.tsamp();
    const auto ffac         = std::max<std::uint64_t>(config.fscrunch, 1U);
    const auto tfac         = std::max<std::uint64_t>(config.tscrunch, 1U);
    const auto nchans_out   = input.nchans() / ffac;
    const auto nsamples_est = nsamples_in / tfac;
    const double obs_sec    = static_cast<double>(nsamples_in) * tsamp_in;

    ctx.debug("psrio digifil: input {} ({} file{}, {}, {} samples, {} obs, "
              "{} ch, tsamp {} us)",
              input.is_psrfits() ? "PSRFITS" : "SIGPROC", input.nfiles(),
              input.nfiles() == 1U ? "" : "s",
              format_bytes(input.total_input_bytes()), nsamples_in,
              format_duration(obs_sec), input.nchans(), tsamp_in * 1.0e6);
    if (input.nfiles() > 1U) {
        std::string_view order_note;
        if (ctx.force_contiguity) {
            order_note = " (forced contiguous)";
        } else if (input.reordered()) {
            order_note = " (sorted chronologically)";
        }
        ctx.debug("psrio digifil: range {} .. {}{}",
                  input.paths().front().filename().string(),
                  input.paths().back().filename().string(), order_note);
    }
    ctx.debug("psrio digifil: output {} ({} ch, {} bits, tsamp {} us, "
              "~{} samples)",
              output.string(), nchans_out, config.nbits,
              tsamp_in * static_cast<double>(tfac) * 1.0e6, nsamples_est);
    ctx.debug("psrio digifil: rescale={}s (constant={}, clip={}), "
              "scloffs={}, threads={}, io_buffer={} MiB, drop_cache={}, "
              "cont={}",
              config.rescale_seconds, config.rescale_constant,
              config.rescale_clip, config.scloffs, default_threads(),
              read_ahead_blocks *
                  (::psrio::detail::ReadAhead::kBlockBytes >> 20U),
              config.drop_cache, ctx.force_contiguity);

    const bool use_staging = !opts.no_atomic;
    const auto write_path  = use_staging ? staging_path_for(output) : output;
    if (use_staging) {
        require_distinct(input.paths(), write_path);
    }

    reporter.begin_stream(nsamples_in, tsamp_in, input.total_input_bytes(),
                          input.nfiles());
    std::size_t last_file_logged = 0;
    digifil::Summary summary{};
    try {
        summary = digifil::run(
            input.block(), input.header(), write_path, config, 0, std::nullopt,
            [&](const digifil::Progress& prog) {
                const auto file_idx = input.file_index_at(
                    prog.samples_in_done > 0U ? prog.samples_in_done - 1U : 0U);
                if (ctx.verbosity >= 2 && !ctx.quiet &&
                    file_idx != last_file_logged &&
                    file_idx < input.paths().size()) {
                    reporter.clear_line();
                    ctx.err << fmt::format(
                        "psrio digifil: file {}/{}: {}\n", file_idx + 1U,
                        input.nfiles(),
                        input.paths()[file_idx].filename().string());
                    last_file_logged = file_idx;
                }
                reporter.update_stream(prog.samples_in_done, file_idx);
            });
        if (use_staging) {
            std::filesystem::rename(write_path, output);
        }
    } catch (...) {
        reporter.clear_line();
        if (use_staging) {
            std::error_code rm_ec;
            std::filesystem::remove(write_path, rm_ec);
        }
        throw;
    }

    const double elapsed = reporter.finish();
    const double avg_rate =
        elapsed > 1.0e-3
            ? static_cast<double>(input.total_input_bytes()) / elapsed
            : 0.0;
    ctx.debug("psrio digifil: {} input samples per block", summary.gulp);
    ctx.info("psrio digifil: wrote {} samples at {} bits from {} input samples "
             "to {} in {} ({})",
             summary.nsamples_out, config.nbits, summary.nsamples_in,
             output.string(), format_duration(elapsed), format_rate(avg_rate));
}

} // namespace

void add_digifil(CLI::App& app, Context& ctx) {
    const auto opts = std::make_shared<DigifilOptions>();
    auto& config    = opts->config;
    auto* sub       = app.add_subcommand(
        "digifil", "Rescale and redigitise a filterbank, as dspsr's digifil");
    sub->add_option("inputs", opts->inputs,
                    "Input files, glob patterns (e.g. \"/dir/*.sf\"), or "
                    "@filelist.txt; multiple files are sorted by time")
        ->expected(0, -1);
    sub->add_option("--file-list", opts->file_list,
                    "Text file listing input paths (one per line, # comments)");
    sub->add_option("-o,--output", opts->output,
                    "Output filterbank (default: first input's stem with .fil, "
                    "in the working directory)");
    sub->add_option("-b,--nbits", config.nbits,
                    "Output bits: 1, 2, 4, 8, 16, or -32 for float")
        ->capture_default_str()
        ->check(CLI::IsMember({1, 2, 4, 8, 16, -32}));
    sub->add_option("-t,--tscrunch", config.tscrunch,
                    "Time samples summed into one")
        ->capture_default_str();
    sub->add_option("-f,--fscrunch", config.fscrunch,
                    "Channels summed into one; must divide nchans")
        ->capture_default_str();
    sub->add_option("-I,--rescale-seconds", config.rescale_seconds,
                    "Rescale interval in seconds; 0 disables rescaling")
        ->capture_default_str();
    sub->add_flag("-c,--rescale-constant", config.rescale_constant,
                  "Keep the offset and scale from the first interval");
    sub->add_option("--rescale-clip", config.rescale_clip,
                    "Robust rescaling: leave values more than this many "
                    "sigma from a channel's mean out of its next statistics "
                    "(0: dspsr's plain mean and variance)")
        ->capture_default_str()
        ->check(CLI::NonNegativeNumber);
    sub->add_option("-s,--scale", config.scale_fac,
                    "Extra scale applied before digitising")
        ->capture_default_str();
    sub->add_option("-P,--poln-select", config.poln_select,
                    "Keep only this IF (-1 keeps all); for PSRFITS, this "
                    "stored product (0 AA, 1 BB, ...)")
        ->capture_default_str();
    sub->add_flag("--scloffs", config.scloffs,
                  "PSRFITS: apply DAT_SCL, DAT_OFFS and ZERO_OFF");
    sub->add_option("-d,--npol", opts->npol,
                    "Output polarisations; only 1 (Stokes I or -P) is "
                    "supported")
        ->capture_default_str();
    sub->add_option("--gulp", opts->gulp,
                    "Input samples per block (default: dspsr's block for "
                    "the RAM limit)")
        ->check(CLI::PositiveNumber);
    sub->add_option("-U,--max-ram", opts->max_ram_mb,
                    "RAM limit in MB that sets the input block, as dspsr "
                    "(default 256)")
        ->check(CLI::PositiveNumber);
    sub->add_option("--io-buffer-mb", opts->io_buffer_mb,
                    "Read-ahead I/O buffer per file in MiB (independent of "
                    "-U / rescaling blocks)")
        ->capture_default_str()
        ->check(CLI::PositiveNumber);
    sub->add_option("--progress", opts->progress,
                    "Progress display: auto (bar on TTY, periodic log on "
                    "non-TTY), bar, log, or none")
        ->capture_default_str()
        ->check(CLI::IsMember({"auto", "bar", "log", "none"}));
    sub->add_flag("--no-atomic", opts->no_atomic,
                  "Write directly to the output path instead of staging via "
                  "<output>.psrio.part");
    sub->add_option("-B,--block-size", opts->block_size,
                    "Accepted for dspsr compatibility; ignored");
    add_threads_option(*sub);
    sub->callback([opts, &ctx] { run_digifil(*opts, ctx); });
}

} // namespace psrio::cli
