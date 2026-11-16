#include "common.hpp"

#include <CLI/CLI.hpp>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <psrio/detail/exceptions.hpp>
#include <psrio/digifil.hpp>
#include <vector>

namespace psrio::cli {

namespace {

struct DigifilOptions {
    std::vector<std::filesystem::path> inputs;
    std::filesystem::path output;
    digifil::Config config;
    std::optional<std::uint64_t> gulp;
    std::optional<double> max_ram_mb;
    std::optional<double> block_size;
    int npol{1};
};

/// dspsr's default output name: the first input's stem with `.fil`, in the
/// working directory.
std::filesystem::path default_output(const std::filesystem::path& input) {
    auto name = input.filename();
    name.replace_extension(".fil");
    return name;
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
    const auto output =
        opts.output.empty() ? default_output(opts.inputs.front()) : opts.output;
    require_distinct(opts.inputs, output);
    auto config = opts.config;
    config.gulp = opts.gulp;
    if (opts.max_ram_mb.has_value()) {
        config.max_ram_mb = *opts.max_ram_mb;
    }
    config.validate();
    // PSRFITS is decoded as dspsr does: scales and offsets only with
    // --scloffs, no weights, and -P picking a stored product.
    Input input(opts.inputs, digifil::decode_options(config));
    if (input.is_psrfits()) {
        if (!config.gulp.has_value()) {
            config.gulp = digifil::psrfits_block_samples(
                input.stored_nbits(), input.stored_npol(), input.nchans(),
                input.subint_samples(), input.sums_polarisations(),
                config.tscrunch, config.fscrunch, config.max_ram_bytes());
        }
        config.poln_select = -1;
    } else if (config.scloffs) {
        ctx.warn("psrio digifil: --scloffs only applies to PSRFITS input");
    }
    const auto summary =
        digifil::run(input.block(), input.header(), output, config);
    ctx.debug("psrio digifil: {} input samples per block", summary.gulp);
    ctx.info("psrio digifil: wrote {} samples at {} bits from {} input samples "
             "to {}",
             summary.nsamples_out, config.nbits, summary.nsamples_in,
             output.string());
}

} // namespace

void add_digifil(CLI::App& app, Context& ctx) {
    const auto opts = std::make_shared<DigifilOptions>();
    auto& config    = opts->config;
    auto* sub       = app.add_subcommand(
        "digifil", "Rescale and redigitise a filterbank, as dspsr's digifil");
    sub->add_option("inputs", opts->inputs,
                    "Input files in time order: .fil or .sf/.fits, or one .h5")
        ->required()
        ->expected(1, -1)
        ->check(CLI::ExistingFile);
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
    sub->add_option("-B,--block-size", opts->block_size,
                    "Accepted for dspsr compatibility; ignored");
    add_threads_option(*sub);
    sub->callback([opts, &ctx] { run_digifil(*opts, ctx); });
}

} // namespace psrio::cli
