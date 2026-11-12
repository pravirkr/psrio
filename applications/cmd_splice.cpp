#include "common.hpp"

#include <CLI/CLI.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fmt/format.h>
#include <memory>
#include <psrio/detail/exceptions.hpp>
#include <span>
#include <utility>
#include <vector>

namespace psrio::cli {

namespace {

struct SpliceOptions {
    std::vector<std::filesystem::path> inputs;
    std::filesystem::path output;
    std::uint64_t gulp{kDefaultGulp};
    bool keep_order{false};
};

struct Band {
    std::unique_ptr<Input> input;
    std::vector<double> freqs;
    std::size_t row_bytes{0};
};

void check_compatible(const Band& ref, const Band& band) {
    const auto& a      = ref.input->header();
    const auto& b      = band.input->header();
    const auto& name   = band.input->path();
    const double tsamp = ref.input->tsamp();
    if (std::abs(band.input->tsamp() - tsamp) > 1e-9 * tsamp) {
        throw ValidationError(fmt::format("psrio: {} has tsamp {} s, not {} s",
                                          name.string(), band.input->tsamp(),
                                          tsamp));
    }
    if (b.nbits != a.nbits ||
        b.samples_are_signed() != a.samples_are_signed()) {
        throw ValidationError(
            fmt::format("psrio: {} has {}-bit samples, not {}-bit",
                        name.string(), b.nbits, a.nbits));
    }
    if (b.nifs != a.nifs) {
        throw ValidationError(fmt::format("psrio: {} has {} IFs, not {}",
                                          name.string(), b.nifs, a.nifs));
    }
    if (a.tstart.has_value() != b.tstart.has_value() ||
        (a.tstart &&
         std::abs(*a.tstart - *b.tstart) * 86400.0 >= 0.5 * tsamp)) {
        throw ValidationError(
            fmt::format("psrio: {} starts at a different time than {}",
                        name.string(), ref.input->path().string()));
    }
}

void run_splice(const SpliceOptions& opts, const Context& ctx) {
    std::vector<Band> bands;
    for (const auto& path : opts.inputs) {
        require_distinct(path, opts.output);
        Band band;
        band.input = std::make_unique<Input>(path);
        (void)band.input->block();
        const auto& header = band.input->header();
        band.freqs         = channel_frequencies(header);
        const auto bits    = static_cast<std::uint64_t>(header.nchans) *
                             static_cast<std::uint64_t>(header.nbits);
        if (bits % 8U != 0U) {
            throw ValidationError(fmt::format(
                "psrio: {} channels of {} bits in {} are not whole bytes",
                header.nchans, header.nbits, path.string()));
        }
        band.row_bytes = static_cast<std::size_t>(bits / 8U);
        if (!bands.empty()) {
            check_compatible(bands.front(), band);
        }
        bands.push_back(std::move(band));
    }
    if (!opts.keep_order) {
        std::ranges::stable_sort(bands, [](const Band& a, const Band& b) {
            return a.freqs.front() > b.freqs.front();
        });
    }

    std::uint64_t nsamples = bands.front().input->nsamples();
    for (const auto& band : bands) {
        if (band.input->nsamples() != nsamples) {
            ctx.warn("{} has {} samples; the output stops at the shortest",
                     band.input->path().string(), band.input->nsamples());
        }
        nsamples = std::min(nsamples, band.input->nsamples());
    }

    std::vector<double> freqs;
    for (const auto& band : bands) {
        freqs.insert(freqs.end(), band.freqs.begin(), band.freqs.end());
    }
    auto header = bands.front().input->header();
    set_channel_frequencies(header, freqs);
    if (!header.frequency_table.empty()) {
        ctx.info("psrio splice: the bands are not contiguous; writing a "
                 "frequency table");
    }
    Output output(opts.output, header);
    auto& writer = output.sigproc_writer();

    const auto nifs     = static_cast<std::size_t>(header.nifs);
    const auto gulp     = std::max<std::uint64_t>(opts.gulp, 1U);
    std::size_t out_row = 0;
    for (const auto& band : bands) {
        out_row += band.row_bytes;
    }
    std::vector<std::vector<std::byte>> blocks(bands.size());
    for (std::size_t i = 0; i < bands.size(); ++i) {
        blocks[i].resize(static_cast<std::size_t>(gulp) * nifs *
                         bands[i].row_bytes);
    }
    std::vector<std::byte> interleaved(static_cast<std::size_t>(gulp) * nifs *
                                       out_row);

    std::uint64_t done = 0;
    while (done < nsamples) {
        const auto want = std::min(gulp, nsamples - done);
        for (std::size_t i = 0; i < bands.size(); ++i) {
            const auto got =
                bands[i].input->block().read_block(want, blocks[i]);
            if (got != want) {
                throw IoError(fmt::format("psrio: short read from {}",
                                          bands[i].input->path().string()));
            }
        }
        const auto nrows = static_cast<std::size_t>(want) * nifs;
        std::byte* dst   = interleaved.data();
        for (std::size_t row = 0; row < nrows; ++row) {
            for (std::size_t i = 0; i < bands.size(); ++i) {
                const auto width = bands[i].row_bytes;
                std::memcpy(dst, blocks[i].data() + (row * width), width);
                dst += width;
            }
        }
        writer.write_block(
            want,
            std::span<const std::byte>(interleaved).first(nrows * out_row));
        done += want;
    }
    output.close();
    ctx.info("psrio splice: wrote {} samples x {} channels from {} files to {}",
             output.nsamples_written(), freqs.size(), bands.size(),
             opts.output.string());
}

} // namespace

void add_splice(CLI::App& app, Context& ctx) {
    const auto opts = std::make_shared<SpliceOptions>();
    auto* sub       = app.add_subcommand(
        "splice",
        "Join sub-band filterbanks across frequency (SIGPROC 'splice')");
    sub->add_option("inputs", opts->inputs,
                    "Filterbanks of one observation, one per sub-band")
        ->required()
        ->expected(2, -1)
        ->check(CLI::ExistingFile);
    sub->add_option("-o,--output", opts->output, "Output filterbank")
        ->required();
    sub->add_option("--gulp", opts->gulp, "Time samples per block")
        ->capture_default_str()
        ->check(CLI::PositiveNumber);
    sub->add_flag("--keep-order", opts->keep_order,
                  "Keep the input order instead of sorting by descending "
                  "frequency");
    sub->callback([opts, &ctx] { run_splice(*opts, ctx); });
}

} // namespace psrio::cli
