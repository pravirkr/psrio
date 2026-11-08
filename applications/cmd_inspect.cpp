#include "common.hpp"

#include <CLI/CLI.hpp>
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <limits>
#include <memory>
#include <optional>
#include <psrio/baseband.hpp>
#include <psrio/detail/exceptions.hpp>
#include <psrio/formats/dada.hpp>
#include <psrio/formats/guppi.hpp>
#include <psrio/formats/sigproc/header.hpp>
#include <psrio/header.hpp>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace psrio::cli {

namespace {

struct InspectOptions {
    std::vector<std::filesystem::path> inputs;
    std::string key;
    bool stats{false};
    std::uint64_t stats_samples{4096};
    bool json{false};
};

struct Field {
    std::string key;
    std::string label;
    std::string text;
    bool numeric{false};
};

class Fields {
public:
    void add(std::string key, std::string label, std::string_view value) {
        m_fields.push_back({
            .key     = std::move(key),
            .label   = std::move(label),
            .text    = std::string(value),
            .numeric = false,
        });
    }

    template <typename T>
        requires std::integral<T> || std::floating_point<T>
    void add(std::string key, std::string label, T value) {
        std::string text;
        if constexpr (std::floating_point<T>) {
            text = std::isfinite(value) ? fmt::format("{}", value) : "null";
        } else {
            text = fmt::format("{}", value);
        }
        m_fields.push_back({
            .key     = std::move(key),
            .label   = std::move(label),
            .text    = std::move(text),
            .numeric = true,
        });
    }

    template <typename T>
    void
    add(std::string key, std::string label, const std::optional<T>& value) {
        if (value.has_value()) {
            add(std::move(key), std::move(label), *value);
        }
    }

    [[nodiscard]] const std::vector<Field>& items() const noexcept {
        return m_fields;
    }

private:
    std::vector<Field> m_fields;
};

constexpr std::array<std::pair<std::string_view, std::string_view>, 9>
    kKeyAliases{
        {
            {"source", "source_name"},
            {"backend", "machine"},
            {"refdm", "dm"},
            {"nsamp", "nsamples"},
            {"fmid", "fcenter"},
            {"src_raj", "ra"},
            {"src_dej", "dec"},
            {"nchan", "nchans"},
            {"nbit", "nbits"},
        },
};

std::string json_string(std::string_view text) {
    std::string out = "\"";
    for (const char c : text) {
        switch (c) {
        case '"':
            out += "\\\"";
            break;
        case '\\':
            out += "\\\\";
            break;
        case '\n':
            out += "\\n";
            break;
        case '\t':
            out += "\\t";
            break;
        default:
            if (static_cast<unsigned char>(c) < 0x20U) {
                out += fmt::format("\\u{:04x}", static_cast<unsigned>(c));
            } else {
                out += c;
            }
        }
    }
    out += '"';
    return out;
}

void add_frequencies(Fields& fields, const sigproc::FilterbankHeader& h) {
    fields.add("fch1", "Frequency of channel 1 (MHz)", h.fch1);
    fields.add("foff", "Channel bandwidth (MHz)", h.foff);
    if (!h.frequency_table.empty()) {
        fields.add("nfreq_table", "Frequency table entries",
                   h.frequency_table.size());
    }
    if ((h.fch1 && h.foff) || !h.frequency_table.empty()) {
        fields.add("ftop", "Top of band (MHz)", h.ftop());
        fields.add("fbottom", "Bottom of band (MHz)", h.fbottom());
        fields.add("fcenter", "Centre frequency (MHz)", h.fcenter());
        fields.add("bandwidth", "Total bandwidth (MHz)", h.bandwidth());
    }
}

void sigproc_fields(Fields& fields,
                    const sigproc::FilterbankHeader& h,
                    std::uint64_t nsamples,
                    std::string_view format) {
    fields.add("format", "Format", format);
    if (h.header_bytes > 0U) {
        fields.add("header_bytes", "Header size (bytes)", h.header_bytes);
        fields.add("data_bytes", "Data size (bytes)", h.data_bytes);
    }
    fields.add("telescope", "Telescope", h.telescope_name());
    fields.add("machine", "Datataking machine", h.machine_name());
    fields.add("data_type", "Data type", h.data_type_name());
    if (h.source_name) {
        fields.add("source_name", "Source name", *h.source_name);
    }
    if (h.rawdatafile) {
        fields.add("rawdatafile", "Raw data file", *h.rawdatafile);
    }
    if (h.src_raj) {
        fields.add("ra", "Source RA (J2000)", h.ra_string());
    }
    if (h.src_dej) {
        fields.add("dec", "Source DEC (J2000)", h.dec_string());
    }
    fields.add("az_start", "Start azimuth (deg)", h.az_start);
    fields.add("za_start", "Start zenith angle (deg)", h.za_start);
    if (h.tstart) {
        fields.add("tstart", "Start time (MJD)", h.tstart);
        fields.add("date", "Start date (UTC)", h.gregorian_date());
    }
    fields.add("tsamp", "Sample time (s)", h.tsamp);
    fields.add("nsamples", "Number of samples", nsamples);
    if (h.tsamp) {
        fields.add("tobs", "Observation length (s)",
                   static_cast<double>(nsamples) * *h.tsamp);
    }
    fields.add("nbits", "Number of bits per sample", h.nbits);
    if (h.nbits == 8) {
        fields.add("signed", "Signed samples", h.samples_are_signed() ? 1 : 0);
    }
    fields.add("nifs", "Number of IFs", h.nifs);
    fields.add("nchans", "Number of channels", h.nchans);
    add_frequencies(fields, h);
    fields.add("dm", "Reference DM (pc/cm^3)", h.refdm);
    fields.add("period", "Folding period (s)", h.period);
    fields.add("ibeam", "Beam number", h.ibeam);
    fields.add("nbeams", "Number of beams", h.nbeams);
    fields.add("barycentric", "Barycentric", h.barycentric);
    fields.add("pulsarcentric", "Pulsarcentric", h.pulsarcentric);
    if (!h.extra.empty()) {
        std::string names;
        for (const auto& extra : h.extra) {
            names += names.empty() ? "" : ",";
            names += extra.name;
        }
        fields.add("extra_keys", "Unrecognised keys", names);
    }
}

void presto_fields(Fields& fields, const Header& h) {
    fields.add("format", "Format", "PRESTO");
    fields.add("telescope", "Telescope", h.telescope);
    fields.add("machine", "Datataking machine", h.backend);
    fields.add("data_type", "Data type", h.data_type);
    fields.add("source_name", "Source name", h.source);
    fields.add("ra", "Source RA (J2000)", h.ra);
    fields.add("dec", "Source DEC (J2000)", h.dec);
    fields.add("tstart", "Start time (MJD)", h.tstart);
    fields.add("date", "Start date (UTC)", h.gregorian_date());
    fields.add("tsamp", "Sample time (s)", h.tsamp);
    fields.add("nsamples", "Number of samples", h.nsamples);
    fields.add("tobs", "Observation length (s)", h.tobs());
    fields.add("nbits", "Number of bits per sample", 32);
    fields.add("nchans", "Number of channels", h.nchans);
    fields.add("fch1", "Frequency of channel 1 (MHz)", h.fch1);
    fields.add("foff", "Channel bandwidth (MHz)", h.foff);
    fields.add("bandwidth", "Total bandwidth (MHz)", h.bandwidth());
    fields.add("dm", "Reference DM (pc/cm^3)", h.dm);
    fields.add("barycentric", "Barycentric", h.barycentric);
}

void baseband_fields(Fields& fields,
                     const BasebandHeader& h,
                     std::string_view format) {
    fields.add("format", "Format", format);
    fields.add("telescope", "Telescope", h.telescope);
    if (!h.backend.empty()) {
        fields.add("machine", "Datataking machine", h.backend);
    }
    fields.add("source_name", "Source name", h.source);
    fields.add("tstart", "Start time (MJD)", h.tstart);
    fields.add("tsamp", "Sample time (s)", h.tsamp);
    if (h.has_nsamples) {
        fields.add("nsamples", "Number of samples", h.nsamples);
        fields.add("tobs", "Observation length (s)",
                   static_cast<double>(h.nsamples) * h.tsamp);
    }
    fields.add("nbits", "Number of bits per sample", h.nbit);
    fields.add("signed", "Signed samples", h.samples_signed ? 1 : 0);
    fields.add("ndim", "Components per value", h.ndim);
    fields.add("npol", "Number of polarisations", h.npol);
    fields.add("nchans", "Number of channels", h.nchan);
    fields.add("nants", "Number of antennas", h.nants);
    fields.add("fch1", "Frequency of channel 1 (MHz)", h.fch1);
    fields.add("foff", "Channel bandwidth (MHz)", h.foff);
    fields.add("fcenter", "Centre frequency (MHz)", h.center_frequency);
    fields.add("bandwidth", "Total bandwidth (MHz)", h.bandwidth);
}

void stats_fields(Fields& fields, Input& input, std::uint64_t max_samples) {
    const auto nsamp = std::min(max_samples, input.nsamples());
    std::vector<float> data(static_cast<std::size_t>(nsamp * input.row()));
    const auto got = input.read_floats(nsamp, data);
    data.resize(static_cast<std::size_t>(got * input.row()));

    double sum           = 0.0;
    double sumsq         = 0.0;
    double lo            = std::numeric_limits<double>::infinity();
    double hi            = -std::numeric_limits<double>::infinity();
    std::uint64_t zeros  = 0;
    std::uint64_t nans   = 0;
    std::uint64_t finite = 0;
    for (const float value : data) {
        if (std::isnan(value)) {
            ++nans;
            continue;
        }
        const auto v = static_cast<double>(value);
        zeros += value == 0.0F ? 1U : 0U;
        if (std::isfinite(v)) {
            ++finite;
            sum += v;
            sumsq += v * v;
        }
        lo = std::min(lo, v);
        hi = std::max(hi, v);
    }
    const auto total = static_cast<double>(data.size());
    const auto n     = static_cast<double>(finite);
    const double mean =
        finite > 0U ? sum / n : std::numeric_limits<double>::quiet_NaN();
    const double rms = finite > 0U ? std::sqrt(sumsq / n)
                                   : std::numeric_limits<double>::quiet_NaN();
    const double var = finite > 0U ? std::max(0.0, (sumsq / n) - (mean * mean))
                                   : std::numeric_limits<double>::quiet_NaN();
    fields.add("stats_samples", "Statistics over samples", got);
    fields.add("min", "Minimum",
               data.empty() || nans == data.size()
                   ? std::numeric_limits<double>::quiet_NaN()
                   : lo);
    fields.add("max", "Maximum",
               data.empty() || nans == data.size()
                   ? std::numeric_limits<double>::quiet_NaN()
                   : hi);
    fields.add("mean", "Mean", mean);
    fields.add("std", "Standard deviation", std::sqrt(var));
    fields.add("rms", "RMS", rms);
    fields.add("zero_fraction", "Fraction of zeros",
               total > 0.0 ? static_cast<double>(zeros) / total : 0.0);
    fields.add("nan_fraction", "Fraction of NaNs",
               total > 0.0 ? static_cast<double>(nans) / total : 0.0);
}

std::string lower(std::string_view text) {
    std::string out(text);
    std::ranges::transform(out, out.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return out;
}

Fields collect(const InspectOptions& opts) {
    Fields fields;
    if (opts.inputs.empty()) {
        throw ValidationError("psrio: an input file is required");
    }
    const auto& path = opts.inputs.front();
    fields.add("filename", "Data file", path.string());
    if (opts.inputs.size() > 1U) {
        fields.add("nfiles", "Files", opts.inputs.size());
        Input input(opts.inputs);
        std::string_view format = "SIGPROC";
        if (input.is_psrfits()) {
            format = "PSRFITS";
        }
        sigproc_fields(fields, input.header(), input.nsamples(), format);
        if (input.is_psrfits()) {
            fields.add("stored_nbits", "Stored bits per sample",
                       input.stored_nbits());
        }
        if (opts.stats) {
            stats_fields(fields, input, opts.stats_samples);
        }
        return fields;
    }
    const auto ext = lower(path.extension().string());

    if (ext == ".dada" || ext == ".raw" || ext == ".guppi") {
        if (opts.stats) {
            throw ValidationError(
                "psrio: --stats is not available for baseband data");
        }
        if (ext == ".dada") {
            const formats::dada::DadaReader reader(path);
            baseband_fields(fields, reader.header(), "DADA");
        } else {
            const formats::guppi::GuppiReader reader(path);
            baseband_fields(fields, reader.header(), "GUPPI");
        }
        return fields;
    }

    Input input(path);
    if (is_presto_path(path)) {
        presto_fields(fields, TimeSeriesReader(path).header());
    } else {
        std::string_view format = "SIGPROC";
        if (input.is_psrfits()) {
            format = "PSRFITS";
        } else if (input.header().header_bytes == 0U) {
            format = "FBH5";
        }
        sigproc_fields(fields, input.header(), input.nsamples(), format);
        if (input.is_psrfits()) {
            fields.add("stored_nbits", "Stored bits per sample",
                       input.stored_nbits());
        }
    }
    if (opts.stats) {
        stats_fields(fields, input, opts.stats_samples);
    }
    return fields;
}

void print_text(const Context& ctx, const Fields& fields) {
    std::size_t width = 0;
    for (const auto& field : fields.items()) {
        width = std::max(width, field.label.size());
    }
    for (const auto& field : fields.items()) {
        ctx.out << fmt::format("{:<{}} : {}\n", field.label, width, field.text);
    }
}

void print_json(const Context& ctx, const Fields& fields) {
    ctx.out << "{\n";
    const auto& items = fields.items();
    for (std::size_t i = 0; i < items.size(); ++i) {
        const auto& field = items[i];
        ctx.out << "  " << json_string(field.key) << ": "
                << (field.numeric ? field.text : json_string(field.text))
                << (i + 1 < items.size() ? ",\n" : "\n");
    }
    ctx.out << "}\n";
}

void run_inspect(const InspectOptions& opts, const Context& ctx) {
    const auto fields = collect(opts);
    if (opts.key.empty()) {
        if (opts.json) {
            print_json(ctx, fields);
        } else {
            print_text(ctx, fields);
        }
        return;
    }
    auto key = lower(opts.key);
    for (const auto& [alias, name] : kKeyAliases) {
        if (key == alias) {
            key = name;
        }
    }
    for (const auto& field : fields.items()) {
        if (field.key == key) {
            ctx.out << field.text << '\n';
            return;
        }
    }
    std::string known;
    for (const auto& field : fields.items()) {
        known += known.empty() ? "" : ", ";
        known += field.key;
    }
    throw ValidationError(
        fmt::format("psrio: key '{}' is not in this file; available keys: {}",
                    opts.key, known));
}

} // namespace

void add_inspect(CLI::App& app, Context& ctx) {
    const auto opts = std::make_shared<InspectOptions>();
    auto* sub       = app.add_subcommand(
        "inspect", "Print the header of a data file (SIGPROC 'header')");
    sub->alias("header");
    sub->add_option("inputs", opts->inputs,
                    "Data files in time order. One file may be .fil, .tim, "
                    ".dat/.inf, .sf/.fits, .dada, .raw, or .h5; several files "
                    "must be .fil or .sf/.fits")
        ->required()
        ->expected(1, -1)
        ->check(CLI::ExistingFile);
    sub->add_option("-k,--key", opts->key,
                    "Print only the value of this key (e.g. tsamp, nchans)");
    sub->add_flag("--stats", opts->stats,
                  "Add statistics of the first samples");
    sub->add_option("--stats-samples", opts->stats_samples,
                    "Time samples used by --stats")
        ->capture_default_str()
        ->check(CLI::PositiveNumber);
    sub->add_flag("--json", opts->json, "Print the header as JSON")
        ->excludes("--key");
    sub->callback([opts, &ctx] { run_inspect(*opts, ctx); });
}

} // namespace psrio::cli
