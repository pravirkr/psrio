#pragma once

#include "psrio/astro.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/mmap.hpp"
#include "psrio/formats/sigproc/header.hpp"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <iterator>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace psrio {

/// Format-specific or unmapped vendor metadata key.
struct ExtraKey {
    std::string name;
    std::string type;
    std::variant<std::string, double, std::int64_t> value;
};

/**
 * @brief Unified observational metadata header for time-series and filterbank
 * data.
 *
 * Normalizes physical, celestial, timing, and telescope parameters across
 * SIGPROC, PRESTO, and other formats into a single idiomatic C++ struct.
 */
struct Header {
    // Physical & Timing
    double tsamp{0.0};
    double tstart{0.0};
    std::uint64_t nsamples{0};
    int nbits{32};
    int nchans{1};
    int nifs{1};
    double fch1{0.0};
    double foff{0.0};
    double dm{0.0};
    int barycentric{0};

    // Source & Instrument
    std::string source{"Fake"};
    std::string telescope{"Fake"};
    std::string backend{"FAKE"};
    std::string filename;
    std::string data_type{"time series"};

    // Celestial Coordinates
    double raj{0.0}; // Packed sexagesimal HHMMSS.ss
    double dej{0.0}; // Packed sexagesimal [+|-]DDMMSS.ss
    std::string ra{"00:00:00.0000"};
    std::string dec{"+00:00:00.0000"};

    // Pointing & Beam
    double az_start{0.0};
    double za_start{0.0};
    double period{0.0};
    int ibeam{0};
    int nbeams{0};

    // Auxiliary / unknown vendor keys
    std::vector<ExtraKey> extra_keys;

    // Derived properties & astronomical convenience helpers
    [[nodiscard]] double tobs() const noexcept {
        return static_cast<double>(nsamples) * tsamp;
    }

    [[nodiscard]] double observation_duration() const noexcept {
        return tobs();
    }

    [[nodiscard]] std::string duration_string() const {
        return astro::format_duration(tobs());
    }

    [[nodiscard]] std::string gregorian_date() const {
        return astro::mjd_to_gregorian(tstart);
    }

    [[nodiscard]] double ra_hours() const noexcept {
        return astro::ra_to_hours(raj);
    }

    [[nodiscard]] double ra_degrees() const noexcept {
        return astro::ra_to_degrees(raj);
    }

    [[nodiscard]] double dec_degrees() const noexcept {
        return astro::dec_to_degrees(dej);
    }

    [[nodiscard]] double bandwidth() const noexcept {
        return std::abs(foff) * static_cast<double>(nchans);
    }

    [[nodiscard]] double ftop() const noexcept { return fch1 - (0.5 * foff); }

    [[nodiscard]] double fbottom() const noexcept {
        return ftop() + (foff * static_cast<double>(nchans));
    }

    [[nodiscard]] double fcenter() const noexcept {
        return 0.5 * (ftop() + fbottom());
    }

    [[nodiscard]] double center_frequency() const noexcept { return fcenter(); }

    [[nodiscard]] std::string basename() const {
        return std::filesystem::path(filename).stem().string();
    }

    [[nodiscard]] std::string extension() const {
        return std::filesystem::path(filename).extension().string();
    }

    // Factory methods
    [[nodiscard]] static Header from_sigproc(const std::filesystem::path& path);
    /// Build from an already parsed SIGPROC header, without touching the file.
    [[nodiscard]] static Header
    from_sigproc(const formats::sigproc::FilterbankHeader& raw,
                 std::string_view filename = "");
    [[nodiscard]] static Header from_inffile(const std::filesystem::path& path);
    [[nodiscard]] static Header from_file(const std::filesystem::path& path);

    /// SIGPROC header with the same metadata.
    ///
    /// data_type is 2 for a single-channel time series and 1 otherwise.
    /// telescope and backend map to ids only when the name is listed in the
    /// SIGPROC tables. refdm is set when dm is non-zero, period when positive,
    /// nbeams when positive. extra_keys are not carried over.
    [[nodiscard]] formats::sigproc::FilterbankHeader to_sigproc() const;

    /// Write a PRESTO `.inf` file in the field order PRESTO's readinf expects.
    ///
    /// The "Central freq of low channel" is the lowest channel centre and the
    /// channel bandwidth is positive, whichever way foff points. nchans is the
    /// number of channels summed into the series.
    /// @throws IoError if the file cannot be written.
    void make_inf(const std::filesystem::path& path) const;
};

inline Header Header::from_sigproc(const std::filesystem::path& path) {
    const detail::MappedFile file(path);
    const auto path_str = path.string();
    return from_sigproc(
        formats::sigproc::FilterbankHeader::parse(file.bytes(), path_str),
        path_str);
}

inline Header
Header::from_sigproc(const formats::sigproc::FilterbankHeader& raw,
                     std::string_view filename) {
    Header hdr;
    hdr.filename    = std::string(filename);
    hdr.tsamp       = raw.tsamp.value_or(0.0);
    hdr.tstart      = raw.tstart.value_or(0.0);
    hdr.nbits       = raw.nbits;
    hdr.nchans      = raw.nchans;
    hdr.nifs        = raw.nifs;
    hdr.fch1        = raw.fch1.value_or(0.0);
    hdr.foff        = raw.foff.value_or(0.0);
    hdr.dm          = raw.refdm.value_or(0.0);
    hdr.barycentric = raw.barycentric.value_or(0);
    hdr.source      = raw.source_name.value_or("Fake");
    hdr.telescope   = std::string(raw.telescope_name());
    hdr.backend     = std::string(raw.machine_name());
    hdr.raj         = raw.src_raj.value_or(0.0);
    hdr.dej         = raw.src_dej.value_or(0.0);
    hdr.ra          = raw.ra_string();
    hdr.dec         = raw.dec_string();
    hdr.az_start    = raw.az_start.value_or(0.0);
    hdr.za_start    = raw.za_start.value_or(0.0);
    hdr.period      = raw.period.value_or(0.0);
    hdr.ibeam       = raw.ibeam.value_or(0);
    hdr.data_type   = (raw.data_type.has_value() && *raw.data_type == 2)
                          ? "time series"
                          : std::string(raw.data_type_name());

    if (raw.nbits <= 4 && raw.nchans == 1) {
        hdr.nsamples =
            (raw.data_bytes * 8U) / static_cast<std::uint64_t>(raw.nbits);
    } else {
        try {
            hdr.nsamples = raw.nsamples();
        } catch (...) {
            if (raw.nbits >= 8 && raw.nchans == 1) {
                hdr.nsamples = raw.data_bytes /
                               (static_cast<std::uint64_t>(raw.nbits) / 8U);
            } else {
                hdr.nsamples = 0;
            }
        }
    }

    for (const auto& u : raw.extra) {
        ExtraKey ek;
        ek.name = std::string(u.name);
        ek.type = "raw";
        hdr.extra_keys.push_back(std::move(ek));
    }

    return hdr;
}

inline formats::sigproc::FilterbankHeader Header::to_sigproc() const {
    namespace sp = formats::sigproc;
    sp::FilterbankHeader raw;
    raw.nchans       = nchans;
    raw.nifs         = nifs;
    raw.nbits        = nbits;
    raw.tsamp        = tsamp;
    raw.tstart       = tstart;
    raw.fch1         = fch1;
    raw.foff         = foff;
    raw.data_type    = (nchans == 1 && data_type == "time series") ? 2 : 1;
    raw.telescope_id = sp::telescope_id_from_name(telescope);
    raw.machine_id   = sp::machine_id_from_name(backend);
    if (!source.empty()) {
        raw.source_name = source;
    }
    raw.src_raj     = raj;
    raw.src_dej     = dej;
    raw.az_start    = az_start;
    raw.za_start    = za_start;
    raw.barycentric = barycentric;
    raw.ibeam       = ibeam;
    if (dm != 0.0) {
        raw.refdm = dm;
    }
    if (period > 0.0) {
        raw.period = period;
    }
    if (nbeams > 0) {
        raw.nbeams = nbeams;
    }
    return raw;
}

inline void Header::make_inf(const std::filesystem::path& path) const {
    std::ofstream inf_out(path);
    if (!inf_out.is_open()) {
        throw IoError(
            fmt::format("psrio: cannot create .inf file '{}'", path.string()));
    }
    const int channels   = nchans > 0 ? nchans : 1;
    const double chan_bw = std::abs(foff);
    const double low_centre =
        foff < 0.0 ? fch1 + (foff * static_cast<double>(channels - 1)) : fch1;

    std::string text;
    const auto field = [&text](std::string_view label, const auto& value) {
        fmt::format_to(std::back_inserter(text), " {:<39}=  {}\n", label,
                       value);
    };
    field("Data file name without suffix", path.stem().string());
    field("Telescope used", telescope);
    field("Instrument used", backend);
    field("Object being observed", source);
    field("J2000 Right Ascension (hh:mm:ss.ssss)",
          ra.empty() ? "00:00:00.0000" : ra);
    field("J2000 Declination     (dd:mm:ss.ssss)",
          dec.empty() ? "00:00:00.0000" : dec);
    field("Data observed by", "Unknown");
    field("Epoch of observation (MJD)", fmt::format("{:.15f}", tstart));
    field("Barycentered?           (1 yes, 0 no)", barycentric);
    field("Number of bins in the time series", nsamples);
    field("Width of each time series bin (sec)", fmt::format("{:.15g}", tsamp));
    field("Any breaks in the data? (1 yes, 0 no)", 0);
    field("Type of observation (EM band)", "Radio");
    field("Beam diameter (arcsec)", 0);
    field("Dispersion measure (cm-3 pc)", fmt::format("{:.12g}", dm));
    field("Central freq of low channel (MHz)",
          fmt::format("{:.12g}", low_centre));
    field("Total bandwidth (MHz)",
          fmt::format("{:.12g}", chan_bw * static_cast<double>(channels)));
    field("Number of channels", channels);
    field("Channel bandwidth (MHz)", fmt::format("{:.12g}", chan_bw));
    field("Data analyzed by", "psrio");
    text.append(" Any additional notes:\n    Written by psrio.\n\n");

    inf_out << text;
    inf_out.close();
    if (!inf_out) {
        throw IoError(fmt::format("psrio: failed to write .inf file '{}'",
                                  path.string()));
    }
}

} // namespace psrio
