#pragma once

#include "psrio/astro.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/mmap.hpp"
#include "psrio/formats/sigproc/header.hpp"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <string>
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
    [[nodiscard]] static Header from_inffile(const std::filesystem::path& path);
    [[nodiscard]] static Header from_file(const std::filesystem::path& path);

    // Serialization
    void make_inf(const std::filesystem::path& path) const;
};

inline Header Header::from_sigproc(const std::filesystem::path& path) {
    const detail::MappedFile file(path);
    const auto raw =
        formats::sigproc::FilterbankHeader::parse(file.bytes(), path.string());
    Header hdr;
    hdr.filename    = path.string();
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

inline void Header::make_inf(const std::filesystem::path& path) const {
    std::ofstream inf_out(path);
    if (!inf_out.is_open()) {
        throw IoError(
            std::format("psrio: cannot create .inf file '{}'", path.string()));
    }
    const std::string base = path.stem().string();
    inf_out << std::format(
        " Data file name without the suffix          =  {}\n", base);
    inf_out << std::format(
        " Telescope used                             =  {}\n", telescope);
    inf_out << std::format(
        " Instrument used                            =  {}\n", backend);
    inf_out << std::format(
        " Object being observed                      =  {}\n", source);
    inf_out << std::format(
        " J2000 Right Ascension (hh:mm:ss.ssss)      =  {}\n", ra);
    inf_out << std::format(
        " J2000 Declination     (dd:mm:ss.ssss)      =  {}\n", dec);
    inf_out << std::format(
        " Epoch of observation (MJD)                 =  {:.12f}\n", tstart);
    inf_out << std::format(
        " Barycentered?           (1=yes, 0=no)      =  {}\n", barycentric);
    inf_out << std::format(
        " Number of bins in the data                 =  {}\n", nsamples);
    inf_out << std::format(
        " Width of each time bin (s)                 =  {:.12e}\n", tsamp);
    inf_out << std::format(
        " Dispersion measure (cm-3 pc)               =  {:.5f}\n", dm);
    inf_out << std::format(
        " Central freq of low channel (Mhz)          =  {:.6f}\n", fch1);
    inf_out << std::format(
        " Total bandwidth (Mhz)                      =  {:.6f}\n", bandwidth());
    inf_out << std::format(
        " Number of channels                         =  {}\n", nchans);
    inf_out << std::format(
        " Channel bandwidth (Mhz)                    =  {:.6f}\n", foff);
}

} // namespace psrio
