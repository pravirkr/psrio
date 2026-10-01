#pragma once

#include "psrio/astro.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/header.hpp"

#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>

namespace psrio::formats::presto {

namespace detail {

[[nodiscard]] inline std::string_view trim_view(std::string_view str) {
    constexpr std::string_view kWhitespace = " \t\r\n";
    const auto first                       = str.find_first_not_of(kWhitespace);
    if (first == std::string_view::npos) {
        return {};
    }
    const auto last = str.find_last_not_of(kWhitespace);
    return str.substr(first, last - first + 1);
}

[[nodiscard]] inline double parse_double(std::string_view str,
                                         double default_val = 0.0) {
    const auto s = trim_view(str);
    if (s.empty()) {
        return default_val;
    }
    double val{0.0};
    // NOLINTNEXTLINE(bugprone-suspicious-stringview-data-usage)
    const auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), val);
    if (ec == std::errc{}) {
        return val;
    }
    // Fallback for scientific notation with lowercase/uppercase e or platform
    // differences
    char* end_ptr = nullptr;
    const std::string null_terminated(s);
    val = std::strtod(null_terminated.c_str(), &end_ptr);
    if (end_ptr != null_terminated.c_str()) {
        return val;
    }
    return default_val;
}

[[nodiscard]] inline std::uint64_t parse_u64(std::string_view str,
                                             std::uint64_t default_val = 0) {
    const auto s = trim_view(str);
    if (s.empty()) {
        return default_val;
    }
    std::uint64_t val{0};
    // NOLINTNEXTLINE(bugprone-suspicious-stringview-data-usage)
    const auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), val);
    return (ec == std::errc{}) ? val : default_val;
}

[[nodiscard]] inline int parse_int(std::string_view str, int default_val = 0) {
    const auto s = trim_view(str);
    if (s.empty()) {
        return default_val;
    }
    int val{0};
    // NOLINTNEXTLINE(bugprone-suspicious-stringview-data-usage)
    const auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), val);
    return (ec == std::errc{}) ? val : default_val;
}

[[nodiscard]] inline bool starts_with_ci(std::string_view str,
                                         std::string_view prefix) noexcept {
    if (str.size() < prefix.size()) {
        return false;
    }
    for (std::size_t i = 0; i < prefix.size(); ++i) {
        const char c1 = (str[i] >= 'A' && str[i] <= 'Z')
                            ? static_cast<char>(str[i] + ('a' - 'A'))
                            : str[i];
        const char c2 = (prefix[i] >= 'A' && prefix[i] <= 'Z')
                            ? static_cast<char>(prefix[i] + ('a' - 'A'))
                            : prefix[i];
        if (c1 != c2) {
            return false;
        }
    }
    return true;
}

} // namespace detail

/**
 * @brief Parse a PRESTO ASCII `.inf` file contents into a Header.
 *
 * @param content The text contents of the `.inf` file.
 * @param fname Optional filename to store in the header.
 * @return Parsed Header.
 * @throws FormatError If required fields cannot be parsed.
 */
[[nodiscard]] inline Header parse_inf(std::string_view content,
                                      std::string_view fname = "") {
    Header hdr;
    hdr.filename  = std::string(fname);
    hdr.data_type = "time series";
    hdr.nchans    = 1;
    hdr.nifs      = 1;
    hdr.nbits     = 32; // PRESTO .dat files are 32-bit floats

    std::size_t pos = 0;
    while (pos < content.size()) {
        const auto next_newline = content.find('\n', pos);
        const auto line_end     = (next_newline == std::string_view::npos)
                                      ? content.size()
                                      : next_newline;
        const auto line =
            detail::trim_view(content.substr(pos, line_end - pos));
        pos = (next_newline == std::string_view::npos) ? content.size()
                                                       : next_newline + 1;

        if (line.empty() || line.front() == '#') {
            continue;
        }

        const auto eq_pos = line.rfind('=');
        if (eq_pos == std::string_view::npos) {
            continue;
        }

        const auto key = detail::trim_view(line.substr(0, eq_pos));
        const auto val = detail::trim_view(line.substr(eq_pos + 1));

        if (detail::starts_with_ci(key, "Data file name without")) {
            if (hdr.filename.empty()) {
                hdr.filename = std::string(val);
            }
        } else if (detail::starts_with_ci(key, "Telescope used")) {
            hdr.telescope = std::string(val);
        } else if (detail::starts_with_ci(key, "Instrument used")) {
            hdr.backend = std::string(val);
        } else if (detail::starts_with_ci(key, "Object being observed")) {
            hdr.source = std::string(val);
        } else if (detail::starts_with_ci(key, "J2000 Right Ascension")) {
            hdr.ra = std::string(val);
            try {
                hdr.raj = astro::string_to_packed(val, "RA");
            } catch (...) {
                hdr.raj = 0.0;
            }
        } else if (detail::starts_with_ci(key, "J2000 Declination")) {
            hdr.dec = std::string(val);
            try {
                hdr.dej = astro::string_to_packed(val, "Dec");
            } catch (...) {
                hdr.dej = 0.0;
            }
        } else if (detail::starts_with_ci(key, "Epoch of observation")) {
            hdr.tstart = detail::parse_double(val);
        } else if (detail::starts_with_ci(key, "Barycentered?")) {
            hdr.barycentric = detail::parse_int(val);
        } else if (detail::starts_with_ci(key, "Number of bins")) {
            hdr.nsamples = detail::parse_u64(val);
        } else if (detail::starts_with_ci(key, "Width of each")) {
            hdr.tsamp = detail::parse_double(val);
        } else if (detail::starts_with_ci(key, "Dispersion measure")) {
            hdr.dm = detail::parse_double(val);
        } else if (detail::starts_with_ci(key, "Central freq of low channel")) {
            hdr.fch1 = detail::parse_double(val);
        } else if (detail::starts_with_ci(key, "Channel bandwidth")) {
            hdr.foff = detail::parse_double(val);
        }
    }

    return hdr;
}

/**
 * @brief Parse a PRESTO ASCII `.inf` file from disk.
 *
 * @param path Path to the `.inf` file.
 * @return Parsed Header.
 * @throws IoError If the file cannot be opened.
 */
[[nodiscard]] inline Header parse_inf_file(const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw IoError(std::format("psrio: cannot open PRESTO .inf file '{}'",
                                  path.string()));
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    return parse_inf(buffer.str(), path.stem().string());
}

} // namespace psrio::formats::presto

namespace psrio {

inline Header Header::from_inffile(const std::filesystem::path& path) {
    return formats::presto::parse_inf_file(path);
}

} // namespace psrio
