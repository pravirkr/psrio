#pragma once

#include <array>
#include <charconv>
#include <cmath>
#include <format>
#include <numbers>
#include <stdexcept>
#include <string>
#include <string_view>

namespace psrio::astro {

/// Dispersion constant in MHz^2 pc^-1 cm^3 s.
/// Standard value widely used in pulsar and FRB astronomy:
/// dt = kDispersionConstant * DM * (1 / (f_low^2) - 1 / (f_high^2))
inline constexpr double kDispersionConstant = 4.148808e3;

/// Sign-and-magnitude representation of a `[+|-]DD:MM:SS.sss` coordinate.
/// For Right Ascension, major is hours in [0, 24).
/// For Declination, major is degrees in [0, 90].
struct Sexagesimal {
    bool negative{false};
    unsigned major{0};
    unsigned minutes{0};
    double seconds{0.0};
};

namespace detail {

inline constexpr std::string_view kWhitespace = " \t\n\v\f\r";

[[nodiscard]] constexpr std::string_view trim(std::string_view str) {
    const auto first = str.find_first_not_of(kWhitespace);
    if (first == std::string_view::npos) {
        return {};
    }
    const auto last = str.find_last_not_of(kWhitespace);
    return str.substr(first, last - first + 1);
}

[[nodiscard]] constexpr bool is_digit(char chr) noexcept {
    return chr >= '0' && chr <= '9';
}

[[nodiscard]] constexpr bool is_plain_decimal(std::string_view field) noexcept {
    if (field.empty() || !is_digit(field.front())) {
        return false;
    }
    bool seen_dot = false;
    for (const char chr : field) {
        if (chr == '.') {
            if (seen_dot) {
                return false;
            }
            seen_dot = true;
        } else if (!is_digit(chr)) {
            return false;
        }
    }
    return true;
}

} // namespace detail

/// Parse a sexagesimal coordinate string `[+|-]DD:MM:SS.sss`.
/// Surrounding whitespace is permitted.
/// @throws std::invalid_argument on malformed input or out-of-range fields.
[[nodiscard]] inline Sexagesimal
parse_sexagesimal(std::string_view input,
                  std::string_view what = "coordinate") {
    const auto invalid = [&] {
        return std::invalid_argument(std::format(
            "psrio: invalid {} format: '{}' (expected [+|-]dd:mm:ss.sss)", what,
            input));
    };

    std::string_view str = detail::trim(input);
    if (str.empty()) {
        throw invalid();
    }

    Sexagesimal out;
    if (str.front() == '+' || str.front() == '-') {
        out.negative = (str.front() == '-');
        str.remove_prefix(1);
    }

    const char* const end = str.data() + str.size();
    // NOLINTNEXTLINE(bugprone-suspicious-stringview-data-usage)
    const auto [p1, e1] = std::from_chars(str.data(), end, out.major);
    if (e1 != std::errc{} || p1 == end || *p1 != ':') {
        throw invalid();
    }

    const auto [p2, e2] = std::from_chars(p1 + 1, end, out.minutes);
    if (e2 != std::errc{} || p2 == end || *p2 != ':') {
        throw invalid();
    }

    const std::string_view sec_field(p2 + 1, end);
    if (!detail::is_plain_decimal(sec_field)) {
        throw invalid();
    }
    // NOLINTNEXTLINE(bugprone-suspicious-stringview-data-usage)
    const auto [p3, e3] = std::from_chars(sec_field.data(), end, out.seconds,
                                          std::chars_format::fixed);
    if (e3 != std::errc{} || p3 != end) {
        throw invalid();
    }

    if (out.minutes >= 60 || out.seconds >= 60.0) {
        throw invalid();
    }
    return out;
}

/// Convert a Sexagesimal coordinate struct to SIGPROC packed coordinate format
/// (hhmmss.ss or [+|-]ddmmss.ss).
[[nodiscard]] constexpr double
sexagesimal_to_packed(const Sexagesimal& sex) noexcept {
    const double val = (static_cast<double>(sex.major) * 10000.0) +
                       (static_cast<double>(sex.minutes) * 100.0) + sex.seconds;
    return sex.negative ? -val : val;
}

/// Parse a sexagesimal coordinate string and return it in SIGPROC packed
/// format.
[[nodiscard]] inline double
string_to_packed(std::string_view str, std::string_view what = "coordinate") {
    return sexagesimal_to_packed(parse_sexagesimal(str, what));
}

/// Convert hours, minutes, and seconds to radians.
[[nodiscard]] inline double
hms_to_rad(int hour, int minutes, double sec) noexcept {
    const double sign = hour < 0 ? -1.0 : 1.0;
    const double total_seconds =
        (std::abs(hour) * 3600.0) + (std::abs(minutes) * 60.0) + std::abs(sec);
    return sign * total_seconds * std::numbers::pi_v<double> / 43200.0;
}

/// Convert degrees, arcminutes, and arcseconds to radians.
[[nodiscard]] inline double
dms_to_rad(int deg, int minutes, double sec) noexcept {
    const double sign =
        (deg < 0 || (deg == 0 && (minutes < 0 || sec < 0))) ? -1.0 : 1.0;
    const double total_seconds =
        (std::abs(deg) * 3600.0) + (std::abs(minutes) * 60.0) + std::abs(sec);
    return sign * total_seconds * std::numbers::pi_v<double> / 648000.0;
}

/// Convert decimal degrees to SIGPROC packed coordinate format (ddmmss.ss).
[[nodiscard]] inline double deg_to_dms(double angle) noexcept {
    const int sign   = (std::signbit(angle) || angle < 0.0) ? -1 : 1;
    angle            = std::abs(angle);
    const auto deg   = static_cast<int>(angle);
    angle            = (angle - static_cast<double>(deg)) * 60.0;
    const auto min   = static_cast<int>(angle);
    const double sec = (angle - static_cast<double>(min)) * 60.0;
    return static_cast<double>(sign) *
           ((static_cast<double>(deg) * 10000.0) +
            (static_cast<double>(min) * 100.0) + sec);
}

/// Convert RA string `HH:MM:SS.sss` to radians.
/// @throws std::invalid_argument if outside [00:00:00, 24:00:00).
[[nodiscard]] inline double ra_to_rad(std::string_view ra_string) {
    const Sexagesimal sex = parse_sexagesimal(ra_string, "right ascension");
    if (sex.negative || sex.major >= 24) {
        throw std::invalid_argument(std::format(
            "psrio: right ascension out of range: '{}' (expected 00:00:00 "
            "to 23:59:59.999...)",
            ra_string));
    }
    return hms_to_rad(static_cast<int>(sex.major),
                      static_cast<int>(sex.minutes), sex.seconds);
}

/// Convert Declination string `[+|-]DD:MM:SS.sss` to radians.
/// @throws std::invalid_argument if outside [-90:00:00, +90:00:00].
[[nodiscard]] inline double dec_to_rad(std::string_view dec_string) {
    const Sexagesimal sex = parse_sexagesimal(dec_string, "declination");
    if (sex.major > 90 ||
        (sex.major == 90 && (sex.minutes != 0 || sex.seconds != 0.0))) {
        throw std::invalid_argument(std::format(
            "psrio: declination out of range: '{}' (expected -90:00:00 "
            "to +90:00:00)",
            dec_string));
    }
    const double magnitude =
        dms_to_rad(static_cast<int>(sex.major), static_cast<int>(sex.minutes),
                   sex.seconds);
    return sex.negative ? -magnitude : magnitude;
}

/// Format packed SIGPROC RA (hhmmss.ss) into "HH:MM:SS.ssss".
[[nodiscard]] inline std::string ra_to_string(double src_raj) {
    const double angle = std::abs(src_raj);
    const int hh       = static_cast<int>(angle) / 10000;
    const int mm       = (static_cast<int>(angle) / 100) % 100;
    const double ss    = angle - (static_cast<double>(mm) * 100.0) -
                         (static_cast<double>(hh) * 10000.0);
    return std::format("{:02d}:{:02d}:{:07.4f}", hh, mm, ss);
}

/// Format packed SIGPROC Declination (ddmmss.ss) into "[+|-]DD:MM:SS.ssss".
[[nodiscard]] inline std::string dec_to_string(double src_dej) {
    const bool negative = std::signbit(src_dej) || src_dej < 0.0;
    const double angle  = std::abs(src_dej);
    const int dd        = static_cast<int>(angle) / 10000;
    const int mm        = (static_cast<int>(angle) / 100) % 100;
    const double ss     = angle - (static_cast<double>(mm) * 100.0) -
                          (static_cast<double>(dd) * 10000.0);
    return std::format("{}{:02d}:{:02d}:{:07.4f}", negative ? "-" : "+", dd, mm,
                       ss);
}

/// Convert packed SIGPROC RA (hhmmss.ss) to decimal hours in [0, 24).
[[nodiscard]] inline double ra_to_hours(double src_raj) noexcept {
    const double angle = std::abs(src_raj);
    const int hh       = static_cast<int>(angle) / 10000;
    const int mm       = (static_cast<int>(angle) / 100) % 100;
    const double ss    = angle - (static_cast<double>(mm) * 100.0) -
                         (static_cast<double>(hh) * 10000.0);
    return static_cast<double>(hh) + (static_cast<double>(mm) / 60.0) +
           (ss / 3600.0);
}

/// Convert packed SIGPROC RA (hhmmss.ss) to decimal degrees in [0, 360).
[[nodiscard]] inline double ra_to_degrees(double src_raj) noexcept {
    return ra_to_hours(src_raj) * 15.0;
}

/// Convert packed SIGPROC RA (hhmmss.ss) to radians in [0, 2*pi).
[[nodiscard]] inline double ra_to_radians(double src_raj) noexcept {
    return ra_to_hours(src_raj) * (std::numbers::pi_v<double> / 12.0);
}

/// Convert packed SIGPROC Declination (ddmmss.ss) to decimal degrees in [-90,
/// +90].
[[nodiscard]] inline double dec_to_degrees(double src_dej) noexcept {
    const bool negative = std::signbit(src_dej) || src_dej < 0.0;
    const double angle  = std::abs(src_dej);
    const int dd        = static_cast<int>(angle) / 10000;
    const int mm        = (static_cast<int>(angle) / 100) % 100;
    const double ss     = angle - (static_cast<double>(mm) * 100.0) -
                          (static_cast<double>(dd) * 10000.0);
    const double deg    = static_cast<double>(dd) +
                          (static_cast<double>(mm) / 60.0) + (ss / 3600.0);
    return negative ? -deg : deg;
}

/// Convert packed SIGPROC Declination (ddmmss.ss) to radians in [-pi/2, +pi/2].
[[nodiscard]] inline double dec_to_radians(double src_dej) noexcept {
    return dec_to_degrees(src_dej) * (std::numbers::pi_v<double> / 180.0);
}

/// Convert integer Modified Julian Date (MJD) to Gregorian calendar date
/// "YYYY-MM-DD".
[[nodiscard]] inline std::string mjd_to_gregorian(int mjd) {
    const int jd    = mjd + 2400001;
    const int a     = jd + 32044;
    const int b     = ((4 * a) + 3) / 146097;
    const int c     = a - ((146097 * b) / 4);
    const int d     = ((4 * c) + 3) / 1461;
    const int e     = c - ((1461 * d) / 4);
    const int m     = ((5 * e) + 2) / 153;
    const int day   = e - (((153 * m) + 2) / 5) + 1;
    const int month = m + 3 - (12 * (m / 10));
    const int year  = (100 * b) + d - 4800 + (m / 10);
    return std::format("{:04d}-{:02d}-{:02d}", year, month, day);
}

/// Convert double MJD (fractional day) to Gregorian calendar date "YYYY-MM-DD".
[[nodiscard]] inline std::string mjd_to_gregorian(double mjd) {
    return mjd_to_gregorian(static_cast<int>(std::floor(mjd)));
}

/// Format a duration in seconds into a human-readable string (e.g. "45.0
/// seconds", "12.3 hours").
[[nodiscard]] inline std::string format_duration(double duration_seconds) {
    constexpr std::array kUnits{"seconds", "minutes", "hours", "days"};
    constexpr std::array kFactors{60.0, 60.0, 24.0};

    std::size_t index = 0;
    double val        = duration_seconds;
    for (const double factor : kFactors) {
        if (val < factor) {
            break;
        }
        val /= factor;
        ++index;
    }
    return std::format("{:.1f} {}", val, kUnits[index]);
}

/// Calculate dispersion delay in seconds between @p freq_mhz and @p
/// ref_freq_mhz. When @p ref_freq_mhz is 0.0, the delay is relative to infinite
/// frequency: dt = kDispersionConstant * dm / (freq_mhz^2).
[[nodiscard]] inline double dispersion_delay(
    double dm, double freq_mhz, double ref_freq_mhz = 0.0) noexcept {
    if (freq_mhz <= 0.0) {
        return 0.0;
    }
    const double inv_f2 = 1.0 / (freq_mhz * freq_mhz);
    if (ref_freq_mhz <= 0.0) {
        return kDispersionConstant * dm * inv_f2;
    }
    const double inv_ref2 = 1.0 / (ref_freq_mhz * ref_freq_mhz);
    return kDispersionConstant * dm * (inv_f2 - inv_ref2);
}

} // namespace psrio::astro
