#pragma once

#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/mjd.hpp"
#include "psrio/formats/psrfits/fits_writer.hpp"
#include "psrio/formats/psrfits/table.hpp"

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <exception>
#include <filesystem>
#include <fmt/format.h>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace psrio::formats::psrfits {

/// `DAT_SCL` and `DAT_OFFS` of one quantised profile.
struct ProfileScale {
    float scale{0.0F};
    float offset{0.0F};
};

/**
 * @brief Scale and offset that store one profile as int16.
 *
 * The profile's range [min, max] maps linearly onto the integers
 * [-32751, 32750], the range PSRCHIVE uses for fold-mode DATA, so that
 * `amp ≈ data * scale + offset`:
 * - `offset = (min * 32750 + max * 32751) / 65501`, evaluated in double and
 *   stored as float;
 * - `scale = (max - min) / 65501` (the difference taken in float), or 0 for
 *   a flat profile (range below 100 * FLT_MIN).
 *
 * Matching PSRCHIVE's precision choices here is what makes DAT_SCL,
 * DAT_OFFS and DATA byte-identical to dspsr output.
 *
 * @return The scale and offset to store in `DAT_SCL` and `DAT_OFFS`.
 */
inline ProfileScale profile_scale(std::span<const float> amps) {
    constexpr double kLowest  = 32751.0; // magnitude of the lowest level
    constexpr double kHighest = 32750.0;
    constexpr double kLevels  = kLowest + kHighest;
    float lo                  = 0.0F;
    float hi                  = 0.0F;
    if (!amps.empty()) {
        lo = amps.front();
        hi = amps.front();
    }
    for (const float value : amps) {
        // A NaN never replaces the running extremes.
        hi = std::max(hi, value);
        lo = std::min(lo, value);
    }
    ProfileScale result;
    result.offset = static_cast<float>(((static_cast<double>(lo) * kHighest) +
                                        (static_cast<double>(hi) * kLowest)) /
                                       kLevels);
    if (std::fabs(lo - hi) > (100.0 * FLT_MIN)) {
        result.scale =
            static_cast<float>(static_cast<double>(hi - lo) / kLevels);
    }
    return result;
}

/**
 * @brief Quantise @p amps with a scale and offset from profile_scale().
 *
 * Each value is `round((amp - offset) / scale)` in float, clamped to int16.
 * A flat profile (scale 0) stores 0 where it equals the offset and the
 * nearest int16 limit otherwise, so the result is always defined.
 */
inline void quantise_with(std::span<const float> amps,
                          const ProfileScale& scale,
                          std::span<std::int16_t> out) {
    if (out.size() != amps.size()) {
        throw ValidationError("psrio: quantised profile size mismatch");
    }
    constexpr auto kLow =
        static_cast<float>(std::numeric_limits<std::int16_t>::min());
    constexpr auto kHigh =
        static_cast<float>(std::numeric_limits<std::int16_t>::max());
    for (std::size_t bin = 0; bin < amps.size(); ++bin) {
        const float diff = amps[bin] - scale.offset;
        float value      = 0.0F;
        if (scale.scale != 0.0F) {
            value = std::round(diff / scale.scale);
        } else if (diff != 0.0F && !std::isnan(diff)) {
            value = diff > 0.0F ? kHigh : kLow;
        } else {
            out[bin] = 0;
            continue;
        }
        if (std::isnan(value)) {
            out[bin] = 0;
            continue;
        }
        value    = std::min(std::max(value, kLow), kHigh);
        out[bin] = static_cast<std::int16_t>(value);
    }
}

inline ProfileScale quantise_profile(std::span<const float> amps,
                                     std::span<std::int16_t> out) {
    if (out.size() != amps.size()) {
        throw ValidationError("psrio: quantised profile size mismatch");
    }
    const auto scale = profile_scale(amps);
    quantise_with(amps, scale, out);
    return scale;
}

/// Everything write_fold_archive() needs besides the profiles.
struct FoldArchiveSpec {
    /// Primary header of the search-mode input; supplies the observation
    /// keywords copied into the archive.
    detail::Cards input_primary;
    /// SUBINT header of the input; supplies TBIN, NBITS, NSBLK, ZERO_OFF,
    /// SIGNINT, and NCH_STRT.
    detail::Cards input_subint;
    /// Start of the first folded sample (STT_IMJD/SMJD/OFFS).
    ::psrio::detail::Mjd start;
    /// Folding period and length of the single sub-integration, seconds.
    double period{0.0};
    double tsubint{0.0};
    /// Sub-integration epoch minus start, seconds.
    double offs_sub{0.0};
    std::uint32_t nbin{0};
    std::uint32_t nchan{0};
    std::uint32_t npol{0};
    std::string pol_type;
    double dm{0.0};
    bool dedispersed{false};
    double centre_frequency{0.0};
    double bandwidth{0.0};
    /// Centre frequency of each channel (DAT_FREQ), MHz.
    std::vector<double> frequencies;
    /// Recorded in the HISTORY table's PROC_CMD (truncated to 256).
    std::string command;
    /// Creation time for DATE and DATE_PRO.
    std::time_t created{0};
};

namespace detail {

enum class KeyKind : std::uint8_t {
    kString,
    kInteger,
    kDouble,
    kStar,  ///< always the PSRFITS "unset" string '*'
    kEmpty, ///< always an empty string
};

struct PrimaryKey {
    std::string_view key;
    KeyKind kind;
    std::string_view comment;
};

[[nodiscard]] constexpr PrimaryKey
primary_key(std::string_view key, KeyKind kind, std::string_view comment) {
    return PrimaryKey{
        .key     = key,
        .kind    = kind,
        .comment = comment,
    };
}

/// Primary keywords of a PSRFITS fold-mode archive, after SIMPLE..EXTEND,
/// in the order PSRCHIVE writes them so headers can be compared card by
/// card.
inline constexpr auto kPrimaryKeys = std::to_array<PrimaryKey>({
    primary_key("OBSERVER", KeyKind::kString, "Observer"),
    primary_key("PROJID", KeyKind::kString, "Project"),
    primary_key("TELESCOP", KeyKind::kString, "Telescope"),
    primary_key("ANT_X", KeyKind::kDouble, "[m] ITRF antenna position, X"),
    primary_key("ANT_Y", KeyKind::kDouble, "[m] ITRF antenna position, Y"),
    primary_key("ANT_Z", KeyKind::kDouble, "[m] ITRF antenna position, Z"),
    primary_key("FRONTEND", KeyKind::kString, "Receiver"),
    primary_key("IBEAM", KeyKind::kString, "Beam of a multibeam receiver"),
    primary_key("NRCVR", KeyKind::kInteger, "Receiver polarisation channels"),
    primary_key("FD_POLN", KeyKind::kString, "Feed basis, LIN or CIRC"),
    primary_key("FD_HAND", KeyKind::kInteger, "Feed handedness, +1 or -1"),
    primary_key("FD_SANG", KeyKind::kDouble, "[deg] Feed symmetry angle"),
    primary_key(
        "FD_XYPH", KeyKind::kDouble, "[deg] Cross phase of the noise diode"),
    primary_key("BACKEND", KeyKind::kString, "Backend"),
    primary_key("BECONFIG", KeyKind::kString, "Backend configuration"),
    primary_key(
        "BE_PHASE", KeyKind::kInteger, "Backend cross-product sign, 0 unknown"),
    primary_key(
        "BE_DCC", KeyKind::kInteger, "1 if downconversion is conjugated"),
    primary_key("BE_DELAY", KeyKind::kDouble, "[s] Backend delay"),
    primary_key("TCYCLE", KeyKind::kDouble, "[s] Online integration cycle"),
    primary_key(
        "OBS_MODE", KeyKind::kString, "Observation mode: PSR, CAL or SEARCH"),
    primary_key("DATE-OBS", KeyKind::kString, "UTC start of the observation"),
    primary_key(
        "OBSFREQ", KeyKind::kDouble, "[MHz] Observing centre frequency"),
    primary_key("OBSBW", KeyKind::kDouble, "[MHz] Observing bandwidth"),
    primary_key("OBSNCHAN", KeyKind::kInteger, "Channels recorded"),
    primary_key(
        "CHAN_DM", KeyKind::kDouble, "[cm-3 pc] DM of any online dedispersion"),
    primary_key("PNT_ID", KeyKind::kString, "Pointing centre"),
    primary_key("SRC_NAME", KeyKind::kString, "Source"),
    primary_key("COORD_MD", KeyKind::kString, "Coordinate system"),
    primary_key("EQUINOX", KeyKind::kDouble, "Equinox of the coordinates"),
    primary_key("RA", KeyKind::kString, "Right ascension, hh:mm:ss"),
    primary_key("DEC", KeyKind::kString, "Declination, dd:mm:ss"),
    primary_key("BMAJ", KeyKind::kDouble, "[deg] Beam major axis"),
    primary_key("BMIN", KeyKind::kDouble, "[deg] Beam minor axis"),
    primary_key("BPA", KeyKind::kDouble, "[deg] Beam orientation"),
    primary_key("STT_CRD1", KeyKind::kString, "First coordinate at start"),
    primary_key("STT_CRD2", KeyKind::kString, "Second coordinate at start"),
    primary_key("TRK_MODE", KeyKind::kString, "Telescope tracking mode"),
    primary_key("STP_CRD1", KeyKind::kString, "First coordinate at stop"),
    primary_key("STP_CRD2", KeyKind::kString, "Second coordinate at stop"),
    primary_key("SCANLEN", KeyKind::kStar, "[s] Planned scan length"),
    primary_key("FD_MODE", KeyKind::kString, "Feed rotation mode"),
    primary_key("FA_REQ", KeyKind::kDouble, "[deg] Requested feed angle"),
    primary_key("CAL_MODE", KeyKind::kEmpty, "Noise diode mode"),
    primary_key("CAL_FREQ", KeyKind::kStar, "[Hz] Noise diode switching rate"),
    primary_key("CAL_DCYC", KeyKind::kStar, "Noise diode duty cycle"),
    primary_key("CAL_PHS", KeyKind::kStar, "Noise diode phase at start"),
    primary_key("CAL_NPHS", KeyKind::kStar, "Noise diode states per cycle"),
});

/// Reformat a sexagesimal angle string with three decimal places in the
/// last field (the precision PSRCHIVE stores). Returns @p text unchanged
/// if it does not parse.
[[nodiscard]] inline std::string sexagesimal_3dp(std::string_view text,
                                                 bool is_signed) {
    text                  = trim(text);
    std::string_view body = text;
    char sign             = '+';
    if (!body.empty() && (body.front() == '-' || body.front() == '+')) {
        sign = body.front();
        body.remove_prefix(1);
    }
    std::array<double, 3> parts{};
    std::size_t field = 0;
    std::size_t start = 0;
    for (std::size_t index = 0; index <= body.size(); ++index) {
        if (index == body.size() || body[index] == ':') {
            if (field >= parts.size()) {
                return std::string{text};
            }
            const auto token = body.substr(start, index - start);
            try {
                parts[field] = std::stod(std::string{token});
            } catch (const std::exception&) {
                return std::string{text};
            }
            ++field;
            start = index + 1;
        }
    }
    if (field != 3U) {
        return std::string{text};
    }
    auto milli       = static_cast<std::int64_t>(std::llround(
        ((parts[0] * 3600.0) + (parts[1] * 60.0) + parts[2]) * 1000.0));
    const auto hours = milli / 3600000;
    milli -= hours * 3600000;
    const auto minutes = milli / 60000;
    milli -= minutes * 60000;
    const auto seconds = static_cast<double>(milli) / 1000.0;
    if (is_signed) {
        return fmt::format("{}{:02d}:{:02d}:{:06.3f}", sign, hours, minutes,
                           seconds);
    }
    return fmt::format("{:02d}:{:02d}:{:06.3f}", hours, minutes, seconds);
}

[[nodiscard]] inline Card primary_card(const PrimaryKey& spec,
                                       const Cards& input) {
    const std::string key{spec.key};
    const std::string comment{spec.comment};
    const auto value = input.find(spec.key);
    switch (spec.kind) {
    case KeyKind::kStar:
        return Card::string(key, "*", comment);
    case KeyKind::kEmpty:
        return Card::string(key, "", comment);
    case KeyKind::kString:
        return Card::string(key, value.value_or(""), comment);
    case KeyKind::kInteger:
        if (const auto number = input.find_number(spec.key)) {
            return Card::integer(key, std::llround(*number), comment);
        }
        return Card::string(key, "*", comment);
    case KeyKind::kDouble:
        if (const auto number = input.find_number(spec.key)) {
            return Card::real(key, *number, comment);
        }
        return Card::string(key, "*", comment);
    }
    return Card::string(key, "*", comment);
}

[[nodiscard]] inline std::string utc_stamp(std::time_t when) {
    std::tm parts{};
    gmtime_r(&when, &parts);
    return fmt::format("{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}",
                       parts.tm_year + 1900, parts.tm_mon + 1, parts.tm_mday,
                       parts.tm_hour, parts.tm_min, parts.tm_sec);
}

/// Local time in `ctime()` layout without the newline (HISTORY DATE_PRO).
[[nodiscard]] inline std::string local_ctime(std::time_t when) {
    std::tm parts{};
    localtime_r(&when, &parts);
    static constexpr std::array<std::string_view, 7> kDays{
        "Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat",
    };
    static constexpr std::array<std::string_view, 12> kMonths{
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    };
    return fmt::format("{} {} {:2d} {:02d}:{:02d}:{:02d} {}",
                       kDays.at(static_cast<std::size_t>(parts.tm_wday)),
                       kMonths.at(static_cast<std::size_t>(parts.tm_mon)),
                       parts.tm_mday, parts.tm_hour, parts.tm_min, parts.tm_sec,
                       parts.tm_year + 1900);
}

[[nodiscard]] inline std::vector<Card>
primary_header(const FoldArchiveSpec& spec) {
    const auto& input = spec.input_primary;
    std::vector<Card> cards;
    cards.push_back(Card::logical("SIMPLE", true, "Standard FITS file"));
    cards.push_back(Card::integer("BITPIX", 8, "No primary data array"));
    cards.push_back(Card::integer("NAXIS", 0, "No primary data array"));
    cards.push_back(Card::logical("EXTEND", true, "Extensions follow"));
    cards.push_back(Card::note(
        "COMMENT", "  PSRFITS fold-mode archive written by psrio cutout."));
    cards.push_back(Card::string("HDRVER", "7.1", "PSRFITS header version"));
    cards.push_back(
        Card::string("FITSTYPE", "PSRFITS", "Pulsar FITS convention"));
    cards.push_back(Card::string("DATE", utc_stamp(spec.created),
                                 "UTC time this file was written"));
    for (const auto& key : kPrimaryKeys) {
        if (key.key == "OBS_MODE") {
            cards.push_back(
                Card::string("OBS_MODE", "PSR", std::string{key.comment}));
        } else if (key.key == "COORD_MD") {
            cards.push_back(
                Card::string("COORD_MD", "J2000", std::string{key.comment}));
        } else if (key.key == "RA" || key.key == "DEC") {
            const bool dec   = key.key == "DEC";
            const auto value = input.find(key.key);
            cards.push_back(Card::string(
                std::string{key.key},
                value.has_value() ? sexagesimal_3dp(*value, dec) : "",
                std::string{key.comment}));
        } else {
            cards.push_back(primary_card(key, input));
        }
    }
    cards.push_back(Card::integer("STT_IMJD", spec.start.day(),
                                  "Start, whole UTC days (MJD)"));
    cards.push_back(Card::integer("STT_SMJD", spec.start.second(),
                                  "[s] Start, whole seconds of the day"));
    cards.push_back(Card::real("STT_OFFS", spec.start.fraction(),
                               "[s] Start, fraction of a second"));
    cards.push_back(
        primary_card(primary_key("STT_LST", KeyKind::kDouble,
                                 "[s] Local sidereal time at start"),
                     input));
    return cards;
}

/// Fixed-width text cell, space padded.
inline void put_text(std::byte* out, std::size_t width, std::string_view text) {
    for (std::size_t index = 0; index < width; ++index) {
        const char character = index < text.size() ? text[index] : ' ';
        out[index]           = std::byte{static_cast<unsigned char>(character)};
    }
}

} // namespace detail

/**
 * @brief Write a one-sub-integration fold-mode PSRFITS archive.
 *
 * The HDUs, keywords and columns are those of a dspsr archive made from
 * search-mode PSRFITS input, so the two can be compared directly: PRIMARY
 * (observation keywords copied from the input), HISTORY
 * (one row), and SUBINT (INDEXVAL, TSUBINT, OFFS_SUB, PERIOD, AUX_DM,
 * AUX_RM, DAT_FREQ, DAT_WTS, DAT_OFFS, DAT_SCL, DATA). DATA is int16,
 * quantised by quantise_profile(); every weight is 1.
 *
 * @param amps Profiles ordered `[pol][chan][bin]`.
 * @throws ValidationError if @p amps or the spec is inconsistent.
 * @throws IoError if the file cannot be written.
 */
inline void write_fold_archive(const std::filesystem::path& path,
                               const FoldArchiveSpec& spec,
                               std::span<const float> amps) {
    const auto nbin  = static_cast<std::size_t>(spec.nbin);
    const auto nchan = static_cast<std::size_t>(spec.nchan);
    const auto npol  = static_cast<std::size_t>(spec.npol);
    const auto nprof = nchan * npol;
    if (nbin == 0U || nchan == 0U || npol == 0U) {
        throw ValidationError("psrio: an archive needs nbin, nchan, and npol");
    }
    if (amps.size() != nprof * nbin) {
        throw ValidationError(
            fmt::format("psrio: archive needs {} amplitudes, got {}",
                        nprof * nbin, amps.size()));
    }
    if (spec.frequencies.size() != nchan) {
        throw ValidationError("psrio: archive needs one frequency per channel");
    }

    FitsWriter out(path);
    out.write_header(detail::primary_header(spec));

    // HISTORY: one row describing this archive.
    {
        const std::vector<TableColumn> columns{
            column("DATE_PRO", "24A", "", "", "When this row was written"),
            column("PROC_CMD", "256A", "", "", "Command that wrote the file"),
            column("SCALE", "8A", "", "", "Amplitude units"),
            column("POL_TYPE", "8A", "", "", "Stored polarisation products"),
            column("NSUB", "1J", "", "", "Sub-integrations"),
            column("NPOL", "1I", "", "", "Polarisation products"),
            column("NBIN", "1I", "", "", "Phase bins"),
            column("NBIN_PRD", "1I", "", "", "Phase bins per period"),
            column("TBIN", "1D", "s", "", "Bin width"),
            column("CTR_FREQ", "1D", "MHz", "", "Band centre"),
            column("NCHAN", "1J", "", "", "Frequency channels"),
            column("CHAN_BW", "1D", "MHz", "", "Channel width"),
            column("REF_FREQ", "1D", "MHz", "", "Frequency of zero delay"),
            column("DM", "1D", "CM-3", "", "Dispersion measure"),
            column("RM", "1D", "RAD", "", "Rotation measure"),
            column("PR_CORR", "1I", "", "", "1 if parallactic angle corrected"),
            column("FD_CORR", "1I", "", "", "1 if feed corrected"),
            column("BE_CORR", "1I", "", "", "1 if backend corrected"),
            column("RM_CORR", "1I", "", "", "1 if Faraday rotation removed"),
            column("DEDISP", "1I", "", "", "1 if dedispersed"),
            column("DDS_MTHD", "32A", "", "", "How dispersion was removed"),
            column("SC_MTHD", "32A", "", "", "Scattered power method"),
            column("CAL_MTHD", "32A", "", "", "Polarisation calibration"),
            column("CAL_FILE", "256A", "", "", "Calibrator file"),
            column("RFI_MTHD", "32A", "", "", "RFI removal"),
            column("RM_MODEL", "32A", "", "", "Extra Faraday rotation model"),
            column("AUX_RM_C", "1I", "", "", "1 if extra rotation removed"),
            column("DM_MODEL", "32A", "", "", "Extra dispersion model"),
            column("AUX_DM_C", "1I", "", "", "1 if extra dispersion removed"),
        };
        out.write_header(
            bintable_header("HISTORY", columns, 1, {}, "Processing history"));
        std::vector<std::byte> row(
            static_cast<std::size_t>(row_width(columns)));
        std::byte* at   = row.data();
        const auto text = [&](std::size_t width, std::string_view value) {
            detail::put_text(at, width, value);
            at += width;
        };
        const auto i16 = [&](std::int16_t value) {
            store_be(at, value);
            at += 2;
        };
        const auto i32 = [&](std::int32_t value) {
            store_be(at, value);
            at += 4;
        };
        const auto f64 = [&](double value) {
            store_be(at, value);
            at += 8;
        };
        text(24, detail::local_ctime(spec.created));
        text(256,
             std::string_view(spec.command)
                 .substr(0, std::min<std::size_t>(spec.command.size(), 256U)));
        std::byte* const scale_cell = at;
        at += 8; // SCALE, filled in last (see below)
        text(8, spec.pol_type);
        i32(1);
        i16(static_cast<std::int16_t>(npol));
        // HISTORY NBIN is int16. Above 32767 bins, dspsr archives hold 32767
        // there and zeros in every later cell and in SCALE; we write the
        // same so such files compare equal.
        constexpr auto kInt16Max = std::numeric_limits<std::int16_t>::max();
        const bool overflow      = std::cmp_greater(nbin, kInt16Max);
        i16(overflow ? kInt16Max : static_cast<std::int16_t>(nbin));
        if (!overflow) {
            i16(static_cast<std::int16_t>(nbin));
            f64(spec.period / static_cast<double>(nbin));
            f64(spec.centre_frequency);
            i32(static_cast<std::int32_t>(nchan));
            f64(spec.bandwidth /
                static_cast<double>(static_cast<float>(nchan)));
            f64(0.0);
            f64(spec.dm);
            f64(0.0);
            i16(0);
            i16(0);
            i16(0);
            i16(0);
            i16(spec.dedispersed ? 1 : 0);
            text(32, "UNSET");
            text(32, "NONE");
            text(32, "NONE");
            text(256, "NONE");
            text(32, "NONE");
            text(32, "NONE");
            i16(0);
            text(32, "NONE");
            i16(0);
            detail::put_text(scale_cell, 8, "FluxDen");
        }
        out.write_data(row);
    }

    // SUBINT: one row.
    {
        const auto& sub = spec.input_subint;
        std::vector<Card> extra;
        extra.push_back(Card::string(
            "EPOCHS", "VALID", "OFFS_SUB marks the sub-integration epoch"));
        extra.push_back(
            Card::string("INT_TYPE", "TIME", "Rows are indexed by time"));
        extra.push_back(Card::string("INT_UNIT", "SEC", "Row index unit"));
        extra.push_back(Card::string("SCALE", "FluxDen", "Amplitude units"));
        extra.push_back(Card::string("POL_TYPE", spec.pol_type,
                                     "Stored polarisation products"));
        extra.push_back(Card::integer("NPOL", static_cast<std::int64_t>(npol),
                                      "Polarisation products"));
        const auto copy_real = [&](std::string_view key, std::string comment) {
            if (const auto value = sub.find_number(key)) {
                extra.push_back(
                    Card::real(std::string{key}, *value, std::move(comment)));
            } else {
                extra.push_back(
                    Card::string(std::string{key}, "*", std::move(comment)));
            }
        };
        const auto copy_int = [&](std::string_view key, std::string comment) {
            if (const auto value = sub.find_number(key)) {
                extra.push_back(Card::integer(std::string{key},
                                              std::llround(*value),
                                              std::move(comment)));
            } else {
                extra.push_back(
                    Card::string(std::string{key}, "*", std::move(comment)));
            }
        };
        copy_real("TBIN", "[s] Input sample interval");
        extra.push_back(Card::integer("NBIN", static_cast<std::int64_t>(nbin),
                                      "Phase bins"));
        extra.push_back(Card::string("NBIN_PRD", "*",
                                     "Phase bins per period (gated data)"));
        extra.push_back(
            Card::string("PHS_OFFS", "*", "Phase of bin 0 (gated data)"));
        copy_int("NBITS", "Input bits per sample");
        copy_real("ZERO_OFF", "Input zero offset");
        copy_int("SIGNINT", "1 if input samples are signed");
        extra.push_back(
            Card::string("NSUBOFFS", "*", "Row offset in a split observation"));
        extra.push_back(Card::integer("NCHAN", static_cast<std::int64_t>(nchan),
                                      "Frequency channels"));
        extra.push_back(Card::real("CHAN_BW",
                                   spec.bandwidth / static_cast<double>(nchan),
                                   "[MHz] Channel width"));
        extra.push_back(
            Card::real("REFFREQ", spec.centre_frequency, "[MHz] Band centre"));
        extra.push_back(
            Card::real("DM", spec.dm, "[cm-3 pc] Dispersion measure"));
        extra.push_back(Card::real("RM", 0.0, "[rad m-2] Rotation measure"));
        extra.push_back(
            Card::string("NCHNOFFS", "*", "Channel offset in a split band"));
        copy_int("NSBLK", "Input samples per row");
        extra.push_back(Card::string("NSTOT", "*", "Input samples in total"));
        extra.push_back(Card::integer("NINFO", 0, "Extra information columns"));
        if (sub.contains("NCH_STRT")) {
            copy_int("NCH_STRT", "First channel in this file");
        }

        const std::vector<TableColumn> columns{
            column("INDEXVAL", "1D", "", "", "Row index when not by time"),
            column("TSUBINT", "1D", "s", "", "Sub-integration length"),
            column("OFFS_SUB", "1D", "s", "", "Epoch minus start time"),
            column("PERIOD", "1D", "s", "", "Folding period"),
            column("AUX_DM", "1D", "CM-3", "", "Extra dispersion measure"),
            column("AUX_RM", "1D", "RAD", "", "Extra rotation measure"),
            column("DAT_FREQ", fmt::format("{}D", nchan), "MHz", "",
                   "Channel centre frequencies"),
            column("DAT_WTS", fmt::format("{}E", nchan), "", "",
                   "Channel weights"),
            column("DAT_OFFS", fmt::format("{}E", nprof), "", "",
                   "Profile offsets"),
            column("DAT_SCL", fmt::format("{}E", nprof), "", "",
                   "Profile scales: value = DATA * scale + offset"),
            column("DATA", fmt::format("{}I", nprof * nbin), "Jy",
                   fmt::format("({},{},{})", nbin, nchan, npol),
                   "Quantised profiles"),
        };
        out.write_header(
            bintable_header("SUBINT", columns, 1, extra, "Folded profiles"));

        // Scales and offsets of every profile; they precede DATA.
        std::vector<float> offsets(nprof);
        std::vector<float> scales(nprof);
        for (std::size_t prof = 0; prof < nprof; ++prof) {
            const auto scale = profile_scale(amps.subspan(prof * nbin, nbin));
            offsets[prof]    = scale.offset;
            scales[prof]     = scale.scale;
        }

        std::vector<std::byte> head(std::size_t{48} + (nchan * 8U) +
                                    (nchan * 4U) + (nprof * 8U));
        std::byte* at = head.data();
        // INDEXVAL, TSUBINT, OFFS_SUB, PERIOD, AUX_DM, AUX_RM.
        const std::array<double, 6> scalars{
            0.0, spec.tsubint, spec.offs_sub, spec.period, 0.0, 0.0,
        };
        for (const double value : scalars) {
            store_be(at, value);
            at += 8;
        }
        for (const double freq : spec.frequencies) {
            store_be(at, freq);
            at += 8;
        }
        for (std::size_t chan = 0; chan < nchan; ++chan) {
            store_be(at, 1.0F);
            at += 4;
        }
        for (const float value : offsets) {
            store_be(at, value);
            at += 4;
        }
        for (const float value : scales) {
            store_be(at, value);
            at += 4;
        }
        out.write_data(head);

        // DATA: quantise and write a few profiles at a time, so no second
        // copy of the archive is held.
        constexpr std::size_t kChunkValues = std::size_t{1} << 21U;
        const auto per_chunk = std::max<std::size_t>(1U, kChunkValues / nbin);
        std::vector<std::int16_t> quantised(std::min(nprof, per_chunk) * nbin);
        std::vector<std::byte> buffer(quantised.size() * 2U);
        for (std::size_t first = 0; first < nprof; first += per_chunk) {
            const auto count  = std::min(per_chunk, nprof - first);
            const auto values = count * nbin;
            for (std::size_t prof = first; prof < first + count; ++prof) {
                quantise_with(amps.subspan(prof * nbin, nbin),
                              ProfileScale{
                                  .scale  = scales[prof],
                                  .offset = offsets[prof],
                              },
                              std::span<std::int16_t>(quantised).subspan(
                                  (prof - first) * nbin, nbin));
            }
            for (std::size_t index = 0; index < values; ++index) {
                store_be(buffer.data() + (index * 2U), quantised[index]);
            }
            out.write_data(
                std::span<const std::byte>(buffer).first(values * 2U));
        }
    }
    out.close();
}

} // namespace psrio::formats::psrfits
