#pragma once

#include "psrio/astro.hpp"
#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <format>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace psrio::formats::sigproc {

/// A header key this library does not interpret.
/// @p raw is the value bytes only, not the length-prefixed name.
struct ExtraKey {
    std::string name;
    std::vector<std::byte> raw;
};

/// Count of time samples for FilterbankReader::read.
struct SampleCount {
    std::uint64_t value{0};
};

/// Count of raw payload bytes for FilterbankReader::read and view.
/// Must be a positive multiple of FilterbankHeader::bytes_per_sample().
struct ByteCount {
    std::uint64_t value{0};
};

struct NameId {
    std::int32_t id{0};
    std::string_view name;
};

/// Telescope ids recognized in pulsar and FRB archives.
inline constexpr auto kTelescopes = std::to_array<NameId>({
    {.id = 0, .name = "Fake"},       {.id = 1, .name = "Arecibo"},
    {.id = 2, .name = "Ooty"},       {.id = 3, .name = "Nancay"},
    {.id = 4, .name = "Parkes"},     {.id = 5, .name = "Jodrell"},
    {.id = 6, .name = "GBT"},        {.id = 7, .name = "GMRT"},
    {.id = 8, .name = "Effelsberg"}, {.id = 9, .name = "Effelsberg LOFAR"},
    {.id = 10, .name = "SRT"},       {.id = 11, .name = "LOFAR"},
    {.id = 12, .name = "VLA"},       {.id = 20, .name = "CHIME"},
    {.id = 21, .name = "FAST"},      {.id = 30, .name = "MWA"},
    {.id = 40, .name = "NC"},        {.id = 41, .name = "NGNC"},
    {.id = 64, .name = "MeerKAT"},   {.id = 65, .name = "KAT-7"},
    {.id = 82, .name = "eMerlin"},   {.id = 1916, .name = "I-LOFAR"},
});

/// Backend ids. Id 10 is "BPSR"; DSPSR uses "ARTEMIS" for
/// that id when the telescope is not Parkes.
inline constexpr auto kMachines = std::to_array<NameId>({
    {.id = 0, .name = "FAKE"},       {.id = 1, .name = "PSPM"},
    {.id = 2, .name = "WAPP"},       {.id = 3, .name = "AOFTM"},
    {.id = 4, .name = "BPP"},        {.id = 5, .name = "OOTY"},
    {.id = 6, .name = "SCAMP"},      {.id = 7, .name = "GMRTFB"},
    {.id = 8, .name = "PULSAR2000"}, {.id = 9, .name = "PARSPEC"},
    {.id = 10, .name = "BPSR"},      {.id = 11, .name = "COBALT"},
    {.id = 14, .name = "GMRTNEW"},   {.id = 20, .name = "CHIME"},
    {.id = 30, .name = "MWA-VCS"},   {.id = 31, .name = "MWAX-VCS"},
    {.id = 32, .name = "MWAX-RTB"},  {.id = 40, .name = "ADU"},
    {.id = 41, .name = "iTPM"},      {.id = 83, .name = "ROACH"},
});

[[nodiscard]] inline std::string_view name_from(std::span<const NameId> table,
                                                std::int32_t id) noexcept {
    for (const NameId& entry : table) {
        if (entry.id == id) {
            return entry.name;
        }
    }
    return "unknown";
}

[[nodiscard]] inline std::string_view telescope_name(std::int32_t id) noexcept {
    return name_from(kTelescopes, id);
}

[[nodiscard]] inline std::string_view machine_name(std::int32_t id) noexcept {
    return name_from(kMachines, id);
}

inline constexpr auto kDataTypes = std::to_array<NameId>({
    {.id = 0, .name = "raw data"},
    {.id = 1, .name = "filterbank"},
    {.id = 2, .name = "time series"},
    {.id = 3, .name = "pulse profiles"},
    {.id = 4, .name = "amplitude spectrum"},
    {.id = 5, .name = "complex spectrum"},
    {.id = 6, .name = "dedispersed subbands"},
});

/// SIGPROC data_type names. Ids outside 0–6 are "unknown".
[[nodiscard]] inline std::string_view data_type_name(std::int32_t id) noexcept {
    return name_from(kDataTypes, id);
}

/// Parsed SIGPROC filterbank header.
///
/// On disk the header is little-endian, from the length-prefixed string
/// HEADER_START through HEADER_END. Integer keys are int32. The key "signed"
/// is one int8. Floating keys are IEEE-754 binary64. A missing optional key
/// stays disengaged, which is distinct from a stored zero. nifs defaults to
/// 1 when its key is absent.
///
/// nsamples() is the number of complete time samples a reader can return.
/// A missing nsamples key, or a stored 0, means "use the file length". Any
/// other stored value that disagrees with the file is reported in warnings()
/// and the smaller count is used. The library does not write that warning
/// to stderr.
struct FilterbankHeader {
    /// Number of frequency channels. Required to be positive by validate().
    std::int32_t nchans{0};

    /// Number of IFs. Defaults to 1 when the key is absent.
    std::int32_t nifs{1};

    /// Bits per sample. Supported values are 1, 2, 4, 8, 16, and 32.
    std::int32_t nbits{0};

    /// Sample interval, in seconds.
    std::optional<double> tsamp;

    /// Modified Julian Date of the first sample.
    std::optional<double> tstart;

    /// Centre frequency of the first channel, in MHz.
    std::optional<double> fch1;

    /// Signed channel spacing, in MHz. Negative when the band steps downward.
    std::optional<double> foff;

    /// Explicit channel centres, in MHz, from fchannel entries.
    /// Empty when the file uses fch1 and foff. channel_frequencies() prefers
    /// this table when it is non-empty.
    std::vector<double> frequency_table;

    /// SIGPROC data_type. 0 is raw, 1 is filterbank.
    std::optional<std::int32_t> data_type;

    /// SIGPROC telescope_id. 0 is a real id ("Fake"), so absence is disengaged.
    std::optional<std::int32_t> telescope_id;

    /// SIGPROC machine_id.
    std::optional<std::int32_t> machine_id;

    /// Source name. The following on-disk token is a string, not a keyword.
    std::optional<std::string> source_name;

    /// Original raw file name. The following on-disk token is a string.
    std::optional<std::string> rawdatafile;

    /// Right ascension packed as hhmmss.ss (J2000).
    std::optional<double> src_raj;

    /// Declination packed as ddmmss.ss (J2000).
    std::optional<double> src_dej;

    /// Azimuth at the start of the scan, in degrees.
    std::optional<double> az_start;

    /// Zenith angle at the start of the scan, in degrees.
    std::optional<double> za_start;

    /// Non-zero when the data are barycentric. Stored as the on-disk int32.
    std::optional<std::int32_t> barycentric;

    /// Non-zero when the data are pulsarcentric. Stored as the on-disk int32.
    std::optional<std::int32_t> pulsarcentric;

    /// Beam index.
    std::optional<std::int32_t> ibeam;

    /// Number of beams.
    std::optional<std::int32_t> nbeams;

    /// Reference dispersion measure, in pc cm^-3.
    std::optional<double> refdm;

    /// Folding period, in seconds.
    std::optional<double> period;

    /// Legacy profile bin count. Not used while reading filterbank data.
    std::optional<std::int32_t> nbins;

    /// Legacy pulse count. Not used while reading filterbank data.
    std::optional<std::int32_t> npuls;

    /// On-disk "signed" byte. Any non-zero value means 8-bit samples are int8.
    std::optional<std::int8_t> signed_data;

    /// Legacy nsamples key, as a uint32. Zero and absence both mean unset.
    std::optional<std::uint32_t> declared_nsamples;

    /// Bytes from the start of the file through HEADER_END.
    std::uint64_t header_bytes{0};

    /// Bytes after the header, including any trailing partial sample.
    std::uint64_t data_bytes{0};

    /// header_bytes + data_bytes.
    std::uint64_t file_bytes{0};

    /// Keys outside the SIGPROC schema, in file order.
    std::vector<ExtraKey> extra;

    /// Parse a SIGPROC header from the leading bytes of a file image.
    /// @throws FormatError if the image is empty, truncated, or not SIGPROC.
    [[nodiscard]] static FilterbankHeader
    parse(std::span<const std::byte> file);

    /// Parse as parse(file), including @p path in error messages.
    [[nodiscard]] static FilterbankHeader parse(std::span<const std::byte> file,
                                                std::string_view path);

    /// Check that this header describes readable filterbank data.
    /// A declared nsamples that disagrees with the file is not an error.
    /// @throws ValidationError if geometry, frequency, tsamp, or data_type
    ///         cannot be read as filterbank data.
    void validate() const;

    /// Bytes in one time sample: nifs * nchans * nbits / 8.
    /// @throws ValidationError if that product is not a positive whole number
    /// of bytes.
    [[nodiscard]] std::uint64_t bytes_per_sample() const;

    /// Complete time samples that fit in data_bytes, ignoring
    /// declared_nsamples.
    [[nodiscard]] std::uint64_t samples_in_file() const;

    /// Time samples a reader will return. See the class note.
    [[nodiscard]] std::uint64_t nsamples() const;

    /// Bytes after the last readable sample. Includes a partial trailing
    /// sample and any complete samples past a smaller declared count.
    [[nodiscard]] std::uint64_t trailing_bytes() const;

    /// Absolute bandwidth |foff| * nchans, in MHz.
    /// @throws ValidationError if foff is absent or nchans is not positive.
    [[nodiscard]] double bandwidth() const;

    /// Upper edge of the first channel, fch1 - 0.5 * foff, in MHz.
    [[nodiscard]] double ftop() const;

    /// Lower edge of the last channel, in MHz.
    [[nodiscard]] double fbottom() const;

    /// Centre frequency of the uniform grid, in MHz.
    [[nodiscard]] double fcenter() const;

    /// Write one centre frequency per channel, in MHz, into @p out.
    /// Uses frequency_table when it is non-empty, otherwise fch1 + i * foff.
    void channel_frequencies(std::span<double> out) const;

    /// Name for telescope_id, or "unknown" when the id is absent or unlisted.
    [[nodiscard]] std::string_view telescope_name() const noexcept;

    /// Name for machine_id, or "unknown" when the id is absent or unlisted.
    [[nodiscard]] std::string_view machine_name() const noexcept;

    /// Name for data_type, or "unknown" when the id is absent or unlisted.
    [[nodiscard]] std::string_view data_type_name() const noexcept;

    /// True when the "signed" key is present and non-zero.
    [[nodiscard]] bool samples_are_signed() const noexcept {
        return signed_data.has_value() && *signed_data != 0;
    }

    /// Right Ascension formatted as "HH:MM:SS.ssss". Empty string if absent.
    [[nodiscard]] std::string ra_string() const;

    /// Declination formatted as "[+/-]DD:MM:SS.ssss". Empty string if absent.
    [[nodiscard]] std::string dec_string() const;

    /// Right Ascension in decimal hours [0, 24).
    [[nodiscard]] std::optional<double> ra_hours() const noexcept;

    /// Right Ascension in decimal degrees [0, 360).
    [[nodiscard]] std::optional<double> ra_degrees() const noexcept;

    /// Right Ascension in radians [0, 2*pi).
    [[nodiscard]] std::optional<double> ra_radians() const noexcept;

    /// Declination in decimal degrees [-90, +90].
    [[nodiscard]] std::optional<double> dec_degrees() const noexcept;

    /// Declination in radians [-pi/2, +pi/2].
    [[nodiscard]] std::optional<double> dec_radians() const noexcept;

    /// Total observation duration in seconds: nsamples() * tsamp.
    [[nodiscard]] double observation_duration() const;

    /// Human-readable duration string (e.g. "120.0 seconds", "2.5 hours").
    [[nodiscard]] std::string duration_string() const;

    /// Calendar date of the observation start as "YYYY-MM-DD".
    [[nodiscard]] std::string gregorian_date() const;

    /// MJD epoch at sample offset @p sample_index.
    [[nodiscard]] double mjd_after_samples(std::uint64_t sample_index) const;

    /// Calculate dispersion delays in seconds for every frequency channel.
    /// @param dm Dispersion measure in pc cm^-3.
    /// @param out Destination span, must have size equal to nchans.
    /// @param ref_freq_mhz Reference frequency in MHz (0.0 defaults to highest
    /// frequency channel).
    void dispersion_delays(double dm,
                           std::span<double> out,
                           double ref_freq_mhz = 0.0) const;

    /// Warnings recorded by parse(). Currently only an nsamples disagreement.
    [[nodiscard]] std::span<const std::string> warnings() const noexcept {
        return m_warnings;
    }

private:
    [[nodiscard]] std::optional<std::uint64_t>
    try_bytes_per_sample() const noexcept;

    void note_sample_count();

    std::vector<std::string> m_warnings;
};

namespace detail {

inline constexpr auto kSchemaTokens = std::to_array<std::string_view>({
    "HEADER_START",  "HEADER_END",  "FREQUENCY_START",
    "FREQUENCY_END", "rawdatafile", "source_name",
    "az_start",      "za_start",    "src_raj",
    "src_dej",       "tstart",      "tsamp",
    "period",        "fch1",        "foff",
    "fchannel",      "nchans",      "telescope_id",
    "machine_id",    "data_type",   "ibeam",
    "nbeams",        "nbits",       "barycentric",
    "pulsarcentric", "nbins",       "nsamples",
    "nifs",          "npuls",       "refdm",
    "signed",
});

[[nodiscard]] inline bool is_schema_token(std::string_view token) noexcept {
    return std::ranges::any_of(
        kSchemaTokens,
        [token](const std::string_view entry) { return entry == token; });
}

[[nodiscard]] inline std::string
error_at(std::string_view path, std::uint64_t offset, std::string_view detail) {
    std::string message = "psrio: ";
    if (!path.empty()) {
        message.append(path);
        message.append(": ");
    }
    message.append("at byte ");
    message.append(std::to_string(offset));
    message.append(": ");
    message.append(detail);
    return message;
}

[[nodiscard]] inline bool
starts_with_schema_token(std::span<const std::byte> file,
                         std::size_t position) {
    if (position > file.size() || file.size() - position < 4U) {
        return false;
    }
    const auto length = ::psrio::detail::load_little_endian<std::uint32_t>(
        file.data() + position);
    if (length < 1U || length > 4096U) {
        return false;
    }
    if (file.size() - position - 4U < length) {
        return false;
    }
    const auto* const text =
        reinterpret_cast<const char*>(file.data() + position + 4U);
    return is_schema_token(std::string_view{text, length});
}

class ByteCursor {
public:
    ByteCursor(std::span<const std::byte> file, std::string_view path)
        : m_file(file),
          m_path(path) {}

    [[nodiscard]] std::uint64_t position() const noexcept { return m_pos; }

    [[nodiscard]] std::string_view path() const noexcept { return m_path; }

    [[nodiscard]] std::span<const std::byte> file() const noexcept {
        return m_file;
    }

    void seek(std::size_t position) noexcept { m_pos = position; }

    [[nodiscard]] std::string read_string(bool key_token = false) {
        const auto start = m_pos;
        if (m_pos > m_file.size() || m_file.size() - m_pos < 4U) {
            throw FormatError(error_at(m_path, start,
                                       key_token
                                           ? "header ended before HEADER_END"
                                           : "truncated string length"));
        }
        const auto length = ::psrio::detail::load_little_endian<std::uint32_t>(
            m_file.data() + m_pos);
        m_pos += 4U;
        if (length < 1U || length > 4096U) {
            if (start == 0U && length == 0x0C000000U) {
                throw FormatError(error_at(
                    m_path, 0, "big-endian SIGPROC header is not supported"));
            }
            throw FormatError(error_at(m_path, start,
                                       "implausible string length " +
                                           std::to_string(length)));
        }
        if (m_pos > m_file.size() || m_file.size() - m_pos < length) {
            throw FormatError(error_at(m_path, start,
                                       key_token ? "truncated header key"
                                                 : "truncated string"));
        }
        std::string value(length, '\0');
        std::memcpy(value.data(), m_file.data() + m_pos, length);
        m_pos += length;
        return value;
    }

    template <typename T> [[nodiscard]] T read_number(std::string_view key) {
        if (m_pos > m_file.size() || m_file.size() - m_pos < sizeof(T)) {
            throw FormatError(
                error_at(m_path, m_pos,
                         std::string("truncated value for key ").append(key)));
        }
        const T value =
            ::psrio::detail::load_little_endian<T>(m_file.data() + m_pos);
        m_pos += sizeof(T);
        return value;
    }

private:
    std::span<const std::byte> m_file;
    std::string_view m_path;
    std::size_t m_pos{0};
};

inline bool consume_known(ByteCursor& cursor,
                          FilterbankHeader& header,
                          std::string_view key) {
    if (key == "nchans") {
        header.nchans = cursor.read_number<std::int32_t>(key);
        return true;
    }
    if (key == "nifs") {
        header.nifs = cursor.read_number<std::int32_t>(key);
        return true;
    }
    if (key == "nbits") {
        header.nbits = cursor.read_number<std::int32_t>(key);
        return true;
    }
    if (key == "telescope_id") {
        header.telescope_id = cursor.read_number<std::int32_t>(key);
        return true;
    }
    if (key == "machine_id") {
        header.machine_id = cursor.read_number<std::int32_t>(key);
        return true;
    }
    if (key == "data_type") {
        header.data_type = cursor.read_number<std::int32_t>(key);
        return true;
    }
    if (key == "barycentric") {
        header.barycentric = cursor.read_number<std::int32_t>(key);
        return true;
    }
    if (key == "pulsarcentric") {
        header.pulsarcentric = cursor.read_number<std::int32_t>(key);
        return true;
    }
    if (key == "ibeam") {
        header.ibeam = cursor.read_number<std::int32_t>(key);
        return true;
    }
    if (key == "nbeams") {
        header.nbeams = cursor.read_number<std::int32_t>(key);
        return true;
    }
    if (key == "nbins") {
        header.nbins = cursor.read_number<std::int32_t>(key);
        return true;
    }
    if (key == "npuls") {
        header.npuls = cursor.read_number<std::int32_t>(key);
        return true;
    }
    if (key == "nsamples") {
        header.declared_nsamples = cursor.read_number<std::uint32_t>(key);
        return true;
    }
    if (key == "az_start") {
        header.az_start = cursor.read_number<double>(key);
        return true;
    }
    if (key == "za_start") {
        header.za_start = cursor.read_number<double>(key);
        return true;
    }
    if (key == "src_raj") {
        header.src_raj = cursor.read_number<double>(key);
        return true;
    }
    if (key == "src_dej") {
        header.src_dej = cursor.read_number<double>(key);
        return true;
    }
    if (key == "tstart") {
        header.tstart = cursor.read_number<double>(key);
        return true;
    }
    if (key == "tsamp") {
        header.tsamp = cursor.read_number<double>(key);
        return true;
    }
    if (key == "fch1") {
        header.fch1 = cursor.read_number<double>(key);
        return true;
    }
    if (key == "foff") {
        header.foff = cursor.read_number<double>(key);
        return true;
    }
    if (key == "refdm") {
        header.refdm = cursor.read_number<double>(key);
        return true;
    }
    if (key == "period") {
        header.period = cursor.read_number<double>(key);
        return true;
    }
    if (key == "fchannel") {
        header.frequency_table.push_back(cursor.read_number<double>(key));
        return true;
    }
    if (key == "source_name") {
        header.source_name = cursor.read_string();
        return true;
    }
    if (key == "rawdatafile") {
        header.rawdatafile = cursor.read_string();
        return true;
    }
    if (key == "signed") {
        header.signed_data = cursor.read_number<std::int8_t>(key);
        return true;
    }
    return false;
}

inline void skip_unknown(ByteCursor& cursor,
                         FilterbankHeader& header,
                         std::string_view key) {
    const auto origin = static_cast<std::size_t>(cursor.position());
    const auto file   = cursor.file();
    const std::string_view name(key);

    // 1. Probe numeric candidate widths: 4-byte (int/uint), 8-byte (double),
    // 1-byte (char)
    for (const std::size_t width :
         {std::size_t{4}, std::size_t{8}, std::size_t{1}}) {
        if (origin <= file.size() && file.size() - origin >= width) {
            const auto next = origin + width;
            if (starts_with_schema_token(file, next)) {
                ExtraKey extra;
                extra.name = name;
                extra.raw.assign(file.data() + origin, file.data() + next);
                header.extra.push_back(std::move(extra));
                cursor.seek(next);
                return;
            }
        }
    }

    // 2. Probe length-prefixed string candidate: 4-byte length + length bytes
    if (origin <= file.size() && file.size() - origin >= 4U) {
        const auto str_len = ::psrio::detail::load_little_endian<std::uint32_t>(
            file.data() + origin);
        if (str_len >= 1U && str_len <= 4096U) {
            const auto total_width = 4U + static_cast<std::size_t>(str_len);
            if (file.size() - origin >= total_width) {
                const auto next = origin + total_width;
                if (starts_with_schema_token(file, next)) {
                    ExtraKey extra;
                    extra.name = name;
                    extra.raw.assign(file.data() + origin, file.data() + next);
                    header.extra.push_back(std::move(extra));
                    cursor.seek(next);
                    return;
                }
            }
        }
    }

    throw FormatError(error_at(
        cursor.path(), origin,
        std::format("cannot determine the size of unknown header key \"{}\"",
                    name)));
}

} // namespace detail

inline FilterbankHeader
FilterbankHeader::parse(std::span<const std::byte> file) {
    return parse(file, {});
}

inline FilterbankHeader FilterbankHeader::parse(std::span<const std::byte> file,
                                                std::string_view path) {
    if (file.empty()) {
        throw FormatError(detail::error_at(path, 0, "file is empty"));
    }

    detail::ByteCursor cursor(file, path);
    const std::string magic = cursor.read_string();
    if (magic != "HEADER_START") {
        throw FormatError(detail::error_at(path, 0, "missing HEADER_START"));
    }

    FilterbankHeader header;
    while (true) {
        const auto key_at     = cursor.position();
        const std::string key = cursor.read_string(true);
        if (key == "HEADER_END") {
            break;
        }
        if (key == "HEADER_START") {
            throw FormatError(
                detail::error_at(path, key_at, "nested HEADER_START"));
        }
        if (key == "FREQUENCY_START" || key == "FREQUENCY_END") {
            continue;
        }
        if (!detail::consume_known(cursor, header, key)) {
            detail::skip_unknown(cursor, header, key);
        }
    }

    header.header_bytes = cursor.position();
    header.file_bytes   = file.size();
    if (header.header_bytes > header.file_bytes) {
        throw FormatError(detail::error_at(path, header.header_bytes,
                                           "header extends past end of file"));
    }
    header.data_bytes = header.file_bytes - header.header_bytes;
    header.note_sample_count();
    return header;
}

inline std::optional<std::uint64_t>
FilterbankHeader::try_bytes_per_sample() const noexcept {
    if (nchans <= 0 || nifs <= 0) {
        return std::nullopt;
    }
    if (nbits != 1 && nbits != 2 && nbits != 4 && nbits != 8 && nbits != 16 &&
        nbits != 32) {
        return std::nullopt;
    }
    const auto channels = static_cast<std::uint64_t>(nchans);
    const auto ifs      = static_cast<std::uint64_t>(nifs);
    const auto width    = static_cast<std::uint64_t>(nbits);
    if (channels > (std::numeric_limits<std::uint64_t>::max() / ifs)) {
        return std::nullopt;
    }
    const auto per_sample_channels = channels * ifs;
    if (per_sample_channels >
        (std::numeric_limits<std::uint64_t>::max() / width)) {
        return std::nullopt;
    }
    const auto bits = per_sample_channels * width;
    if (bits == 0U || bits % 8U != 0U) {
        return std::nullopt;
    }
    return bits / 8U;
}

inline std::uint64_t FilterbankHeader::bytes_per_sample() const {
    if (nchans <= 0) {
        throw ValidationError("psrio: nchans must be positive");
    }
    if (nifs <= 0) {
        throw ValidationError("psrio: nifs must be positive");
    }
    if (nbits != 1 && nbits != 2 && nbits != 4 && nbits != 8 && nbits != 16 &&
        nbits != 32) {
        throw ValidationError(
            std::format("psrio: nbits {} is not supported", nbits));
    }
    if (const auto stride = try_bytes_per_sample()) {
        return *stride;
    }
    throw ValidationError(
        "psrio: nifs * nchans * nbits does not form a byte stride");
}

inline void FilterbankHeader::note_sample_count() {
    if (!declared_nsamples.has_value() || *declared_nsamples == 0U) {
        return;
    }
    const auto stride = try_bytes_per_sample();
    if (!stride.has_value()) {
        return;
    }
    const auto in_file = data_bytes / *stride;
    if (in_file == *declared_nsamples) {
        return;
    }
    const auto used = std::min<std::uint64_t>(*declared_nsamples, in_file);
    m_warnings.push_back(std::format("header nsamples {} disagrees with {} "
                                     "complete samples in the file; using {}",
                                     *declared_nsamples, in_file, used));
}

inline void FilterbankHeader::validate() const {
    (void)bytes_per_sample();
    if (!frequency_table.empty()) {
        if (frequency_table.size() != static_cast<std::size_t>(nchans)) {
            throw ValidationError(
                "psrio: frequency table length does not equal nchans");
        }
    } else if (!fch1.has_value() || !foff.has_value()) {
        throw ValidationError(
            "psrio: header needs fch1 and foff, or a frequency table");
    }
    if (!tsamp.has_value() || !(*tsamp > 0.0)) {
        throw ValidationError("psrio: tsamp must be present and positive");
    }
    if (data_type.has_value() && *data_type != 0 && *data_type != 1) {
        throw ValidationError(std::format(
            "psrio: data_type {} ({}) is not filterbank data; other SIGPROC "
            "types are read by a separate reader",
            *data_type, sigproc::data_type_name(*data_type)));
    }
}

inline std::uint64_t FilterbankHeader::samples_in_file() const {
    return data_bytes / bytes_per_sample();
}

inline std::uint64_t FilterbankHeader::nsamples() const {
    const auto in_file = samples_in_file();
    if (!declared_nsamples.has_value() || *declared_nsamples == 0U) {
        return in_file;
    }
    return std::min<std::uint64_t>(*declared_nsamples, in_file);
}

inline std::uint64_t FilterbankHeader::trailing_bytes() const {
    const auto used = nsamples() * bytes_per_sample();
    if (used > data_bytes) {
        return 0;
    }
    return data_bytes - used;
}

inline double FilterbankHeader::bandwidth() const {
    if (nchans <= 0 || !foff.has_value()) {
        throw ValidationError("psrio: bandwidth needs nchans and foff");
    }
    return std::fabs(*foff) * static_cast<double>(nchans);
}

inline double FilterbankHeader::ftop() const {
    if (!fch1.has_value() || !foff.has_value()) {
        throw ValidationError("psrio: ftop needs fch1 and foff");
    }
    return *fch1 - (0.5 * (*foff));
}

inline double FilterbankHeader::fbottom() const {
    if (nchans <= 0 || !fch1.has_value() || !foff.has_value()) {
        throw ValidationError("psrio: fbottom needs fch1, foff, and nchans");
    }
    return ftop() + ((*foff) * static_cast<double>(nchans));
}

inline double FilterbankHeader::fcenter() const {
    if (nchans <= 0 || !fch1.has_value() || !foff.has_value()) {
        throw ValidationError("psrio: fcenter needs fch1, foff, and nchans");
    }
    return ftop() + (0.5 * (*foff) * static_cast<double>(nchans));
}

inline void FilterbankHeader::channel_frequencies(std::span<double> out) const {
    if (nchans <= 0 || out.size() != static_cast<std::size_t>(nchans)) {
        throw ValidationError(
            "psrio: channel frequency span must have one entry per channel");
    }
    if (!frequency_table.empty()) {
        if (frequency_table.size() != out.size()) {
            throw ValidationError(
                "psrio: frequency table length does not equal nchans");
        }
        std::ranges::copy(frequency_table, out.begin());
        return;
    }
    if (!fch1.has_value() || !foff.has_value()) {
        throw ValidationError(
            "psrio: fch1 and foff are required to derive channel frequencies");
    }
    for (std::size_t index = 0; index < out.size(); ++index) {
        out[index] = *fch1 + (static_cast<double>(index) * (*foff));
    }
}

inline std::string_view FilterbankHeader::telescope_name() const noexcept {
    if (!telescope_id.has_value()) {
        return "unknown";
    }
    return sigproc::telescope_name(*telescope_id);
}

inline std::string_view FilterbankHeader::machine_name() const noexcept {
    if (!machine_id.has_value()) {
        return "unknown";
    }
    return sigproc::machine_name(*machine_id);
}

inline std::string_view FilterbankHeader::data_type_name() const noexcept {
    if (!data_type.has_value()) {
        return "unknown";
    }
    return sigproc::data_type_name(*data_type);
}

inline std::string FilterbankHeader::ra_string() const {
    if (!src_raj.has_value()) {
        return {};
    }
    return astro::ra_to_string(*src_raj);
}

inline std::string FilterbankHeader::dec_string() const {
    if (!src_dej.has_value()) {
        return {};
    }
    return astro::dec_to_string(*src_dej);
}

inline std::optional<double> FilterbankHeader::ra_hours() const noexcept {
    if (!src_raj.has_value()) {
        return std::nullopt;
    }
    return astro::ra_to_hours(*src_raj);
}

inline std::optional<double> FilterbankHeader::ra_degrees() const noexcept {
    if (!src_raj.has_value()) {
        return std::nullopt;
    }
    return astro::ra_to_degrees(*src_raj);
}

inline std::optional<double> FilterbankHeader::ra_radians() const noexcept {
    if (!src_raj.has_value()) {
        return std::nullopt;
    }
    return astro::ra_to_radians(*src_raj);
}

inline std::optional<double> FilterbankHeader::dec_degrees() const noexcept {
    if (!src_dej.has_value()) {
        return std::nullopt;
    }
    return astro::dec_to_degrees(*src_dej);
}

inline std::optional<double> FilterbankHeader::dec_radians() const noexcept {
    if (!src_dej.has_value()) {
        return std::nullopt;
    }
    return astro::dec_to_radians(*src_dej);
}

inline double FilterbankHeader::observation_duration() const {
    return static_cast<double>(nsamples()) * tsamp.value_or(0.0);
}

inline std::string FilterbankHeader::duration_string() const {
    return astro::format_duration(observation_duration());
}

inline std::string FilterbankHeader::gregorian_date() const {
    if (!tstart.has_value()) {
        return {};
    }
    return astro::mjd_to_gregorian(*tstart);
}

inline double
FilterbankHeader::mjd_after_samples(std::uint64_t sample_index) const {
    return tstart.value_or(0.0) +
           ((static_cast<double>(sample_index) * tsamp.value_or(0.0)) /
            86400.0);
}

inline void FilterbankHeader::dispersion_delays(double dm,
                                                std::span<double> out,
                                                double ref_freq_mhz) const {
    if (nchans <= 0 || out.size() != static_cast<std::size_t>(nchans)) {
        throw ValidationError(
            "psrio: channel delay span must have one entry per channel");
    }
    channel_frequencies(out);
    double ref_freq = ref_freq_mhz;
    if (ref_freq <= 0.0) {
        ref_freq = *std::ranges::max_element(out);
    }
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = astro::dispersion_delay(dm, out[i], ref_freq);
    }
}

} // namespace psrio::formats::sigproc
