#pragma once

#include "psrio/baseband.hpp"
#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <exception>
#include <format>
#include <map>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace psrio::formats::guppi {

inline constexpr std::size_t kRecordBytes        = 80;
inline constexpr std::size_t kMaxHeaderRecords   = 2560;
inline constexpr std::size_t kDirectioAlignBytes = 512;

using HeaderValue = std::variant<std::int64_t, double, std::string>;

/// One GUPPI RAW ASCII header: 80-byte `KEY = VALUE` records through `END`.
///
/// Keys are stored upper-case and matched without regard to case. This is a
/// baseband header. It is not a `psrio::Header` and it is not an intensity
/// block.
class GuppiHeader {
public:
    void set(std::string_view key, HeaderValue value) {
        m_entries[upper(key)] = std::move(value);
    }

    template <typename T> [[nodiscard]] T get(std::string_view key) const {
        const auto it = m_entries.find(upper(key));
        if (it == m_entries.end()) {
            throw FormatError(
                std::format("psrio: GUPPI header missing {}", upper(key)));
        }
        if (!std::holds_alternative<T>(it->second)) {
            throw FormatError(std::format(
                "psrio: GUPPI key {} does not hold the requested type",
                upper(key)));
        }
        return std::get<T>(it->second);
    }

    template <typename T>
    [[nodiscard]] T get(std::string_view key, const T& fallback) const {
        const auto it = m_entries.find(upper(key));
        if (it == m_entries.end()) {
            return fallback;
        }
        if (!std::holds_alternative<T>(it->second)) {
            throw FormatError(std::format(
                "psrio: GUPPI key {} does not hold the requested type",
                upper(key)));
        }
        return std::get<T>(it->second);
    }

    [[nodiscard]] bool has_key(std::string_view key) const {
        return m_entries.contains(upper(key));
    }

    void clear() { m_entries.clear(); }

    [[nodiscard]] std::size_t size() const noexcept { return m_entries.size(); }

    [[nodiscard]] const std::map<std::string, HeaderValue>& entries() const noexcept {
        return m_entries;
    }

    /// Complex polarizations stored in the block. Header `NPOL` of 4 counts as 2.
    [[nodiscard]] std::int64_t complex_npol() const {
        const auto npol = int_or("NPOL", std::int64_t{1});
        return npol < 2 ? std::int64_t{1} : std::int64_t{2};
    }

    [[nodiscard]] std::int64_t int_or(std::string_view key,
                                      std::int64_t fallback) const {
        const auto it = m_entries.find(upper(key));
        if (it == m_entries.end()) {
            return fallback;
        }
        if (std::holds_alternative<std::int64_t>(it->second)) {
            return std::get<std::int64_t>(it->second);
        }
        if (std::holds_alternative<double>(it->second)) {
            const auto value = std::get<double>(it->second);
            if (!std::isfinite(value) || std::trunc(value) != value) {
                throw FormatError(std::format(
                    "psrio: GUPPI key {} is not an integer", upper(key)));
            }
            return static_cast<std::int64_t>(value);
        }
        throw FormatError(std::format(
            "psrio: GUPPI key {} does not hold the requested type", upper(key)));
    }

    [[nodiscard]] double real_or(std::string_view key, double fallback) const {
        const auto it = m_entries.find(upper(key));
        if (it == m_entries.end()) {
            return fallback;
        }
        if (std::holds_alternative<double>(it->second)) {
            return std::get<double>(it->second);
        }
        if (std::holds_alternative<std::int64_t>(it->second)) {
            return static_cast<double>(std::get<std::int64_t>(it->second));
        }
        throw FormatError(std::format(
            "psrio: GUPPI key {} does not hold the requested type", upper(key)));
    }

    [[nodiscard]] std::string string_or(std::string_view key,
                                        std::string fallback) const {
        const auto it = m_entries.find(upper(key));
        if (it == m_entries.end()) {
            return fallback;
        }
        if (!std::holds_alternative<std::string>(it->second)) {
            throw FormatError(std::format(
                "psrio: GUPPI key {} does not hold the requested type",
                upper(key)));
        }
        return std::get<std::string>(it->second);
    }

    /// Parse one header from @p bytes. Returns the number of bytes consumed,
    /// including DIRECTIO padding when that key is a non-zero integer.
    std::uint64_t parse(std::span<const std::byte> bytes) {
        m_entries.clear();
        std::uint64_t pos = 0;
        bool ended        = false;
        for (std::size_t record = 0; record < kMaxHeaderRecords; ++record) {
            if (pos + kRecordBytes > bytes.size()) {
                throw FormatError("psrio: GUPPI header ended before END");
            }
            std::string text(kRecordBytes, '\0');
            for (std::size_t index = 0; index < kRecordBytes; ++index) {
                text[index] = static_cast<char>(std::to_integer<unsigned char>(
                    bytes[static_cast<std::size_t>(pos) + index]));
            }
            pos += kRecordBytes;
            if (text.starts_with("END")) {
                ended = true;
                break;
            }
            const auto eq = text.find('=');
            if (eq == std::string::npos) {
                continue;
            }
            add_keyword(upper(trim(text.substr(0, eq))),
                        trim(text.substr(eq + 1)));
        }
        if (!ended) {
            throw FormatError(
                "psrio: GUPPI header exceeded 2560 records without END");
        }
        if (has_key("DIRECTIO") && get<std::int64_t>("DIRECTIO") != 0) {
            const auto remainder = pos % kDirectioAlignBytes;
            const auto padding   = remainder == 0U
                                       ? std::uint64_t{0}
                                       : kDirectioAlignBytes - remainder;
            if (pos + padding > bytes.size()) {
                throw FormatError(
                    "psrio: GUPPI DIRECTIO padding extends past the buffer");
            }
            pos += padding;
        }
        return pos;
    }

    /// Channels per antenna: `OBSNCHAN / NANTS`.
    [[nodiscard]] std::int64_t antnchan() const {
        const auto nants    = get<std::int64_t>("NANTS", std::int64_t{1});
        const auto obsnchan = get<std::int64_t>("OBSNCHAN");
        if (nants == 0 || obsnchan % nants != 0) {
            throw ValidationError("psrio: NANTS must divide OBSNCHAN");
        }
        return obsnchan / nants;
    }

    /// Time samples in one payload block.
    [[nodiscard]] std::int64_t ntime() const {
        const auto npol     = complex_npol();
        const auto nbits    = int_or("NBITS", std::int64_t{8});
        const auto blocsize = get<std::int64_t>("BLOCSIZE");
        const auto obsnchan = get<std::int64_t>("OBSNCHAN");
        if (nbits != 2 && nbits != 4 && nbits != 8 && nbits != 16) {
            throw ValidationError("psrio: GUPPI NBITS must be 2, 4, 8, or 16");
        }
        if (blocsize < 0 || obsnchan <= 0 || npol <= 0) {
            throw ValidationError("psrio: GUPPI block geometry is invalid");
        }
        const auto denom = 2 * obsnchan * npol * nbits;
        if (denom == 0 || (8 * blocsize) % denom != 0) {
            throw ValidationError("psrio: GUPPI BLOCSIZE is not a whole number "
                                  "of complex samples");
        }
        return (8 * blocsize) / denom;
    }

    /// Complex samples in the block (`npol * ntime * OBSNCHAN`).
    [[nodiscard]] std::int64_t blocksize() const {
        const auto nants    = get<std::int64_t>("NANTS", std::int64_t{1});
        const auto obsnchan = get<std::int64_t>("OBSNCHAN");
        const auto samples  = ntime();
        const auto npol     = complex_npol();
        if (nants > 1) {
            if (obsnchan % nants != 0) {
                throw ValidationError("psrio: NANTS must divide OBSNCHAN");
            }
            return npol * samples * (obsnchan / nants) * nants;
        }
        return npol * samples * obsnchan;
    }

    /// Array description such as `Array<Complex<Int8>, 3>`.
    [[nodiscard]] std::string block_type() const {
        const auto ndims =
            get<std::int64_t>("NANTS", std::int64_t{1}) > 1 ? 4 : 3;
        const auto nbits = get<std::int64_t>("NBITS", std::int64_t{8});
        if (nbits != 8 && nbits != 16) {
            throw ValidationError("psrio: GUPPI NBITS must be 8 or 16");
        }
        const char* element = nbits == 8 ? "Complex<Int8>" : "Complex<Int16>";
        return std::format("Array<{}, {}>", element, ndims);
    }

    /// Centre frequency of 1-based channel @p channel, in MHz.
    [[nodiscard]] double channel_frequency(std::int64_t channel) const {
        const auto n       = antnchan();
        if (!has_key("OBSFREQ")) {
            throw FormatError("psrio: GUPPI header missing OBSFREQ");
        }
        const auto obsfreq = real_or("OBSFREQ", 0.0);
        const auto chan_bw = channel_spacing();
        if (n <= 0) {
            throw ValidationError(
                "psrio: GUPPI antenna channel count is empty");
        }
        const auto rel = (((channel - 1) % n) + n) % n;
        return obsfreq - (static_cast<double>(n) * chan_bw / 2.0) +
               (chan_bw * (static_cast<double>(rel) + 0.5));
    }

    /// Frequencies from @p start_channel through @p end_channel, 1-based.
    /// `end_channel == -1` means the last antenna channel.
    [[nodiscard]] std::vector<double>
    channel_frequencies(std::int64_t start_channel,
                        std::int64_t end_channel = -1) const {
        if (end_channel != -1 && start_channel > end_channel) {
            throw ValidationError("psrio: GUPPI channel range is inverted");
        }
        const auto n    = antnchan();
        const auto last = end_channel == -1 ? n : end_channel;
        if (start_channel > last) {
            throw ValidationError("psrio: GUPPI channel range is inverted");
        }
        std::vector<double> freqs;
        freqs.reserve(static_cast<std::size_t>(last - start_channel + 1));
        for (std::int64_t channel = start_channel; channel <= last; ++channel) {
            freqs.push_back(channel_frequency(channel));
        }
        return freqs;
    }

    /// `CHAN_BW` when present, otherwise `OBSBW` divided by the channels on one antenna.
    [[nodiscard]] double channel_spacing() const {
        if (!has_key("CHAN_BW") && !has_key("OBSBW")) {
            throw FormatError("psrio: GUPPI header missing CHAN_BW");
        }
        if (has_key("CHAN_BW")) {
            return real_or("CHAN_BW", 0.0);
        }
        const auto n = antnchan();
        if (n == 0) {
            throw ValidationError("psrio: GUPPI antenna channel count is empty");
        }
        return real_or("OBSBW", 0.0) / static_cast<double>(n);
    }

    /// Fill a stream header from this block. The caller sets the file name,
    /// sample count, and byte cursor.
    [[nodiscard]] BasebandHeader to_baseband_header() const {
        BasebandHeader out;
        out.format   = BasebandFormat::kGuppiRaw;
        out.npol     = static_cast<std::uint64_t>(complex_npol());
        out.nants    = static_cast<std::uint64_t>(int_or("NANTS", std::int64_t{1}));
        out.nchan    = static_cast<std::uint64_t>(antnchan());
        out.nbit     = static_cast<int>(int_or("NBITS", std::int64_t{8}));
        out.ndim     = 2;
        out.samples_signed = true;
        out.msb_first = (out.nbit == 4);
        const auto fmt = string_or("PKTFMT", "");
        out.order = equals_ignore_case(fmt, "SIMPLE") ? BasebandOrder::kTimeMajor
                                                      : BasebandOrder::kChannelMajor;
        out.tsamp = real_or("TBIN", 0.0);
        out.overlap = static_cast<std::uint64_t>(std::max<std::int64_t>(
            int_or("OVERLAP", std::int64_t{0}), std::int64_t{0}));
        out.source    = string_or("SRC_NAME", "Unknown");
        out.telescope = string_or("TELESCOP", "Unknown");
        out.backend   = string_or("BACKEND", "");
        const auto chan_bw = channel_spacing();
        out.foff       = chan_bw;
        out.fch1       = channel_frequency(1);
        out.bandwidth  = std::abs(chan_bw) * static_cast<double>(out.nchan);
        out.center_frequency =
            out.fch1 + (out.foff * (static_cast<double>(out.nchan) - 1.0) / 2.0);
        out.tstart    = start_mjd();
        out.utc_start = out.tstart == 0.0 ? 0 : astro::mjd_to_time(out.tstart);
        out.extra.reserve(m_entries.size());
        for (const auto& [key, value] : m_entries) {
            ExtraKey item;
            item.name = key;
            if (std::holds_alternative<std::string>(value)) {
                item.type  = "string";
                item.value = std::get<std::string>(value);
            } else if (std::holds_alternative<double>(value)) {
                item.type  = "double";
                item.value = std::get<double>(value);
            } else {
                item.type  = "int";
                item.value = std::get<std::int64_t>(value);
            }
            out.extra.push_back(std::move(item));
        }
        return out;
    }

private:
    std::map<std::string, HeaderValue> m_entries;

    static std::string upper(std::string_view key) {
        std::string out(key);
        for (char& character : out) {
            character = static_cast<char>(
                std::toupper(static_cast<unsigned char>(character)));
        }
        return out;
    }

    static std::string trim(std::string_view text) {
        const auto first = text.find_first_not_of(" \t");
        if (first == std::string_view::npos) {
            return {};
        }
        const auto last = text.find_last_not_of(" \t");
        return std::string{text.substr(first, last - first + 1)};
    }

    void add_keyword(const std::string& key, std::string value) {
        if (!value.empty() && (value.front() == '\'' || value.front() == '"')) {
            const char quote = value.front();
            const auto end_quote = value.find(quote, 1);
            if (end_quote != std::string::npos) {
                m_entries[key] = trim(value.substr(1, end_quote - 1));
                return;
            }
        }
        const auto slash = value.find('/');
        if (slash != std::string::npos) {
            value = trim(value.substr(0, slash));
        }
        if (value.empty()) {
            m_entries[key] = std::string{};
            return;
        }
        const bool floating = value.find('.') != std::string::npos ||
                              value.find('e') != std::string::npos ||
                              value.find('E') != std::string::npos;
        try {
            std::size_t processed = 0;
            if (floating) {
                const double dval = std::stod(value, &processed);
                if (processed == value.size()) {
                    m_entries[key] = dval;
                    return;
                }
            } else {
                const auto ival =
                    static_cast<std::int64_t>(std::stoll(value, &processed));
                if (processed == value.size()) {
                    m_entries[key] = ival;
                    return;
                }
            }
        } catch (const std::exception&) {
        }
        m_entries[key] = value;
    }

    [[nodiscard]] double start_mjd() const {
        if (!has_key("STT_IMJD")) {
            return 0.0;
        }
        const auto imjd = static_cast<double>(int_or("STT_IMJD", std::int64_t{0}));
        const auto smjd = static_cast<double>(int_or("STT_SMJD", std::int64_t{0}));
        const auto offs = real_or("STT_OFFS", 0.0);
        double mjd      = imjd + ((smjd + offs) / 86400.0);
        if (has_key("PKTIDX") && has_key("PKTSIZE") && real_or("TBIN", 0.0) > 0.0) {
            const auto pktidx = static_cast<double>(int_or("PKTIDX", std::int64_t{0}));
            const auto pkt    = real_or("PKTSIZE", 0.0);
            const auto bits   = static_cast<double>(antnchan()) *
                                static_cast<double>(int_or("NANTS", std::int64_t{1})) *
                                static_cast<double>(complex_npol()) * 2.0 *
                                static_cast<double>(int_or("NBITS", std::int64_t{8}));
            if (bits > 0.0) {
                const auto seconds = pktidx * pkt * 8.0 * real_or("TBIN", 0.0) / bits;
                mjd += seconds / 86400.0;
            }
        }
        return mjd;
    }

    [[nodiscard]] static bool equals_ignore_case(std::string_view left,
                                                 std::string_view right) {
        if (left.size() != right.size()) {
            return false;
        }
        for (std::size_t index = 0; index < left.size(); ++index) {
            const auto a = static_cast<unsigned char>(left[index]);
            const auto b = static_cast<unsigned char>(right[index]);
            if (std::tolower(a) != std::tolower(b)) {
                return false;
            }
        }
        return true;
    }
};

/// Copy interleaved little-endian complex samples into @p dest.
///
/// @p T is `std::int8_t` or `std::int16_t`. Each complex value is real then
/// imaginary. @p packed must hold `dest.size() * 2 * sizeof(T)` bytes.
template <typename T>
    requires std::same_as<T, std::int8_t> || std::same_as<T, std::int16_t>
inline void unpack_complex(std::span<const std::byte> packed,
                           std::span<std::complex<T>> dest) {
    const auto bytes_per_value = 2U * sizeof(T);
    if (packed.size() != dest.size() * bytes_per_value) {
        throw ValidationError(
            "psrio: GUPPI complex buffer does not match the packed bytes");
    }
    for (std::size_t index = 0; index < dest.size(); ++index) {
        const auto* raw = packed.data() + (index * bytes_per_value);
        const auto real = detail::load_little_endian<T>(raw);
        const auto imag = detail::load_little_endian<T>(raw + sizeof(T));
        dest[index]     = std::complex<T>(real, imag);
    }
}


} // namespace psrio::formats::guppi
