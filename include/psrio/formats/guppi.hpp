#pragma once

#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/mmap.hpp"

#include <algorithm>
#include <cctype>
#include <complex>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <format>
#include <map>
#include <span>
#include <stdexcept>
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

    template <typename T> T get(std::string_view key) const {
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

    template <typename T> T get(std::string_view key, const T& fallback) const {
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
        const auto npol     = get<std::int64_t>("NPOL", std::int64_t{1}) < 2
                                  ? std::int64_t{1}
                                  : std::int64_t{2};
        const auto nbits    = get<std::int64_t>("NBITS", std::int64_t{8});
        const auto blocsize = get<std::int64_t>("BLOCSIZE");
        const auto obsnchan = get<std::int64_t>("OBSNCHAN");
        if (nbits != 8 && nbits != 16) {
            throw ValidationError("psrio: GUPPI NBITS must be 8 or 16");
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
        const auto npol     = get<std::int64_t>("NPOL", std::int64_t{1}) < 2
                                  ? std::int64_t{1}
                                  : std::int64_t{2};
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
        const auto obsfreq = get<double>("OBSFREQ");
        const auto chan_bw = get<double>("CHAN_BW");
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

    void add_keyword(const std::string& key, const std::string& value) {
        if (value.size() >= 2 &&
            (value.front() == '\'' || value.front() == '"')) {
            m_entries[key] = value.substr(1, value.size() - 2);
            return;
        }
        const bool floating = value.find('.') != std::string::npos ||
                              value.find('e') != std::string::npos ||
                              value.find('E') != std::string::npos;
        try {
            if (floating) {
                m_entries[key] = std::stod(value);
            } else {
                m_entries[key] = static_cast<std::int64_t>(std::stoll(value));
            }
            return;
        } catch (const std::invalid_argument&) {
            m_entries[key] = value;
        } catch (const std::out_of_range&) {
            m_entries[key] = value;
        }
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

/// Sequential reader for a GUPPI RAW file.
///
/// The file alternates an ASCII header, optional DIRECTIO padding, and one
/// `BLOCSIZE` payload. `tell` and `seek` use file byte offsets. A seek that
/// lands inside the current payload resumes the byte read. A seek outside
/// that window waits for `read_header`. This reader does not model
/// `concepts::BlockReader`.
class RawReader {
public:
    explicit RawReader(const std::filesystem::path& path) {
        try {
            m_file = detail::MappedFile(path);
            m_file.advise_sequential();
        } catch (const std::exception& ex) {
            throw IoError(ex.what());
        }
        m_bytes = m_file.bytes();
    }

    RawReader(const RawReader&)                = delete;
    RawReader& operator=(const RawReader&)     = delete;
    RawReader(RawReader&&) noexcept            = default;
    RawReader& operator=(RawReader&&) noexcept = default;
    ~RawReader()                               = default;

    /// Byte offset of the next read.
    [[nodiscard]] std::uint64_t tell() const noexcept { return m_pos; }

    /// Move to a file byte offset. The current payload window stays armed
    /// when @p byte_offset still lies inside it.
    void seek(std::uint64_t byte_offset) {
        if (byte_offset > m_bytes.size()) {
            throw ValidationError(
                "psrio: GUPPI seek is past the end of the file");
        }
        m_pos        = byte_offset;
        m_in_payload = m_has_payload && byte_offset >= m_payload_begin &&
                       byte_offset < m_payload_end;
    }

    void rewind() { seek(0); }

    /// Parse the header at the current offset and leave the cursor on the
    /// payload that follows it. The cursor does not move when parsing fails.
    GuppiHeader read_header() {
        if (m_in_payload) {
            throw ValidationError(
                "psrio: GUPPI payload remains; finish it or seek before the "
                "next header");
        }
        if (m_pos > m_bytes.size()) {
            throw FormatError(
                "psrio: GUPPI cursor is past the end of the file");
        }
        GuppiHeader header;
        const auto consumed =
            header.parse(m_bytes.subspan(static_cast<std::size_t>(m_pos)));
        const auto header_end = m_pos + consumed;
        if (!header.has_key("BLOCSIZE")) {
            m_pos           = header_end;
            m_has_payload   = false;
            m_in_payload    = false;
            m_payload_begin = 0;
            m_payload_end   = 0;
            return header;
        }
        const auto bloc = header.get<std::int64_t>("BLOCSIZE");
        if (bloc < 0) {
            throw ValidationError("psrio: GUPPI BLOCSIZE must be >= 0");
        }
        const auto payload = static_cast<std::uint64_t>(bloc);
        if (payload > m_bytes.size() - header_end) {
            throw FormatError(
                "psrio: GUPPI payload extends past the end of the file");
        }
        m_pos           = header_end;
        m_payload_begin = header_end;
        m_payload_end   = header_end + payload;
        m_has_payload   = true;
        m_in_payload    = payload > 0U;
        return header;
    }

    /// Copy the next payload bytes. The return count is short at the end of
    /// the current payload. @p dest must have size @p nbytes.
    std::uint64_t read_bytes(std::uint64_t nbytes, std::span<std::byte> dest) {
        if (dest.size() != nbytes) {
            throw ValidationError(
                "psrio: destination size does not match requested bytes");
        }
        if (!m_in_payload) {
            throw ValidationError(
                "psrio: GUPPI read_bytes needs a header before the payload");
        }
        const auto available = m_payload_end - m_pos;
        const auto to_copy   = std::min(nbytes, available);
        if (to_copy > 0U) {
            std::memcpy(dest.data(),
                        m_bytes.data() + static_cast<std::size_t>(m_pos),
                        static_cast<std::size_t>(to_copy));
        }
        m_pos += to_copy;
        if (m_pos == m_payload_end) {
            m_in_payload = false;
        }
        return to_copy;
    }

    [[nodiscard]] std::vector<std::byte> read_bytes(std::uint64_t nbytes) {
        std::vector<std::byte> out(static_cast<std::size_t>(nbytes));
        const auto actual = read_bytes(nbytes, std::span<std::byte>(out));
        out.resize(static_cast<std::size_t>(actual));
        return out;
    }

    /// Zero-copy view of the next @p nbytes payload bytes.
    [[nodiscard]] std::span<const std::byte> view_bytes(std::uint64_t nbytes) {
        if (!m_in_payload) {
            throw ValidationError(
                "psrio: GUPPI view_bytes needs a header before the payload");
        }
        const auto available = m_payload_end - m_pos;
        if (nbytes > available) {
            throw ValidationError(
                "psrio: GUPPI view extends past the current payload");
        }
        const auto view = m_bytes.subspan(static_cast<std::size_t>(m_pos),
                                          static_cast<std::size_t>(nbytes));
        m_pos += nbytes;
        if (m_pos == m_payload_end) {
            m_in_payload = false;
        }
        return view;
    }

private:
    detail::MappedFile m_file;
    std::span<const std::byte> m_bytes;
    std::uint64_t m_pos{0};
    std::uint64_t m_payload_begin{0};
    std::uint64_t m_payload_end{0};
    bool m_has_payload{false};
    bool m_in_payload{false};
};

} // namespace psrio::formats::guppi
