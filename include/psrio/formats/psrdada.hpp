#pragma once

/// PSRDADA ring support. Opt in with this header; psrio.hpp does not include
/// it. Linking `psrio::psrdada` supplies the PSRDADA include path and library.

#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/skip.hpp"
#include "psrio/detail/packed_bits.hpp"

#include <algorithm>
#include <array>
#include <ascii_header.h>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <dada_hdu.h>
#include <format>
#include <ipcbuf.h>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace psrio::formats::psrdada {

namespace detail {

inline bool equal_ignore_case(std::string_view left, std::string_view right) {
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

inline int require_int(const char* header, const char* key) {
    int value = 0;
    if (ascii_header_get(header, key, "%d", &value) < 0) {
        throw ValidationError(
            std::format("psrio: PSRDADA header missing {}", key));
    }
    return value;
}

inline double
require_real(const char* header, const char* key, const char* alternate) {
    double value = 0.0;
    if (ascii_header_get(header, key, "%lf", &value) < 0 &&
        (alternate == nullptr ||
         ascii_header_get(header, alternate, "%lf", &value) < 0)) {
        throw ValidationError(
            std::format("psrio: PSRDADA header missing {}", key));
    }
    return value;
}

inline std::time_t parse_utc_start(const char* text) {
    std::tm broken{};
    if (strptime(text, "%Y-%m-%d-%H:%M:%S", &broken) == nullptr) {
        throw ValidationError(
            "psrio: PSRDADA UTC_START must be YYYY-MM-DD-HH:MM:SS");
    }
    const auto utc = timegm(&broken);
    if (utc == static_cast<std::time_t>(-1)) {
        throw ValidationError("psrio: PSRDADA UTC_START is not a valid time");
    }
    return utc;
}

inline SampleType type_from_nbit(int nbit) {
    switch (nbit) {
    case 1:
        return SampleType::kUInt1;
    case 2:
        return SampleType::kUInt2;
    case 4:
        return SampleType::kUInt4;
    case 8:
        return SampleType::kUInt8;
    case 16:
        return SampleType::kUInt16;
    case 32:
        return SampleType::kUInt32;
    case -32:
        return SampleType::kFloat32;
    default:
        throw ValidationError(
            "psrio: PSRDADA NBIT is not a supported sample width");
    }
}

/// True when the page just opened by ipcbuf_get_next_read is the last one.
///
/// A last page can be completely full, so its byte count matches every other
/// page. The ring records that on the transfer. Checking ipcbuf_eod before
/// the read reports end-of-data while unread pages are still queued.
inline bool page_is_last(const ipcbuf_t* ring) {
    if (ring == nullptr || ring->sync == nullptr || ring->iread < 0 ||
        ring->xfer >= IPCBUF_XFERS) {
        return false;
    }
    const auto reader = static_cast<unsigned>(ring->iread);
    if (reader >= IPCBUF_READERS) {
        return false;
    }
    const auto slot  = static_cast<unsigned>(ring->xfer);
    const auto* sync = ring->sync;
    return sync->eod[slot] != 0 && sync->r_bufs[reader] == sync->e_buf[slot];
}

} // namespace detail

/// One PSRDADA observation, read as packed intensity samples.
///
/// The constructor connects as a reader and blocks until a header is
/// installed. `NPOL` must be 1, `NDIM` must be absent or 1, and `ORDER` must
/// be absent or `TF`. `NBIT` 32 is an unsigned integer. `NBIT` -32 is float.
/// Channels are not reversed.
///
/// A ring is not a file. A released page cannot be read again, so seek and
/// skip backward throw `ValidationError` when the target is before the oldest
/// byte still held (the current page, plus a carry of less than one sample
/// stride from the page before it). Forward skip discards bytes by reading
/// them. A forward skip that runs past end-of-data releases the pages it
/// crossed, then throws.
///
/// `has_nsamples()` stays false until a read sees end-of-data. `nsamples()`
/// stays 0 while the ring is live. `tell()` is the sample cursor.
class RingReader {
public:
    explicit RingReader(std::uint32_t key);

    RingReader(const RingReader&)                = delete;
    RingReader& operator=(const RingReader&)     = delete;
    RingReader(RingReader&&) noexcept            = default;
    RingReader& operator=(RingReader&&) noexcept = default;
    ~RingReader()                                = default;

    [[nodiscard]] std::uint64_t nchans() const noexcept { return m_nchans; }
    [[nodiscard]] std::uint64_t nifs() const noexcept { return 1U; }
    [[nodiscard]] int nbits() const noexcept { return m_nbits; }
    [[nodiscard]] SampleType sample_type() const noexcept { return m_type; }
    [[nodiscard]] std::uint64_t bytes_per_sample() const noexcept {
        return m_stride;
    }

    /// Samples observed once end-of-data has been seen. Zero while live.
    [[nodiscard]] std::uint64_t nsamples() const noexcept { return m_nsamples; }

    /// False until a read reaches end-of-data.
    [[nodiscard]] bool has_nsamples() const noexcept { return m_has_nsamples; }

    /// Sample interval in seconds. The header stores `TSAMP` in microseconds.
    [[nodiscard]] double tsamp() const noexcept { return m_tsamp; }

    /// Zero. The ring header carries `UTC_START` rather than an MJD.
    [[nodiscard]] double tstart() const noexcept { return 0.0; }

    [[nodiscard]] std::string_view source_name() const noexcept {
        return m_source_name;
    }

    [[nodiscard]] std::string_view telescope() const noexcept {
        return m_telescope;
    }

    [[nodiscard]] double raj() const noexcept { return 0.0; }
    [[nodiscard]] double dej() const noexcept { return 0.0; }

    /// Centre frequency of channel 0, in MHz.
    [[nodiscard]] double fch1() const noexcept { return m_fch1; }

    /// Signed channel spacing, in MHz.
    [[nodiscard]] double foff() const noexcept { return m_foff; }

    [[nodiscard]] double bandwidth() const noexcept {
        return std::abs(m_foff) * static_cast<double>(m_nchans);
    }

    [[nodiscard]] double center_frequency() const noexcept {
        return m_fch1 + (m_foff * (static_cast<double>(m_nchans) - 1.0) / 2.0);
    }

    /// Beam index from the header, or 0 when `BEAM` is absent.
    [[nodiscard]] int beam() const noexcept { return m_beam; }

    /// Spectra per second.
    [[nodiscard]] double spectra_rate() const noexcept {
        return m_spectra_rate;
    }

    /// `UTC_START` as POSIX seconds.
    [[nodiscard]] std::time_t utc_start() const noexcept { return m_utc; }

    [[nodiscard]] std::uint64_t tell() const noexcept {
        return m_stride == 0U ? 0U : m_byte / m_stride;
    }

    void set_fswap(bool enable) noexcept { m_apply_fswap = enable; }
    [[nodiscard]] bool fswap_enabled() const noexcept { return m_apply_fswap; }

    /// Move to @p sample. @p sample may equal nsamples() after end-of-data.
    /// @throws ValidationError if @p sample is outside the bytes still held,
    ///         or past the observation once its length is known.
    void seek(std::uint64_t sample);

    void rewind() { seek(0U); }

    /// Move the cursor by @p delta time samples. Negative moves backward.
    /// @throws ValidationError if the landing index is before the held bytes
    ///         or past the readable samples.
    void skip(std::int64_t delta);

    /// Copy up to @p count packed time samples into @p dest.
    ///
    /// @p dest must hold at least `count * bytes_per_sample()` bytes. The
    /// copy is raw ring order. Returns the number of time samples copied,
    /// which is short at end-of-data and zero when the cursor is already there.
    /// @throws ValidationError if @p dest is smaller than the requested block.
    std::uint64_t read_block(std::uint64_t count, std::span<std::byte> dest);

    /// Unpack the next @p count time samples into @p dest as float.
    ///
    /// @p dest must hold `count * nchans` values. Sub-byte samples are
    /// least-significant field first. 32-bit integer samples are not floats.
    /// A short read leaves the unused tail of @p dest untouched.
    /// @throws ValidationError if @p dest has the wrong size.
    std::uint64_t read_samples(std::uint64_t count, std::span<float> dest);

    /// Unpack next @p count time samples into uint8_t @p dest.
    std::uint64_t read_samples(std::uint64_t count, std::span<std::uint8_t> dest);

    /// Unpack next @p count time samples into uint16_t @p dest.
    std::uint64_t read_samples(std::uint64_t count, std::span<std::uint16_t> dest);

    /// Read raw packed bytes directly into caller storage.
    std::uint64_t read_bytes(std::uint64_t nbytes, std::span<std::byte> dest);

    /// Convenience allocating read: copy next @p count packed time samples.
    [[nodiscard]] std::vector<std::byte> read_block(std::uint64_t count);

    /// Convenience allocating read: unpack next @p count time samples into a new vector.
    template <typename T = float>
    [[nodiscard]] std::vector<T> read_samples(std::uint64_t count);

    /// Convenience allocating read: copy next @p nbytes payload bytes.
    [[nodiscard]] std::vector<std::byte> read_bytes(std::uint64_t nbytes);

private:
    class HduGuard {
    public:
        HduGuard()                           = default;
        HduGuard(const HduGuard&)            = delete;
        HduGuard& operator=(const HduGuard&) = delete;

        HduGuard(HduGuard&& other) noexcept
            : m_hdu(std::exchange(other.m_hdu, nullptr)),
              m_header(std::exchange(other.m_header, nullptr)),
              m_data(std::exchange(other.m_data, nullptr)),
              m_connected(std::exchange(other.m_connected, false)),
              m_locked(std::exchange(other.m_locked, false)),
              m_header_open(std::exchange(other.m_header_open, false)),
              m_holding_page(std::exchange(other.m_holding_page, false)) {}

        HduGuard& operator=(HduGuard&& other) noexcept {
            if (this != &other) {
                release();
                m_hdu          = std::exchange(other.m_hdu, nullptr);
                m_header       = std::exchange(other.m_header, nullptr);
                m_data         = std::exchange(other.m_data, nullptr);
                m_connected    = std::exchange(other.m_connected, false);
                m_locked       = std::exchange(other.m_locked, false);
                m_header_open  = std::exchange(other.m_header_open, false);
                m_holding_page = std::exchange(other.m_holding_page, false);
            }
            return *this;
        }

        ~HduGuard() { release(); }

        void open(std::uint32_t key) {
            m_hdu = dada_hdu_create(nullptr);
            if (m_hdu == nullptr) {
                throw IoError("psrio: failed to create a PSRDADA HDU");
            }
            dada_hdu_set_key(m_hdu, static_cast<key_t>(key));
            if (dada_hdu_connect(m_hdu) < 0) {
                throw IoError(std::format(
                    "psrio: failed to connect to PSRDADA ring {:#x}", key));
            }
            m_connected = true;
            if (dada_hdu_lock_read(m_hdu) < 0) {
                throw IoError(std::format(
                    "psrio: failed to lock PSRDADA ring {:#x} for reading",
                    key));
            }
            m_locked = true;
            m_header = m_hdu->header_block;
            m_data   = &m_hdu->data_block->buf;
        }

        void set_header_open(bool open) noexcept { m_header_open = open; }
        void set_holding(bool holding) noexcept { m_holding_page = holding; }

        [[nodiscard]] ipcbuf_t* header() const noexcept { return m_header; }
        [[nodiscard]] ipcbuf_t* data() const noexcept { return m_data; }

    private:
        void release() noexcept {
            if (m_header_open && m_header != nullptr) {
                ipcbuf_mark_cleared(m_header);
                m_header_open = false;
            }
            if (m_holding_page && m_data != nullptr) {
                ipcbuf_mark_cleared(m_data);
                m_holding_page = false;
            }
            if (m_hdu == nullptr) {
                return;
            }
            if (m_locked) {
                dada_hdu_unlock_read(m_hdu);
                m_locked = false;
            }
            if (m_connected) {
                dada_hdu_disconnect(m_hdu);
                m_connected = false;
            }
            dada_hdu_destroy(m_hdu);
            m_hdu    = nullptr;
            m_header = nullptr;
            m_data   = nullptr;
        }

        dada_hdu_t* m_hdu   = nullptr;
        ipcbuf_t* m_header  = nullptr;
        ipcbuf_t* m_data    = nullptr;
        bool m_connected    = false;
        bool m_locked       = false;
        bool m_header_open  = false;
        bool m_holding_page = false;
    };

    void parse_header(const char* header);
    void seek_forward(std::uint64_t target);
    void fetch_page();
    void release_page();
    void retain_open_sample();
    void note_eod() noexcept;
    std::uint64_t pull(std::span<std::byte> dest);
    const std::byte* bytes_at_cursor(std::uint64_t& available) const;
    bool holds(std::uint64_t absolute) const noexcept;
    template <typename T>
    void unpack_to(std::span<const std::byte> packed,
                   std::span<T> dest) const;
    template <typename T>
    std::uint64_t read_samples_impl(std::uint64_t count,
                                    std::span<T> dest);
    void copy_held(std::uint64_t begin,
                   std::uint64_t end,
                   std::vector<std::byte>& out) const;

    [[nodiscard]] std::uint64_t page_end() const noexcept {
        return m_page_origin + m_page_size;
    }

    static constexpr std::size_t kDiscardBytes = 4096U;

    HduGuard m_guard;
    ipcbuf_t* m_data            = nullptr;
    char* m_page                = nullptr;
    std::uint64_t m_page_origin = 0;
    std::uint64_t m_page_size   = 0;
    bool m_page_final           = false;
    std::uint64_t m_next_origin = 0;
    std::vector<std::byte> m_carry;
    std::uint64_t m_carry_origin = 0;
    std::uint64_t m_byte         = 0;

    std::string m_source_name{"Unknown"};
    std::string m_telescope{"Unknown"};
    bool m_apply_fswap{false};

    std::uint64_t m_nchans   = 0;
    int m_nbits              = 0;
    SampleType m_type        = SampleType::kUInt8;
    std::uint64_t m_stride   = 0;
    int m_beam               = 0;
    double m_tsamp           = 0.0;
    double m_fch1            = 0.0;
    double m_foff            = 0.0;
    double m_spectra_rate    = 0.0;
    std::time_t m_utc        = 0;
    bool m_has_nsamples      = false;
    std::uint64_t m_nsamples = 0;
    bool m_eod               = false;
};

inline RingReader::RingReader(std::uint32_t key) {
    m_guard.open(key);
    m_data                     = m_guard.data();
    std::uint64_t header_bytes = 0;
    char* header = ipcbuf_get_next_read(m_guard.header(), &header_bytes);
    if (header == nullptr) {
        throw IoError(
            "psrio: PSRDADA ring closed before a header was available");
    }
    m_guard.set_header_open(true);
    try {
        parse_header(header);
    } catch (...) {
        if (ipcbuf_mark_cleared(m_guard.header()) == 0) {
            m_guard.set_header_open(false);
        }
        throw;
    }
    if (ipcbuf_mark_cleared(m_guard.header()) < 0) {
        throw IoError("psrio: PSRDADA could not release the header page");
    }
    m_guard.set_header_open(false);
}

inline void RingReader::parse_header(const char* header) {
    const int nchan = detail::require_int(header, "NCHAN");
    const int nbit  = detail::require_int(header, "NBIT");
    const int npol  = detail::require_int(header, "NPOL");
    if (nchan <= 0) {
        throw ValidationError("psrio: PSRDADA NCHAN must be positive");
    }
    if (npol != 1) {
        throw ValidationError(
            "psrio: PSRDADA NPOL must be 1 (detected intensity)");
    }
    int ndim = 1;
    if (ascii_header_get(header, "NDIM", "%d", &ndim) >= 0 && ndim != 1) {
        throw ValidationError(
            "psrio: PSRDADA NDIM must be 1 (real intensity samples)");
    }
    std::array<char, 16> order{};
    if (ascii_header_get(header, "ORDER", "%15s", order.data()) >= 0 &&
        !detail::equal_ignore_case(order.data(), "TF")) {
        throw ValidationError(
            "psrio: PSRDADA ORDER must be TF (time then frequency)");
    }

    std::array<char, 64> src{};
    if (ascii_header_get(header, "SRC_NAME", "%63s", src.data()) >= 0 ||
        ascii_header_get(header, "SOURCE", "%63s", src.data()) >= 0) {
        m_source_name = src.data();
    } else {
        m_source_name = "Unknown";
    }
    std::array<char, 64> tel{};
    if (ascii_header_get(header, "TELESCOPE", "%63s", tel.data()) >= 0) {
        m_telescope = tel.data();
    } else {
        m_telescope = "Unknown";
    }

    const double bandwidth = detail::require_real(header, "BANDWIDTH", "BW");
    const double cfreq     = detail::require_real(header, "CFREQ", "FREQ");
    if (!std::isfinite(bandwidth) || !std::isfinite(cfreq)) {
        throw ValidationError(
            "psrio: PSRDADA frequency keywords must be finite");
    }
    double tsamp_us = 0.0;
    if (ascii_header_get(header, "TSAMP", "%lf", &tsamp_us) < 0 ||
        !std::isfinite(tsamp_us) || !(tsamp_us > 0.0)) {
        throw ValidationError(
            "psrio: PSRDADA TSAMP must be positive microseconds");
    }
    std::array<char, 64> utc{};
    if (ascii_header_get(header, "UTC_START", "%63s", utc.data()) < 0) {
        throw ValidationError("psrio: PSRDADA header missing UTC_START");
    }
    int beam = 0;
    if (ascii_header_get(header, "BEAM", "%d", &beam) < 0) {
        beam = 0;
    }
    if (beam < 0) {
        throw ValidationError("psrio: PSRDADA BEAM must be >= 0");
    }

    m_type           = detail::type_from_nbit(nbit);
    m_nbits          = nbit;
    m_nchans         = static_cast<std::uint64_t>(nchan);
    const auto width = static_cast<std::uint64_t>(bits_of(m_type));
    if (width == 0U ||
        m_nchans > std::numeric_limits<std::uint64_t>::max() / width) {
        throw ValidationError("psrio: PSRDADA channel geometry is too large");
    }
    const auto bits = m_nchans * width;
    if (bits % 8U != 0U) {
        throw ValidationError(
            "psrio: PSRDADA NCHAN * NBIT must be a positive multiple of 8");
    }
    m_stride = bits / 8U;
    if (m_stride == 0U) {
        throw ValidationError(
            "psrio: PSRDADA NCHAN * NBIT must be a positive multiple of 8");
    }

    const auto channels  = static_cast<double>(nchan);
    const double spacing = bandwidth / channels;
    const double start   = cfreq - (bandwidth / 2.0);
    m_fch1               = start + (spacing / 2.0);
    m_foff               = spacing;
    m_tsamp              = tsamp_us / 1.0e6;
    m_spectra_rate       = 1.0e6 / tsamp_us;
    m_utc                = detail::parse_utc_start(utc.data());
    m_beam               = beam;
}

inline void RingReader::note_eod() noexcept {
    m_eod = true;
    if (!m_has_nsamples) {
        m_has_nsamples = true;
        m_nsamples     = m_stride == 0U ? 0U : m_byte / m_stride;
    }
}

inline bool RingReader::holds(std::uint64_t absolute) const noexcept {
    if (!m_carry.empty() && absolute >= m_carry_origin &&
        absolute <= m_carry_origin + m_carry.size()) {
        return true;
    }
    if (m_page != nullptr && absolute >= m_page_origin &&
        absolute <= page_end()) {
        return true;
    }
    return absolute == m_byte && m_page == nullptr && m_carry.empty();
}

inline const std::byte*
RingReader::bytes_at_cursor(std::uint64_t& available) const {
    if (!m_carry.empty() && m_byte >= m_carry_origin &&
        m_byte < m_carry_origin + m_carry.size()) {
        const auto offset = static_cast<std::size_t>(m_byte - m_carry_origin);
        available         = m_carry.size() - offset;
        return m_carry.data() + offset;
    }
    if (m_page != nullptr && m_byte >= m_page_origin && m_byte < page_end()) {
        const auto offset = static_cast<std::size_t>(m_byte - m_page_origin);
        available         = m_page_size - offset;
        return reinterpret_cast<const std::byte*>(m_page) + offset;
    }
    available = 0;
    return nullptr;
}

inline void RingReader::copy_held(std::uint64_t begin,
                                  std::uint64_t end,
                                  std::vector<std::byte>& out) const {
    const auto append = [&](const std::byte* data, std::uint64_t origin,
                            std::uint64_t size) {
        if (data == nullptr || size == 0U || begin >= end) {
            return;
        }
        const auto region_end = origin + size;
        const auto lower      = std::max(begin, origin);
        const auto upper      = std::min(end, region_end);
        if (lower >= upper) {
            return;
        }
        const auto* source = data + static_cast<std::size_t>(lower - origin);
        const auto count   = static_cast<std::size_t>(upper - lower);
        out.insert(out.end(), source, source + count);
    };
    append(m_carry.empty() ? nullptr : m_carry.data(), m_carry_origin,
           m_carry.size());
    append(m_page == nullptr ? nullptr
                             : reinterpret_cast<const std::byte*>(m_page),
           m_page_origin, m_page_size);
}

inline void RingReader::retain_open_sample() {
    const auto end = page_end();
    if (m_stride == 0U || end % m_stride == 0U) {
        m_carry.clear();
        return;
    }
    const auto sample_start = end - (end % m_stride);
    const auto nbytes       = end - sample_start;
    std::vector<std::byte> kept;
    kept.reserve(static_cast<std::size_t>(nbytes));
    copy_held(sample_start, end, kept);
    if (kept.size() != nbytes) {
        throw FormatError("psrio: PSRDADA lost bytes at a page boundary");
    }
    m_carry        = std::move(kept);
    m_carry_origin = sample_start;
}

inline void RingReader::release_page() {
    if (m_page == nullptr) {
        return;
    }
    retain_open_sample();
    if (ipcbuf_mark_cleared(m_data) < 0) {
        throw IoError("psrio: PSRDADA could not release a ring page");
    }
    m_page       = nullptr;
    m_page_size  = 0;
    m_page_final = false;
    m_guard.set_holding(false);
}

inline void RingReader::fetch_page() {
    if (m_data == nullptr) {
        throw ValidationError("psrio: PSRDADA reader is empty");
    }
    if (m_eod || m_page != nullptr) {
        return;
    }
    std::uint64_t nbytes = 0;
    char* page           = ipcbuf_get_next_read(m_data, &nbytes);
    if (page == nullptr) {
        if (ipcbuf_eod(m_data) > 0) {
            note_eod();
            return;
        }
        throw IoError("psrio: PSRDADA read failed");
    }
    m_page        = page;
    m_page_origin = m_next_origin;
    m_page_size   = nbytes;
    m_next_origin += nbytes;
    m_page_final = detail::page_is_last(m_data);
    m_guard.set_holding(true);
    if (nbytes != 0U) {
        return;
    }
    if (ipcbuf_mark_cleared(m_data) < 0) {
        throw IoError("psrio: PSRDADA could not release a ring page");
    }
    m_page       = nullptr;
    m_page_size  = 0;
    m_page_final = false;
    m_guard.set_holding(false);
    note_eod();
}

inline std::uint64_t RingReader::pull(std::span<std::byte> dest) {
    std::uint64_t filled = 0;
    while (filled < dest.size()) {
        std::uint64_t have = 0;
        const auto* source = bytes_at_cursor(have);
        if (m_page_final && have > 0U && m_byte < page_end()) {
            const auto remaining = page_end() - m_byte;
            if (remaining < m_stride) {
                note_eod();
                break;
            }
            const auto whole = remaining - (remaining % m_stride);
            have             = std::min(have, whole);
        }
        if (have == 0U) {
            if (m_page != nullptr && m_byte >= page_end()) {
                if (m_page_final) {
                    note_eod();
                    break;
                }
                release_page();
            } else if (m_page != nullptr) {
                throw FormatError(
                    "psrio: PSRDADA cursor is outside the held page");
            }
            if (m_eod) {
                break;
            }
            fetch_page();
            if (m_eod || m_page == nullptr) {
                break;
            }
            continue;
        }
        const auto need  = static_cast<std::uint64_t>(dest.size()) - filled;
        const auto count = std::min(have, need);
        std::memcpy(dest.data() + filled, source,
                    static_cast<std::size_t>(count));
        m_byte += count;
        filled += count;
        if (m_page_final && m_page != nullptr && m_byte <= page_end() &&
            page_end() - m_byte < m_stride) {
            note_eod();
        }
    }
    return filled;
}

inline void RingReader::seek_forward(std::uint64_t target) {
    std::array<std::byte, kDiscardBytes> sink{};
    while (m_byte < target) {
        const auto chunk =
            std::min(target - m_byte, static_cast<std::uint64_t>(sink.size()));
        const auto got = pull(
            std::span<std::byte>(sink.data(), static_cast<std::size_t>(chunk)));
        if (got < chunk) {
            throw ValidationError("psrio: seek is past the readable samples");
        }
    }
}

inline void RingReader::seek(std::uint64_t sample) {
    if (m_stride == 0U) {
        throw ValidationError("psrio: PSRDADA reader is empty");
    }
    if (sample > std::numeric_limits<std::uint64_t>::max() / m_stride) {
        throw ValidationError("psrio: seek is past the readable samples");
    }
    if (m_has_nsamples && sample > m_nsamples) {
        throw ValidationError("psrio: seek is past the readable samples");
    }
    const auto target = sample * m_stride;
    if (target == m_byte) {
        return;
    }
    if (target < m_byte || holds(target)) {
        if (!holds(target)) {
            throw ValidationError(
                "psrio: seek lands before the oldest byte still held in the "
                "current page");
        }
        m_byte = target;
        return;
    }
    seek_forward(target);
}

inline void RingReader::skip(std::int64_t delta) {
    if (delta == 0) {
        return;
    }
    if (m_has_nsamples) {
        seek(::psrio::detail::apply_skip(tell(), delta, m_nsamples));
        return;
    }
    if (delta > 0) {
        const auto step = static_cast<std::uint64_t>(delta);
        if (step > std::numeric_limits<std::uint64_t>::max() - tell()) {
            throw ValidationError(
                "psrio: skip lands past the readable samples");
        }
        seek(tell() + step);
        return;
    }
    if (delta == std::numeric_limits<std::int64_t>::min() ||
        static_cast<std::uint64_t>(-delta) > tell()) {
        throw ValidationError("psrio: skip lands before the first sample");
    }
    seek(tell() - static_cast<std::uint64_t>(-delta));
}

inline std::uint64_t RingReader::read_block(std::uint64_t count,
                                            std::span<std::byte> dest) {
    if (m_stride == 0U) {
        throw ValidationError("psrio: PSRDADA reader is empty");
    }
    if (count > 0U &&
        m_stride > (std::numeric_limits<std::uint64_t>::max() / count)) {
        throw ValidationError("psrio: requested sample block is too large");
    }
    const auto bytes_needed = count * m_stride;
    if (dest.size() < bytes_needed) {
        throw ValidationError(std::format(
            "psrio: destination has {} bytes but the block needs {}",
            dest.size(), bytes_needed));
    }
    if (count == 0U) {
        return 0U;
    }
    const auto got = pull(dest.first(static_cast<std::size_t>(bytes_needed)));
    if (got % m_stride != 0U) {
        throw FormatError("psrio: PSRDADA observation ended inside a sample");
    }
    const auto samples = got / m_stride;
    if (m_apply_fswap && samples > 0U) {
        reverse_channels(dest.first(static_cast<std::size_t>(got)),
                         samples, m_nchans, m_nbits);
    }
    return samples;
}

template <typename T>
inline void RingReader::unpack_to(std::span<const std::byte> packed,
                                  std::span<T> dest) const {
    if constexpr (std::same_as<T, std::uint8_t>) {
        if (m_nbits > 8) {
            throw ValidationError(
                "psrio: uint8_t output requires nbits of 8 or less");
        }
    } else if constexpr (std::same_as<T, std::uint16_t>) {
        if (m_nbits != 16) {
            throw ValidationError(
                "psrio: uint16_t output requires 16-bit samples");
        }
    }

    const auto samples = packed.size() / static_cast<std::size_t>(m_stride);
    const auto values  = samples * static_cast<std::size_t>(m_nchans);
    if (dest.size() < values) {
        throw ValidationError(
            "psrio: destination is smaller than the unpacked samples");
    }
    const auto output = dest.first(values);
    if (m_nbits <= 4) {
        if constexpr (std::same_as<T, std::uint16_t>) {
            throw ValidationError("psrio: uint16_t output requires 16-bit samples");
        } else {
            ::psrio::detail::unpack_sub_byte(packed, output, m_nbits,
                                             BitOrder::kLsbFirst);
        }
    } else if (m_nbits == 8) {
        if constexpr (std::same_as<T, std::uint16_t>) {
            throw ValidationError("psrio: uint16_t output requires 16-bit samples");
        } else {
            ::psrio::detail::unpack_8bit(packed, output, false);
        }
    } else if (m_nbits == 16) {
        if constexpr (std::same_as<T, std::uint8_t>) {
            throw ValidationError(
                "psrio: uint8_t output requires nbits of 8 or less");
        } else {
            ::psrio::detail::unpack_16le(packed, output);
        }
    } else if (m_type == SampleType::kUInt32) {
        if constexpr (std::same_as<T, float>) {
            for (std::size_t index = 0; index < values; ++index) {
                output[index] = static_cast<float>(
                    ::psrio::detail::load_little_endian<std::uint32_t>(
                        packed.data() + (index * sizeof(std::uint32_t))));
            }
        } else {
            throw ValidationError("psrio: 32-bit samples unpack to float");
        }
    } else if constexpr (std::same_as<T, float>) {
        ::psrio::detail::unpack_32le(packed, output);
    } else {
        throw FormatError("psrio: PSRDADA sample type cannot be unpacked");
    }
}

template <typename T>
inline std::uint64_t RingReader::read_samples_impl(std::uint64_t count,
                                                   std::span<T> dest) {
    if (m_stride == 0U || m_nchans == 0U) {
        throw ValidationError("psrio: PSRDADA reader is empty");
    }
    if (count > 0U &&
        m_nchans > (std::numeric_limits<std::uint64_t>::max() / count)) {
        throw ValidationError("psrio: requested sample block is too large");
    }
    const auto expected = count * m_nchans;
    if (expected != dest.size()) {
        throw ValidationError(std::format(
            "psrio: destination has {} values but the request needs {}",
            dest.size(), expected));
    }
    if (count == 0U) {
        return 0U;
    }
    if (m_stride > std::numeric_limits<std::uint64_t>::max() / count) {
        throw ValidationError("psrio: requested sample block is too large");
    }
    const auto bytes_needed = count * m_stride;
    std::uint64_t have      = 0;
    const auto* source      = bytes_at_cursor(have);
    if (have == 0U && m_page == nullptr && !m_eod) {
        fetch_page();
        source = bytes_at_cursor(have);
    }
    auto whole = have;
    if (m_page_final && m_page != nullptr && m_byte < page_end()) {
        const auto remaining = page_end() - m_byte;
        whole = remaining < m_stride ? 0U : remaining - (remaining % m_stride);
        whole = std::min(whole, have);
    }
    if (source != nullptr && whole >= bytes_needed) {
        unpack_to(std::span<const std::byte>(
                      source, static_cast<std::size_t>(bytes_needed)),
                  dest);
        m_byte += bytes_needed;
        if (m_page_final && m_page != nullptr && m_byte <= page_end() &&
            page_end() - m_byte < m_stride) {
            note_eod();
        }
        return count;
    }

    std::vector<std::byte> packed(static_cast<std::size_t>(bytes_needed));
    const auto got = pull(packed);
    if (got % m_stride != 0U) {
        throw FormatError("psrio: PSRDADA observation ended inside a sample");
    }
    const auto nread = got / m_stride;
    if (nread == 0U) {
        return 0U;
    }
    unpack_to(std::span<const std::byte>(
                  packed.data(), static_cast<std::size_t>(nread * m_stride)),
              dest.first(static_cast<std::size_t>(nread * m_nchans)));
    return nread;
}

inline std::uint64_t RingReader::read_samples(std::uint64_t count,
                                              std::span<float> dest) {
    return read_samples_impl(count, dest);
}

inline std::uint64_t RingReader::read_samples(std::uint64_t count,
                                              std::span<std::uint8_t> dest) {
    return read_samples_impl(count, dest);
}

inline std::uint64_t RingReader::read_samples(std::uint64_t count,
                                              std::span<std::uint16_t> dest) {
    return read_samples_impl(count, dest);
}

inline std::uint64_t RingReader::read_bytes(std::uint64_t nbytes,
                                            std::span<std::byte> dest) {
    if (dest.size() != nbytes) {
        throw ValidationError(std::format(
            "psrio: destination has {} bytes but the request needs {}",
            dest.size(), nbytes));
    }
    if (nbytes == 0U) {
        return 0U;
    }
    if (m_stride == 0U || nbytes % m_stride != 0U) {
        throw ValidationError(
            "psrio: byte request must be a positive multiple of the sample stride");
    }
    const auto samples = nbytes / m_stride;
    const auto read_count = read_block(samples, dest);
    return read_count * m_stride;
}

inline std::vector<std::byte> RingReader::read_block(std::uint64_t count) {
    const auto bytes_needed = count * m_stride;
    std::vector<std::byte> out(static_cast<std::size_t>(bytes_needed));
    const auto actual = read_block(count, std::span<std::byte>(out));
    out.resize(static_cast<std::size_t>(actual * m_stride));
    return out;
}

inline std::vector<std::byte> RingReader::read_bytes(std::uint64_t nbytes) {
    std::vector<std::byte> out(static_cast<std::size_t>(nbytes));
    const auto actual = read_bytes(nbytes, std::span<std::byte>(out));
    out.resize(static_cast<std::size_t>(actual));
    return out;
}

template <typename T>
inline std::vector<T> RingReader::read_samples(std::uint64_t count) {
    const auto values_needed = count * m_nchans;
    std::vector<T> out(static_cast<std::size_t>(values_needed));
    const auto actual = read_samples(count, std::span<T>(out));
    out.resize(static_cast<std::size_t>(actual * m_nchans));
    return out;
}

} // namespace psrio::formats::psrdada
