#pragma once

#include "psrio/formats/guppi/header.hpp"

#include "psrio/detail/mmap.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <exception>
#include <filesystem>
#include <format>
#include <limits>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

namespace psrio::formats::guppi {

class GuppiSet;

/// Sequential byte reader for a GUPPI RAW file.
///
/// The file alternates an ASCII header, optional DIRECTIO padding, and one
/// `BLOCSIZE` payload. `tell` and `seek` use file byte offsets. A seek that
/// lands inside the current payload resumes the byte read. A seek outside
/// that window waits for `read_header`. This reader does not model
/// `concepts::BlockReader` or `concepts::BasebandReader`.
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

/// One GUPPI band, possibly split across time-ordered files.
///
/// `read_block` and `read_samples` return canonical time-major voltages.
/// Overlap samples are dropped. A `PKTIDX` jump raises `FormatError`.
class GuppiReader {
public:
    explicit GuppiReader(const std::filesystem::path& path) {
        open_one(path);
        if (m_blocks.empty()) {
            throw FormatError("psrio: GUPPI file has no data blocks");
        }
    }

    /// Time-ordered files of one frequency band.
    explicit GuppiReader(std::span<const std::filesystem::path> files) {
        if (files.empty()) {
            throw ValidationError("psrio: GUPPI file list is empty");
        }
        for (const auto& file : files) {
            open_one(file);
        }
        if (m_blocks.empty()) {
            throw FormatError("psrio: GUPPI file has no data blocks");
        }
    }

    /// Open every file. One `OBSFREQ` is a time sequence. Several are stitched
    /// across frequency.
    [[nodiscard]] static GuppiSet open(std::span<const std::filesystem::path> files);

    GuppiReader(const GuppiReader&)                = delete;
    GuppiReader& operator=(const GuppiReader&)     = delete;
    GuppiReader(GuppiReader&&) noexcept            = default;
    GuppiReader& operator=(GuppiReader&&) noexcept = default;
    ~GuppiReader()                                 = default;

    [[nodiscard]] const BasebandHeader& header() const noexcept {
        return m_header;
    }
    [[nodiscard]] std::uint64_t npol() const noexcept { return m_header.npol; }
    [[nodiscard]] std::uint64_t nchan() const noexcept { return m_header.nchan; }
    [[nodiscard]] std::uint64_t nants() const noexcept { return m_header.nants; }
    [[nodiscard]] int nbit() const noexcept { return m_header.nbit; }
    [[nodiscard]] int ndim() const noexcept { return m_header.ndim; }
    [[nodiscard]] SampleType sample_type() const { return m_header.sample_type(); }
    [[nodiscard]] bool samples_signed() const noexcept {
        return m_header.samples_signed;
    }
    [[nodiscard]] bool msb_first() const noexcept {
        return m_header.msb_first;
    }
    void set_msb_first(bool msb) noexcept {
        m_header.msb_first = msb;
    }
    [[nodiscard]] std::uint64_t bytes_per_sample() const noexcept {
        return m_stride;
    }
    [[nodiscard]] std::uint64_t nsamples() const noexcept {
        return m_header.nsamples;
    }
    [[nodiscard]] bool has_nsamples() const noexcept {
        return m_header.has_nsamples;
    }
    [[nodiscard]] double tsamp() const noexcept { return m_header.tsamp; }
    [[nodiscard]] double tstart() const noexcept { return m_header.tstart; }
    [[nodiscard]] std::time_t utc_start() const noexcept {
        return m_header.utc_start;
    }
    [[nodiscard]] double fch1() const noexcept { return m_header.fch1; }
    [[nodiscard]] double foff() const noexcept { return m_header.foff; }
    [[nodiscard]] double bandwidth() const noexcept {
        return m_header.bandwidth;
    }
    [[nodiscard]] double center_frequency() const noexcept {
        return m_header.center_frequency;
    }
    [[nodiscard]] std::string_view source_name() const noexcept {
        return m_header.source;
    }
    [[nodiscard]] std::string_view telescope() const noexcept {
        return m_header.telescope;
    }
    [[nodiscard]] std::uint64_t tell() const noexcept { return m_sample; }

    void seek(std::uint64_t sample) {
        if (!m_header.has_nsamples) {
            throw ValidationError("psrio: baseband length is not known");
        }
        if (sample > m_header.nsamples) {
            throw ValidationError("psrio: seek is past the readable samples");
        }
        m_sample = sample;
    }

    void rewind() { seek(0); }

    void skip(std::int64_t delta) {
        if (!m_header.has_nsamples) {
            throw ValidationError("psrio: baseband length is not known");
        }
        m_sample = detail::apply_skip(m_sample, delta, m_header.nsamples);
    }

    /// On-disk payload of one block, including overlap samples.
    [[nodiscard]] std::span<const std::byte>
    view_native(std::uint64_t block_index) const {
        if (block_index >= m_blocks.size()) {
            throw ValidationError("psrio: GUPPI block index is out of range");
        }
        const auto& block = m_blocks[block_index];
        const auto bytes  = m_files[block.file_index].mapped.bytes();
        return bytes.subspan(static_cast<std::size_t>(block.payload),
                             static_cast<std::size_t>(block.payload_bytes));
    }

    std::uint64_t read_block(std::uint64_t count, std::span<std::byte> dest) {
        detail::require_block(count, dest, m_stride);
        std::uint64_t done = 0;
        while (done < count) {
            const auto left = available();
            if (left == 0U) {
                break;
            }
            const auto& block = find_block(m_sample);
            const auto into   = m_sample - block.first;
            const auto room   = block.valid - into;
            const auto take   = std::min(count - done, room);
            const auto dest_slice = dest.subspan(
                static_cast<std::size_t>(done * m_stride),
                static_cast<std::size_t>(take * m_stride));
            if (block.is_gap) {
                std::memset(dest_slice.data(), 0, dest_slice.size());
            } else {
                const auto native = m_files[block.file_index].mapped.bytes().subspan(
                    static_cast<std::size_t>(block.payload),
                    static_cast<std::size_t>(block.payload_bytes));
                detail::copy_to_canonical(
                    native, block.native_ntime, into, take, m_layout,
                    dest_slice);
            }
            m_sample += take;
            done += take;
        }
        return done;
    }

    std::uint64_t read_bytes(std::uint64_t nbytes, std::span<std::byte> dest) {
        if (dest.size() != nbytes) {
            throw ValidationError(
                "psrio: destination size does not match requested bytes");
        }
        if (m_stride == 0U || nbytes % m_stride != 0U) {
            throw ValidationError(
                "psrio: byte request must be a multiple of the sample stride");
        }
        return read_block(nbytes / m_stride, dest) * m_stride;
    }

    /// Read raw baseband samples in native file layout (FTPRI for GUPPI raw)
    /// without transposing axes to canonical time-major order.
    std::uint64_t read_native(std::uint64_t count, std::span<std::byte> dest) {
        detail::require_block(count, dest, m_stride);
        const auto take_total = std::min(count, available());
        if (take_total == 0U) {
            return 0U;
        }
        const auto pol_dim_bytes =
            m_header.nchan > 0U ? (m_stride / m_header.nchan) : m_stride;
        std::uint64_t done = 0;
        while (done < count) {
            const auto left = available();
            if (left == 0U) {
                break;
            }
            const auto& block = find_block(m_sample);
            const auto into   = m_sample - block.first;
            const auto room   = block.valid - into;
            const auto take   = std::min(count - done, room);

            if (block.is_gap) {
                if (!m_layout.channel_major) {
                    const auto dest_slice = dest.subspan(
                        static_cast<std::size_t>(done * m_stride),
                        static_cast<std::size_t>(take * m_stride));
                    std::memset(dest_slice.data(), 0, dest_slice.size());
                } else {
                    for (std::uint64_t c = 0; c < m_header.nchan; ++c) {
                        auto* dst_c = dest.data() + static_cast<std::size_t>(
                                          ((c * count) + done) * pol_dim_bytes);
                        std::memset(dst_c, 0,
                                    static_cast<std::size_t>(take * pol_dim_bytes));
                    }
                }
            } else {
                const auto native = m_files[block.file_index].mapped.bytes().subspan(
                    static_cast<std::size_t>(block.payload),
                    static_cast<std::size_t>(block.payload_bytes));
                if (!m_layout.channel_major) {
                    const auto dest_slice = dest.subspan(
                        static_cast<std::size_t>(done * m_stride),
                        static_cast<std::size_t>(take * m_stride));
                    std::memcpy(dest_slice.data(),
                                native.data() + static_cast<std::size_t>(into * m_stride),
                                static_cast<std::size_t>(take * m_stride));
                    if (m_layout.offset_binary) {
                        detail::apply_offset_binary(dest_slice);
                    }
                } else {
                    if (m_header.nchan > 0U && m_stride % m_header.nchan == 0U) {
                        for (std::uint64_t c = 0; c < m_header.nchan; ++c) {
                            const auto* src_c =
                                native.data() +
                                static_cast<std::size_t>(
                                    ((c * block.native_ntime) + into) *
                                    pol_dim_bytes);
                            auto* dst_c = dest.data() +
                                          static_cast<std::size_t>(
                                              ((c * count) + done) *
                                              pol_dim_bytes);
                            std::memcpy(
                                dst_c, src_c,
                                static_cast<std::size_t>(take * pol_dim_bytes));
                            if (m_layout.offset_binary) {
                                detail::apply_offset_binary(std::span<std::byte>(
                                    dst_c,
                                    static_cast<std::size_t>(take * pol_dim_bytes)));
                            }
                        }
                    } else {
                        // Channel stride doesn't align to byte; fallback
                        const auto dest_slice = dest.subspan(
                            static_cast<std::size_t>(done * m_stride),
                            static_cast<std::size_t>(take * m_stride));
                        detail::copy_to_canonical(
                            native, block.native_ntime, into, take, m_layout,
                            dest_slice);
                    }
                }
            }
            m_sample += take;
            done += take;
        }
        return done;
    }

    std::uint64_t read_voltages(std::uint64_t nsamps,
                                VoltageRead how,
                                std::span<std::byte> dest) {
        detail::require_dual_pol_bytes(nsamps, dest, m_header, m_stride);
        const auto take_total = std::min(nsamps, available());
        if (take_total == 0U) {
            return 0U;
        }
        if (how.order != VoltageOrder::kFreqMajor || !m_layout.channel_major ||
            m_header.nchan == 0U || m_stride % m_header.nchan != 0U) {
            const auto needed = static_cast<std::size_t>(take_total * m_stride);
            if (m_scratch_bytes.size() < needed) {
                m_scratch_bytes.resize(needed);
            }
            const std::span<std::byte> native(m_scratch_bytes.data(), needed);
            const auto got = read_block(take_total, native);
            detail::arrange_dual_pol(
                native.first(static_cast<std::size_t>(got * m_stride)),
                got, nsamps, m_header.nchan, m_header.nbit,
                how.frequency_ascending && m_header.foff < 0.0,
                how.order == VoltageOrder::kFreqMajor, dest);
            return got;
        }
        const bool reverse_channels =
            how.frequency_ascending && (m_header.foff < 0.0);
        const auto pol_dim_bytes =
            static_cast<std::size_t>(m_stride / m_header.nchan);
        std::uint64_t done = 0;
        while (done < nsamps) {
            const auto left = available();
            if (left == 0U) {
                break;
            }
            const auto& block = find_block(m_sample);
            const auto into   = m_sample - block.first;
            const auto room   = block.valid - into;
            const auto take   = std::min(nsamps - done, room);

            if (block.is_gap) {
                for (std::uint64_t c = 0; c < m_header.nchan; ++c) {
                    const auto dest_c =
                        reverse_channels ? (m_header.nchan - 1U - c) : c;
                    auto* dst_c = dest.data() +
                                  static_cast<std::size_t>(
                                      ((dest_c * nsamps) + done) * pol_dim_bytes);
                    std::memset(dst_c, 0,
                                static_cast<std::size_t>(take * pol_dim_bytes));
                }
            } else {
                const auto native = m_files[block.file_index].mapped.bytes().subspan(
                    static_cast<std::size_t>(block.payload),
                    static_cast<std::size_t>(block.payload_bytes));
                for (std::uint64_t c = 0; c < m_header.nchan; ++c) {
                    const auto dest_c =
                        reverse_channels ? (m_header.nchan - 1U - c) : c;
                    const auto* src_c =
                        native.data() +
                        static_cast<std::size_t>(
                            ((c * block.native_ntime) + into) * pol_dim_bytes);
                    auto* dst_c = dest.data() +
                                  static_cast<std::size_t>(
                                      ((dest_c * nsamps) + done) * pol_dim_bytes);
                    std::memcpy(dst_c, src_c,
                                static_cast<std::size_t>(take * pol_dim_bytes));
                    if (m_layout.offset_binary) {
                        detail::apply_offset_binary(std::span<std::byte>(
                            dst_c, static_cast<std::size_t>(take * pol_dim_bytes)));
                    }
                }
            }
            m_sample += take;
            done += take;
        }
        return done;
    }

    [[nodiscard]] static std::size_t num_groups() noexcept { return 1; }

    std::uint64_t read_groups(std::uint64_t nsamps,
                              std::span<std::span<std::byte>> dest_groups) {
        if (dest_groups.size() != 1U) {
            throw ValidationError("psrio: dest_groups size must match num_groups()");
        }
        const auto take_total = std::min(nsamps, available());
        if (take_total == 0U) {
            return 0U;
        }
        const auto needed = static_cast<std::size_t>(take_total * m_stride);
        if (dest_groups[0].size() < needed) {
            throw ValidationError("psrio: dest_groups buffer is too small");
        }
        if (!m_layout.channel_major || m_header.nchan == 0U ||
            m_stride % m_header.nchan != 0U) {
            return read_voltages(
                nsamps,
                VoltageRead{.order = VoltageOrder::kFreqMajor,
                            .frequency_ascending = true,},
                dest_groups[0]);
        }
        const bool reverse_channels = (m_header.foff < 0.0);
        const auto pol_dim_bytes =
            static_cast<std::size_t>(m_stride / m_header.nchan);
        std::uint64_t done = 0;
        while (done < nsamps) {
            const auto left = available();
            if (left == 0U) {
                break;
            }
            const auto& block = find_block(m_sample);
            const auto into   = m_sample - block.first;
            const auto room   = block.valid - into;
            const auto take   = std::min(nsamps - done, room);

            if (block.is_gap) {
                for (std::uint64_t c = 0; c < m_header.nchan; ++c) {
                    const auto dest_c =
                        reverse_channels ? (m_header.nchan - 1U - c) : c;
                    auto* dst_c = dest_groups[0].data() +
                                  static_cast<std::size_t>(
                                      ((dest_c * nsamps) + done) * pol_dim_bytes);
                    std::memset(dst_c, 0,
                                static_cast<std::size_t>(take * pol_dim_bytes));
                }
            } else {
                const auto native = m_files[block.file_index].mapped.bytes().subspan(
                    static_cast<std::size_t>(block.payload),
                    static_cast<std::size_t>(block.payload_bytes));
                for (std::uint64_t c = 0; c < m_header.nchan; ++c) {
                    const auto dest_c =
                        reverse_channels ? (m_header.nchan - 1U - c) : c;
                    const auto* src_c =
                        native.data() +
                        static_cast<std::size_t>(
                            ((c * block.native_ntime) + into) * pol_dim_bytes);
                    auto* dst_c = dest_groups[0].data() +
                                  static_cast<std::size_t>(
                                      ((dest_c * nsamps) + done) * pol_dim_bytes);
                    std::memcpy(dst_c, src_c,
                                static_cast<std::size_t>(take * pol_dim_bytes));
                    if (m_layout.offset_binary) {
                        detail::apply_offset_binary(std::span<std::byte>(
                            dst_c, static_cast<std::size_t>(take * pol_dim_bytes)));
                    }
                }
            }
            m_sample += take;
            done += take;
        }
        return done;
    }

    [[nodiscard]] std::vector<std::vector<std::byte>> read_groups(std::uint64_t nsamps) {
        std::vector<std::vector<std::byte>> buffers(1);
        buffers[0].resize(static_cast<std::size_t>(nsamps * m_stride));
        std::array<std::span<std::byte>, 1> group_array{buffers[0]};
        const std::span<std::span<std::byte>> spans(group_array.data(), 1);
        read_groups(nsamps, spans);
        return buffers;
    }

    [[nodiscard]] std::uint64_t dropped_packets() const noexcept {
        return m_dropped_packets;
    }

    [[nodiscard]] std::uint64_t dropped_samples() const noexcept {
        return m_dropped_samples;
    }

    template <typename T>
    std::uint64_t read_samples(std::uint64_t count, std::span<T> dest) {
        detail::require_values(count, std::span<const T>(dest), m_header);
        const auto take = std::min(count, available());
        if (take == 0U) {
            return 0U;
        }
        const auto needed = static_cast<std::size_t>(take * m_stride);
        if (m_scratch_bytes.size() < needed) {
            m_scratch_bytes.resize(needed);
        }
        const std::span<std::byte> packed(m_scratch_bytes.data(), needed);
        const auto got = read_block(take, packed);
        const auto values = got * m_header.values_per_sample();
        auto geometry          = m_layout;
        geometry.channel_major = false;
        geometry.pol_major     = false;
        geometry.offset_binary = false;
        if constexpr (kIsBasebandComplex<T>) {
            detail::unpack_complex(
                packed.first(static_cast<std::size_t>(got * m_stride)),
                dest.first(static_cast<std::size_t>(values)), geometry);
        } else {
            detail::unpack_components(
                packed.first(static_cast<std::size_t>(got * m_stride)),
                dest.first(static_cast<std::size_t>(values)), geometry);
        }
        return got;
    }

    template <typename T>
    [[nodiscard]] std::vector<T> read_samples(std::uint64_t count) {
        const auto values = count * m_header.values_per_sample();
        std::vector<T> out(static_cast<std::size_t>(values));
        const auto actual = read_samples(count, std::span<T>(out));
        out.resize(static_cast<std::size_t>(actual * m_header.values_per_sample()));
        return out;
    }

private:
    struct StoredFile {
        detail::MappedFile mapped;
    };

    struct Block {
        std::uint64_t file_index{0};
        std::uint64_t payload{0};
        std::uint64_t payload_bytes{0};
        std::uint64_t native_ntime{0};
        std::uint64_t valid{0};
        std::uint64_t first{0};
        std::int64_t pktidx{0};
        bool has_pktidx{false};
        bool is_gap{false};
    };

    void open_one(const std::filesystem::path& path) {
        StoredFile stored;
        try {
            stored.mapped = detail::MappedFile(path);
            stored.mapped.advise_sequential();
        } catch (const std::exception& ex) {
            throw IoError(ex.what());
        }
        m_files.push_back(std::move(stored));
        if (m_header.filename.empty()) {
            m_header.filename = path.string();
        }
        index_file(m_files.size() - 1U);
    }

    void index_file(std::uint64_t file_index) {
        const auto bytes = m_files[file_index].mapped.bytes();
        std::uint64_t pos = 0;
        while (pos < bytes.size()) {
            if (bytes.size() - pos < kRecordBytes) {
                throw FormatError(
                    "psrio: GUPPI trailing bytes are not a header");
            }
            GuppiHeader header;
            const auto consumed =
                header.parse(bytes.subspan(static_cast<std::size_t>(pos)));
            if (!header.has_key("BLOCSIZE")) {
                throw FormatError("psrio: GUPPI block is missing BLOCSIZE");
            }
            const auto bloc = header.int_or("BLOCSIZE", std::int64_t{-1});
            if (bloc < 0) {
                throw ValidationError("psrio: GUPPI BLOCSIZE must be >= 0");
            }
            const auto payload = pos + consumed;
            if (static_cast<std::uint64_t>(bloc) > bytes.size() - payload) {
                throw FormatError(
                    "psrio: GUPPI payload extends past the end of the file");
            }
            adopt(header, file_index, payload, static_cast<std::uint64_t>(bloc));
            pos = payload + static_cast<std::uint64_t>(bloc);
        }
    }

    void adopt(const GuppiHeader& header,
               std::uint64_t file_index,
               std::uint64_t payload,
               std::uint64_t bloc) {
        const auto native = header.ntime();
        const auto overlap = header.int_or("OVERLAP", std::int64_t{0});
        if (native <= 0 || overlap < 0 || overlap >= native) {
            throw ValidationError("psrio: GUPPI overlap does not fit in the block");
        }
        if (m_blocks.empty()) {
            const auto filename = m_header.filename;
            m_proto              = header;
            m_header             = header.to_baseband_header();
            m_header.filename    = filename;
            m_header.nsamples    = 0;
            m_header.has_nsamples = true;
            m_layout             = layout();
            m_layout.offset_binary =
                header.has_key("PKTFMT") &&
                equals_ignore_case(header.string_or("PKTFMT", ""), "VDIF");
            m_stride = m_header.bytes_per_sample();
        } else if (!same_band(m_proto, header)) {
            throw ValidationError(
                "psrio: GUPPI blocks do not share a geometry");
        }
        const auto native_u = static_cast<std::uint64_t>(native);
        if (m_stride == 0U || bloc != native_u * m_stride) {
            throw ValidationError("psrio: GUPPI BLOCSIZE does not match the geometry");
        }
        Block block;
        block.file_index    = file_index;
        block.payload       = payload;
        block.payload_bytes = bloc;
        block.native_ntime  = native_u;
        block.valid         = native_u - static_cast<std::uint64_t>(overlap);
        block.first         = m_header.nsamples;
        block.has_pktidx    = header.has_key("PKTIDX");
        if (block.has_pktidx) {
            block.pktidx = header.int_or("PKTIDX", std::int64_t{0});
        }
        if (!m_blocks.empty()) {
            if (!block.has_pktidx || !m_blocks.back().has_pktidx) {
                throw FormatError(
                    "psrio: GUPPI PKTIDX is required across blocks");
            }
            if (block.pktidx < m_blocks.back().pktidx) {
                throw FormatError("psrio: GUPPI PKTIDX moved backward");
            }
            const auto step = static_cast<std::int64_t>(
                packet_step(header, block.valid * m_stride, block.valid));
            const auto expected = m_blocks.back().pktidx + step;
            if (block.pktidx > expected) {
                const auto diff = block.pktidx - expected;
                static constexpr std::uint64_t kMaxPacketJump = 10'000'000;
                if (std::cmp_greater(diff, kMaxPacketJump)) {
                    throw FormatError(std::format(
                        "psrio: GUPPI PKTIDX jump from {} to {} exceeds maximum allowed gap ({})",
                        m_blocks.back().pktidx, block.pktidx, kMaxPacketJump));
                }
                std::uint64_t dropped_samps = 0;
                if (step > 0 && block.valid > 0 &&
                    block.valid % static_cast<std::uint64_t>(step) == 0U) {
                    dropped_samps =
                        static_cast<std::uint64_t>(diff) *
                        (block.valid / static_cast<std::uint64_t>(step));
                } else {
                    const auto width = packet_width(header);
                    if (width > 0 && m_stride > 0) {
                        dropped_samps =
                            (static_cast<std::uint64_t>(diff) * width) / m_stride;
                    } else {
                        dropped_samps = static_cast<std::uint64_t>(diff);
                    }
                }
                m_dropped_packets += static_cast<std::uint64_t>(diff);
                m_dropped_samples += dropped_samps;

                Block gap_block;
                gap_block.first      = m_header.nsamples;
                gap_block.valid      = dropped_samps;
                gap_block.is_gap     = true;
                gap_block.has_pktidx = true;
                gap_block.pktidx     = expected;
                m_header.nsamples += dropped_samps;
                m_blocks.push_back(gap_block);
            } else if (block.pktidx < expected) {
                throw FormatError(std::format(
                    "psrio: GUPPI PKTIDX jumped from {} to {}",
                    m_blocks.back().pktidx, block.pktidx));
            }
        }
        block.first = m_header.nsamples;
        m_header.nsamples += block.valid;
        m_blocks.push_back(block);
    }

    [[nodiscard]] std::uint64_t available() const {
        if (!m_header.has_nsamples) {
            throw ValidationError("psrio: baseband length is not known");
        }
        if (m_sample > m_header.nsamples) {
            throw ValidationError(
                "psrio: reader cursor is past the readable samples");
        }
        return m_header.nsamples - m_sample;
    }

    [[nodiscard]] detail::BasebandLayout layout() const {
        detail::BasebandLayout out;
        out.nants          = m_header.nants;
        out.nchan          = m_header.nchan;
        out.npol           = m_header.npol;
        out.nbit           = m_header.nbit;
        out.ndim           = m_header.ndim;
        out.samples_signed = m_header.samples_signed;
        out.msb_first      = m_header.msb_first;
        out.channel_major  = m_header.order == BasebandOrder::kChannelMajor;
        return out;
    }

    [[nodiscard]] const Block& find_block(std::uint64_t sample) const {
        std::size_t lower = 0;
        std::size_t upper = m_blocks.size();
        while (lower + 1U < upper) {
            const auto mid = lower + ((upper - lower) / 2U);
            if (m_blocks[mid].first <= sample) {
                lower = mid;
            } else {
                upper = mid;
            }
        }
        return m_blocks[lower];
    }

    [[nodiscard]] static bool same_band(const GuppiHeader& first,
                                        const GuppiHeader& other) {
        if (first.int_or("OBSNCHAN", 0) != other.int_or("OBSNCHAN", 0) ||
            first.int_or("NPOL", 1) != other.int_or("NPOL", 1) ||
            first.int_or("NBITS", 8) != other.int_or("NBITS", 8) ||
            first.int_or("NANTS", 1) != other.int_or("NANTS", 1) ||
            first.int_or("BLOCSIZE", 0) != other.int_or("BLOCSIZE", 0) ||
            first.int_or("OVERLAP", 0) != other.int_or("OVERLAP", 0)) {
            return false;
        }
        if (first.real_or("TBIN", 0.0) != other.real_or("TBIN", 0.0) ||
            first.real_or("OBSFREQ", 0.0) != other.real_or("OBSFREQ", 0.0) ||
            first.channel_spacing() != other.channel_spacing()) {
            return false;
        }
        return equals_ignore_case(first.string_or("PKTFMT", ""),
                                  other.string_or("PKTFMT", ""));
    }

    [[nodiscard]] static std::uint64_t packet_step(const GuppiHeader& header,
                                                   std::uint64_t valid_bytes,
                                                   std::uint64_t valid_samples) {
        if (!header.has_key("PKTSIZE")) {
            return valid_samples;
        }
        auto packet = header.int_or("PKTSIZE", std::int64_t{0});
        if (packet <= 0) {
            throw ValidationError("psrio: GUPPI PKTSIZE must be positive");
        }
        const auto format = header.string_or("PKTFMT", "");
        if (equals_ignore_case(format, "VDIF") && header.complex_npol() == 2 &&
            packet <= (std::numeric_limits<std::int64_t>::max() / 2)) {
            packet *= 2;
        }
        if (equals_ignore_case(format, "1SFA") &&
            header.int_or("NBITS", std::int64_t{8}) == 2) {
            if (packet % 4 != 0) {
                throw ValidationError(
                    "psrio: GUPPI 2-bit PKTSIZE is not divisible by 4");
            }
            packet /= 4;
        }
        const auto width = static_cast<std::uint64_t>(packet);
        if (width == 0U || valid_bytes % width != 0U) {
            throw ValidationError(
                "psrio: GUPPI valid block is not a whole number of packets");
        }
        return valid_bytes / width;
    }

    [[nodiscard]] static std::uint64_t packet_width(const GuppiHeader& header) {
        if (!header.has_key("PKTSIZE")) {
            return 0;
        }
        auto packet = header.int_or("PKTSIZE", std::int64_t{0});
        if (packet <= 0) {
            return 0;
        }
        const auto format = header.string_or("PKTFMT", "");
        if (equals_ignore_case(format, "VDIF") && header.complex_npol() == 2 &&
            packet <= (std::numeric_limits<std::int64_t>::max() / 2)) {
            packet *= 2;
        }
        if (equals_ignore_case(format, "1SFA") &&
            header.int_or("NBITS", std::int64_t{8}) == 2) {
            packet /= 4;
        }
        return static_cast<std::uint64_t>(packet);
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

    BasebandHeader m_header;
    std::uint64_t m_sample{0};
    std::uint64_t m_stride{0};
    std::vector<std::byte> m_scratch_bytes;
    std::vector<StoredFile> m_files;
    std::vector<Block> m_blocks;
    GuppiHeader m_proto;
    detail::BasebandLayout m_layout;
    std::uint64_t m_dropped_packets{0};
    std::uint64_t m_dropped_samples{0};
};

/// Type-erased GUPPI set: one band, or several bands stitched in frequency.
class GuppiSet : public BasebandSource {
public:
    [[nodiscard]] static GuppiSet open(std::span<const std::filesystem::path> files);

private:
    explicit GuppiSet(BasebandSource source) : BasebandSource(std::move(source)) {}
};

inline GuppiSet GuppiReader::open(std::span<const std::filesystem::path> files) {
    return GuppiSet::open(files);
}

inline GuppiSet GuppiSet::open(std::span<const std::filesystem::path> files) {
    if (files.empty()) {
        throw ValidationError("psrio: GUPPI file list is empty");
    }
    struct Group {
        double frequency{0.0};
        std::vector<std::filesystem::path> paths;
        std::vector<std::int64_t> packet_index;
    };
    std::vector<Group> groups;
    for (const auto& file : files) {
        RawReader raw(file);
        const auto header = raw.read_header();
        if (!header.has_key("OBSFREQ")) {
            throw FormatError("psrio: GUPPI header missing OBSFREQ");
        }
        const auto frequency = header.real_or("OBSFREQ", 0.0);
        const auto packet    = header.int_or("PKTIDX", std::int64_t{0});
        const auto found = std::ranges::find_if(groups, [&](const Group& group) {
            return group.frequency == frequency;
        });
        if (found == groups.end()) {
            groups.push_back(Group{.frequency = frequency, .paths = {file}, .packet_index = {packet}});
        } else {
            found->paths.push_back(file);
            found->packet_index.push_back(packet);
        }
    }
    std::vector<BasebandSource> bands;
    bands.reserve(groups.size());
    for (auto& group : groups) {
        std::vector<std::size_t> order(group.paths.size());
        for (std::size_t index = 0; index < order.size(); ++index) {
            order[index] = index;
        }
        std::ranges::stable_sort(order, [&](std::size_t left, std::size_t right) {
            return group.packet_index[left] < group.packet_index[right];
        });
        std::vector<std::filesystem::path> ordered;
        ordered.reserve(order.size());
        for (const auto index : order) {
            ordered.push_back(std::move(group.paths[index]));
        }
        bands.emplace_back(GuppiReader(std::span<const std::filesystem::path>(
            ordered.data(), ordered.size())));
    }
    if (bands.size() == 1U) {
        return GuppiSet{std::move(bands.front())};
    }
    return GuppiSet{BasebandSource{FrequencyStitch{std::move(bands)}}};
}

static_assert(concepts::BasebandReader<GuppiReader>);

} // namespace psrio::formats::guppi
