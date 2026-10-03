#pragma once

/// PSRDADA file reader. One ASCII header of `HDR_SIZE` bytes and one payload.
/// Included from `psrio.hpp`. The live ring stays in `formats/psrdada.hpp`.

#include "psrio/baseband.hpp"
#include "psrio/detail/mmap.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <exception>
#include <filesystem>
#include <format>
#include <map>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace psrio::formats::dada {

inline constexpr std::uint64_t kDefaultHeaderBytes = 4096;

using HeaderValue = std::variant<std::int64_t, double, std::string>;

/// One DADA text header.
class DadaHeader {
public:
    /// Parse a header at the front of @p bytes. Returns `HDR_SIZE`.
    std::uint64_t parse(std::span<const std::byte> bytes) {
        m_entries.clear();
        std::uint64_t hdr_size = kDefaultHeaderBytes;
        std::uint64_t pos      = 0;
        bool saw_end           = false;
        while (pos < hdr_size && pos < bytes.size()) {
            if (bytes[static_cast<std::size_t>(pos)] == std::byte{0}) {
                break;
            }
            std::uint64_t end = pos;
            while (end < bytes.size() && end < hdr_size &&
                   bytes[static_cast<std::size_t>(end)] != std::byte{'\n'} &&
                   bytes[static_cast<std::size_t>(end)] != std::byte{0}) {
                ++end;
            }
            if (end >= hdr_size && end < bytes.size() &&
                bytes[static_cast<std::size_t>(end)] != std::byte{'\n'} &&
                bytes[static_cast<std::size_t>(end)] != std::byte{0}) {
                throw FormatError("psrio: DADA header line crosses HDR_SIZE");
            }
            if (end == bytes.size() && end < hdr_size) {
                throw FormatError("psrio: DADA header ended before HDR_SIZE");
            }
            std::string line;
            line.resize(static_cast<std::size_t>(end - pos));
            for (std::uint64_t index = pos; index < end; ++index) {
                line[static_cast<std::size_t>(index - pos)] = static_cast<char>(
                    std::to_integer<unsigned char>(bytes[static_cast<std::size_t>(index)]));
            }
            if (end < bytes.size() &&
                bytes[static_cast<std::size_t>(end)] == std::byte{'\n'}) {
                pos = end + 1;
            } else {
                pos = end;
            }
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            if (line.empty()) {
                continue;
            }
            const auto comment = line.find('#');
            const auto body =
                comment == std::string::npos ? std::string_view{line}
                                             : std::string_view{line}.substr(0, comment);
            if (comment != std::string::npos &&
                line.find("end of header") != std::string::npos) {
                saw_end = true;
                break;
            }
            const auto trimmed = trim(body);
            if (trimmed.empty()) {
                continue;
            }
            const auto gap = trimmed.find_first_of(" \t");
            const auto key = upper(gap == std::string_view::npos
                                       ? trimmed
                                       : trimmed.substr(0, gap));
            const auto value =
                gap == std::string_view::npos ? std::string_view{} : trim(trimmed.substr(gap + 1));
            store(key, value);
            if (key == "HDR_SIZE") {
                const auto size = int_or("HDR_SIZE", std::int64_t{0});
                if (size <= 0) {
                    throw FormatError("psrio: DADA HDR_SIZE must be positive");
                }
                hdr_size = static_cast<std::uint64_t>(size);
            }
        }
        if (bytes.size() < hdr_size) {
            throw FormatError("psrio: DADA header extends past the end of the file");
        }
        if (!saw_end && pos < hdr_size &&
            bytes[static_cast<std::size_t>(pos)] != std::byte{0} &&
            m_entries.empty()) {
            throw FormatError("psrio: DADA header is empty");
        }
        if (!has_key("HEADER") || !equals_ignore_case(string_or("HEADER", ""), "DADA")) {
            throw FormatError("psrio: DADA header is missing HEADER DADA");
        }
        if (has_key("ORDER") && !equals_ignore_case(string_or("ORDER", ""), "TF")) {
            throw FormatError("psrio: DADA ORDER must be TF");
        }
        m_hdr_size = hdr_size;
        return hdr_size;
    }

    [[nodiscard]] bool has_key(std::string_view key) const {
        return m_entries.contains(upper(key));
    }

    [[nodiscard]] std::int64_t int_or(std::string_view key,
                                      std::int64_t fallback) const {
        const auto it = m_entries.find(upper(key));
        if (it == m_entries.end()) {
            return fallback;
        }
        if (!std::holds_alternative<std::int64_t>(it->second)) {
            throw FormatError(std::format(
                "psrio: DADA key {} does not hold an integer", upper(key)));
        }
        return std::get<std::int64_t>(it->second);
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
            "psrio: DADA key {} does not hold a number", upper(key)));
    }

    [[nodiscard]] std::string string_or(std::string_view key,
                                        std::string fallback) const {
        const auto it = m_entries.find(upper(key));
        if (it == m_entries.end()) {
            return fallback;
        }
        if (!std::holds_alternative<std::string>(it->second)) {
            throw FormatError(std::format(
                "psrio: DADA key {} does not hold a string", upper(key)));
        }
        return std::get<std::string>(it->second);
    }

    [[nodiscard]] std::uint64_t header_bytes() const noexcept { return m_hdr_size; }

    [[nodiscard]] BasebandHeader to_baseband_header() const {
        if (!has_key("NBIT") || !has_key("NPOL") || !has_key("TSAMP")) {
            throw FormatError(
                "psrio: DADA header needs NBIT, NPOL, and TSAMP");
        }
        BasebandHeader out;
        out.format = BasebandFormat::kDada;
        out.order  = BasebandOrder::kTimeMajor;
        out.npol   = static_cast<std::uint64_t>(require_positive("NPOL"));
        out.nchan  = static_cast<std::uint64_t>(
            has_key("NCHAN") ? require_positive("NCHAN") : 1);
        out.nants  = 1;
        out.nbit   = static_cast<int>(int_or("NBIT", std::int64_t{8}));
        out.ndim   = static_cast<int>(int_or("NDIM", std::int64_t{1}));
        if (out.ndim != 1 && out.ndim != 2) {
            throw ValidationError("psrio: DADA NDIM must be 1 or 2");
        }
        out.samples_signed = out.nbit < 0 || out.ndim == 2;
        out.tsamp          = real_or("TSAMP", 0.0) * 1.0e-6;
        const auto bandwidth = has_key("BW") ? real_or("BW", 0.0) : 0.0;
        const auto center    = has_key("FREQ") ? real_or("FREQ", 0.0) : 0.0;
        out.foff = out.nchan == 0U ? 0.0 : bandwidth / static_cast<double>(out.nchan);
        out.fch1 = center - (bandwidth / 2.0) + (out.foff / 2.0);
        out.bandwidth        = std::abs(bandwidth);
        out.center_frequency = center;
        out.source           = string_or("SOURCE", "Unknown");
        out.telescope        = string_or("TELESCOPE", "Unknown");
        out.backend          = string_or("INSTRUMENT", "");
        if (has_key("MJD_START")) {
            try {
                out.tstart = std::stod(string_or("MJD_START", "0"));
            } catch (const std::exception&) {
                throw FormatError("psrio: DADA MJD_START is not a number");
            }
        }
        if (has_key("UTC_START")) {
            out.utc_start = parse_utc(string_or("UTC_START", ""));
            if (!has_key("MJD_START")) {
                constexpr double kUnixEpochMjd = 40587.0;
                constexpr double kSecondsInDay = 86400.0;
                out.tstart = kUnixEpochMjd +
                             (static_cast<double>(out.utc_start) / kSecondsInDay);
            }
        } else if (out.tstart != 0.0) {
            out.utc_start = astro::mjd_to_time(out.tstart);
        }
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
    std::uint64_t m_hdr_size{kDefaultHeaderBytes};

    [[nodiscard]] std::int64_t require_positive(const char* key) const {
        const auto value = int_or(key, std::int64_t{0});
        if (value <= 0) {
            throw ValidationError(
                std::format("psrio: DADA {} must be positive", key));
        }
        return value;
    }

    void store(const std::string& key, std::string_view value) {
        if (is_int_key(key)) {
            try {
                m_entries[key] = static_cast<std::int64_t>(std::stoll(std::string{value}));
            } catch (const std::exception&) {
                throw FormatError(
                    std::format("psrio: DADA key {} is not an integer", key));
            }
            return;
        }
        if (is_real_key(key)) {
            try {
                m_entries[key] = std::stod(std::string{value});
            } catch (const std::exception&) {
                throw FormatError(
                    std::format("psrio: DADA key {} is not a number", key));
            }
            return;
        }
        m_entries[key] = std::string{value};
    }

    [[nodiscard]] static bool is_int_key(std::string_view key) {
        return key == "HDR_SIZE" || key == "FILE_SIZE" || key == "FILE_NUMBER" ||
               key == "OBS_OFFSET" || key == "OBS_OVERLAP" || key == "NBIT" ||
               key == "NDIM" || key == "NPOL" || key == "NCHAN" ||
               key == "RESOLUTION" || key == "DSB";
    }

    [[nodiscard]] static bool is_real_key(std::string_view key) {
        return key == "FREQ" || key == "BW" || key == "TSAMP";
    }

    [[nodiscard]] static std::string upper(std::string_view key) {
        std::string out{key};
        for (char& character : out) {
            character = static_cast<char>(
                std::toupper(static_cast<unsigned char>(character)));
        }
        return out;
    }

    [[nodiscard]] static std::string_view trim(std::string_view text) {
        const auto first = text.find_first_not_of(" \t");
        if (first == std::string_view::npos) {
            return {};
        }
        const auto last = text.find_last_not_of(" \t");
        return text.substr(first, last - first + 1);
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

    [[nodiscard]] static std::time_t parse_utc(const std::string& text) {
        std::tm broken{};
        if (strptime(text.c_str(), "%Y-%m-%d-%H:%M:%S", &broken) == nullptr) {
            throw FormatError(
                "psrio: DADA UTC_START must be YYYY-MM-DD-HH:MM:SS");
        }
        const auto utc = timegm(&broken);
        if (utc == static_cast<std::time_t>(-1)) {
            throw FormatError("psrio: DADA UTC_START is not a valid time");
        }
        return utc;
    }
};

/// Memory-mapped DADA file or a time sequence of files ordered by `OBS_OFFSET`.
class DadaReader {
public:
    explicit DadaReader(const std::filesystem::path& path) {
        open_files(std::span<const std::filesystem::path>(&path, 1));
    }

    explicit DadaReader(std::span<const std::filesystem::path> files) {
        open_files(files);
    }

    [[nodiscard]] static DadaReader open(std::span<const std::filesystem::path> files) {
        return DadaReader(files);
    }

    DadaReader(const DadaReader&)                = delete;
    DadaReader& operator=(const DadaReader&)     = delete;
    DadaReader(DadaReader&&) noexcept            = default;
    DadaReader& operator=(DadaReader&&) noexcept = default;
    ~DadaReader()                                = default;

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

    /// Canonical bytes when they are contiguous in one file.
    [[nodiscard]] std::span<const std::byte> view_bytes(std::uint64_t nbytes) {
        if (m_layout.pol_major) {
            throw ValidationError(
                "psrio: DADA view is not contiguous in canonical order");
        }
        if (m_stride == 0U || nbytes % m_stride != 0U) {
            throw ValidationError(
                "psrio: byte request must be a multiple of the sample stride");
        }
        const auto count = nbytes / m_stride;
        if (count > available()) {
            throw ValidationError("psrio: DADA view extends past the file");
        }
        const auto& segment = find_segment(m_sample);
        const auto into     = m_sample - segment.first;
        if (count > segment.nsamples - into) {
            throw ValidationError("psrio: DADA view crosses a file boundary");
        }
        const auto bytes = segment.mapped.bytes().subspan(
            static_cast<std::size_t>(segment.payload + (into * m_stride)),
            static_cast<std::size_t>(nbytes));
        m_sample += count;
        return bytes;
    }

    std::uint64_t read_block(std::uint64_t count, std::span<std::byte> dest) {
        detail::require_block(count, dest, m_stride);
        std::uint64_t done = 0;
        while (done < count) {
            const auto left = available();
            if (left == 0U) {
                break;
            }
            const auto& segment = find_segment(m_sample);
            const auto into     = m_sample - segment.first;
            const auto room     = segment.nsamples - into;
            const auto take     = std::min(count - done, room);
            const auto native   = segment.mapped.bytes().subspan(
                static_cast<std::size_t>(segment.payload),
                static_cast<std::size_t>(segment.payload_bytes));
            detail::copy_to_canonical(
                native, segment.nsamples, into, take, m_layout,
                dest.subspan(static_cast<std::size_t>(done * m_stride),
                             static_cast<std::size_t>(take * m_stride)));
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
        const auto got    = read_block(take, packed);
        const auto values = got * m_header.values_per_sample();
        auto geometry          = m_layout;
        geometry.pol_major     = false;
        geometry.channel_major = false;
        if constexpr (kIsBasebandComplex<T>) {
            detail::unpack_complex(
                packed.first(static_cast<std::size_t>(got * m_stride)),
                dest.first(static_cast<std::size_t>(values)),
                geometry);
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

    std::uint64_t read_voltages(std::uint64_t nsamps,
                                VoltageRead how,
                                std::span<std::byte> dest) {
        detail::require_dual_pol_bytes(nsamps, dest, m_header, m_stride);
        const auto take = std::min(nsamps, available());
        if (take == 0U) {
            return 0U;
        }
        const auto needed = static_cast<std::size_t>(take * m_stride);
        if (m_scratch_bytes.size() < needed) {
            m_scratch_bytes.resize(needed);
        }
        const std::span<std::byte> native(m_scratch_bytes.data(), needed);
        const auto got = read_block(take, native);
        detail::arrange_dual_pol(
            native.first(static_cast<std::size_t>(got * m_stride)),
            got, nsamps, m_header.nchan, m_header.nbit,
            how.frequency_ascending && m_header.foff < 0.0,
            how.order == VoltageOrder::kFreqMajor, dest);
        return got;
    }

    [[nodiscard]] static std::size_t num_groups() noexcept {
        return 1;
    }

    std::uint64_t read_native(std::uint64_t count, std::span<std::byte> dest) {
        return read_block(count, dest);
    }

    std::uint64_t read_groups(std::uint64_t count,
                              std::span<std::span<std::byte>> dest_groups) {
        if (dest_groups.size() != 1U) {
            throw ValidationError("psrio: dest_groups size must match num_groups()");
        }
        return read_voltages(
            count,
            VoltageRead{.order = VoltageOrder::kFreqMajor,
                        .frequency_ascending = true,},
            dest_groups[0]);
    }

    [[nodiscard]] std::vector<std::vector<std::byte>> read_groups(std::uint64_t count) {
        std::vector<std::vector<std::byte>> buffers(1);
        buffers[0].resize(static_cast<std::size_t>(count * m_stride));
        std::array<std::span<std::byte>, 1> group_array{buffers[0]};
        const std::span<std::span<std::byte>> spans(group_array.data(), 1);
        read_groups(count, spans);
        return buffers;
    }

    [[nodiscard]] static std::uint64_t dropped_packets() noexcept { return 0; }
    [[nodiscard]] static std::uint64_t dropped_samples() noexcept { return 0; }

private:
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
        out.nants           = m_header.nants;
        out.nchan           = m_header.nchan;
        out.npol            = m_header.npol;
        out.nbit            = m_header.nbit;
        out.ndim            = m_header.ndim;
        out.samples_signed  = m_header.samples_signed;
        out.msb_first       = m_header.msb_first;
        out.channel_major   = m_header.order == BasebandOrder::kChannelMajor;
        return out;
    }

    struct Segment {
        detail::MappedFile mapped;
        std::uint64_t payload{0};
        std::uint64_t payload_bytes{0};
        std::uint64_t nsamples{0};
        std::uint64_t first{0};
        std::uint64_t obs_offset{0};
        std::uint64_t raw_payload{0};
    };

    void open_files(std::span<const std::filesystem::path> files) {
        if (files.empty()) {
            throw ValidationError("psrio: DADA file list is empty");
        }
        struct Item {
            std::filesystem::path path;
            std::uint64_t offset{0};
        };
        std::vector<Item> items;
        items.reserve(files.size());
        for (const auto& file : files) {
            const detail::MappedFile mapped(file);
            DadaHeader header;
            header.parse(mapped.bytes());
            const auto offset = header.int_or("OBS_OFFSET", std::int64_t{0});
            if (offset < 0) {
                throw ValidationError("psrio: DADA OBS_OFFSET must be >= 0");
            }
            items.push_back(Item{.path = file, .offset = static_cast<std::uint64_t>(offset)});
        }
        std::ranges::stable_sort(items,
                                 [](const Item& left, const Item& right) {
                                     return left.offset < right.offset;
                                 });
        for (const auto& item : items) {
            open_one(item.path);
        }
    }

    void open_one(const std::filesystem::path& path) {
        Segment segment;
        try {
            segment.mapped = detail::MappedFile(path);
            segment.mapped.advise_sequential();
        } catch (const std::exception& ex) {
            throw IoError(ex.what());
        }
        DadaHeader header;
        const auto hdr_bytes = header.parse(segment.mapped.bytes());
        const auto filename  = m_header.filename.empty() ? path.string() : m_header.filename;
        auto built           = header.to_baseband_header();
        built.filename       = filename;
        if (m_segments.empty()) {
            m_header             = std::move(built);
            m_header.nsamples    = 0;
            m_header.has_nsamples = true;
            m_header.overlap     = 0;
            m_stride             = m_header.bytes_per_sample();
            m_layout             = layout();
            m_layout.pol_major   = m_header.npol > 1U && m_header.nchan > 1U;
            m_layout.channel_major = false;
            const auto resolution = header.int_or("RESOLUTION", std::int64_t{1});
            if (resolution > 1 &&
                m_stride % static_cast<std::uint64_t>(resolution) != 0U) {
                throw ValidationError(
                    "psrio: DADA sample stride is not a multiple of RESOLUTION");
            }
        } else if (built.npol != m_header.npol || built.nchan != m_header.nchan ||
                   built.nbit != m_header.nbit || built.ndim != m_header.ndim ||
                   built.tsamp != m_header.tsamp) {
            throw ValidationError("psrio: DADA files do not share a geometry");
        }
        const auto file_size = segment.mapped.bytes().size();
        if (hdr_bytes > file_size) {
            throw FormatError("psrio: DADA header extends past the end of the file");
        }
        const auto rest = file_size - hdr_bytes;
        const auto declared = header.int_or("FILE_SIZE", std::int64_t{0});
        if (declared < 0) {
            throw ValidationError("psrio: DADA FILE_SIZE must be >= 0");
        }
        std::uint64_t raw = rest;
        if (declared > 0) {
            const auto dec = static_cast<std::uint64_t>(declared);
            if (dec > file_size) {
                throw FormatError("psrio: DADA FILE_SIZE extends past the end of the file");
            }
            if (dec == file_size || dec > rest) {
                raw = dec - hdr_bytes;
            } else if (dec == rest) {
                raw = dec;
            } else {
                if (dec >= hdr_bytes && (dec - hdr_bytes) % m_stride == 0 &&
                    dec % m_stride != 0) {
                    raw = dec - hdr_bytes;
                } else {
                    raw = dec;
                }
            }
        }
        const auto overlap = header.int_or("OBS_OVERLAP", std::int64_t{0});
        if (overlap < 0) {
            throw ValidationError("psrio: DADA OBS_OVERLAP must be >= 0");
        }
        const auto skip =
            m_segments.empty() ? std::uint64_t{0}
                               : static_cast<std::uint64_t>(overlap);
        if (skip > raw) {
            throw ValidationError("psrio: DADA OBS_OVERLAP exceeds the payload");
        }
        const auto obs_offset = header.int_or("OBS_OFFSET", std::int64_t{0});
        if (obs_offset < 0) {
            throw ValidationError("psrio: DADA OBS_OFFSET must be >= 0");
        }
        if (!m_segments.empty()) {
            const auto& previous = m_segments.back();
            const auto previous_end = previous.obs_offset + previous.raw_payload;
            if (static_cast<std::uint64_t>(obs_offset) + skip != previous_end) {
                throw FormatError("psrio: DADA files do not meet in OBS_OFFSET");
            }
        }
        const auto valid = raw - skip;
        if (m_stride == 0U || valid % m_stride != 0U) {
            throw ValidationError(
                "psrio: DADA payload is not a whole number of samples");
        }
        segment.payload       = hdr_bytes + skip;
        segment.payload_bytes = valid;
        segment.nsamples      = valid / m_stride;
        segment.first         = m_header.nsamples;
        segment.obs_offset    = static_cast<std::uint64_t>(obs_offset);
        segment.raw_payload   = raw;
        m_header.nsamples += segment.nsamples;
        m_segments.push_back(std::move(segment));
    }

    [[nodiscard]] const Segment& find_segment(std::uint64_t sample) const {
        std::size_t lower = 0;
        std::size_t upper = m_segments.size();
        while (lower + 1U < upper) {
            const auto mid = lower + ((upper - lower) / 2U);
            if (m_segments[mid].first <= sample) {
                lower = mid;
            } else {
                upper = mid;
            }
        }
        return m_segments[lower];
    }

    BasebandHeader m_header;
    std::uint64_t m_sample{0};
    std::uint64_t m_stride{0};
    mutable std::vector<std::byte> m_scratch_bytes;
    std::vector<Segment> m_segments;
    detail::BasebandLayout m_layout;
};

static_assert(concepts::BasebandReader<DadaReader>);

} // namespace psrio::formats::dada
