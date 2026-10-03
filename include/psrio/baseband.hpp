#pragma once

/// Baseband block front. Intensity readers stay on `BlockSource`.
/// GUPPI RAW and DADA model `concepts::BasebandReader` and can be stored in
/// `BasebandSource`. `FrequencyStitch` presents several bands as one channel
/// axis.

#include "psrio/common/types.hpp"
#include "psrio/detail/baseband_layout.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/detail/skip.hpp"
#include "psrio/header.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <limits>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace psrio {

enum class BasebandFormat : std::uint8_t { kGuppiRaw, kDada, kMemory };

/// On-disk axis order. Readers present canonical time-major samples.
enum class BasebandOrder : std::uint8_t { kTimeMajor, kChannelMajor };

/// Axis order for a dual-polarisation export. Outermost axis first.
enum class VoltageOrder : std::uint8_t {
    kTimeMajor, ///< `[time][channel][pol][real, imag]` (`TFPRI`)
    kFreqMajor, ///< `[channel][time][pol][real, imag]` (`FTPRI`)
};

/// How `read_voltages` lays out one coherent-dedispersion block.
struct VoltageRead {
    VoltageOrder order{VoltageOrder::kFreqMajor};
    /// Channel 0 is the lowest frequency. A negative `foff` reverses each band.
    bool frequency_ascending{true};
};

/**
 * @brief Flat metadata for one baseband voltage stream.
 *
 * This is not `psrio::Header` and it is not an intensity block. Format cards
 * that do not have a field are kept in `extra`.
 */
struct BasebandHeader {
    BasebandFormat format{BasebandFormat::kMemory};
    BasebandOrder order{BasebandOrder::kTimeMajor};
    std::uint64_t npol{1};
    std::uint64_t nchan{1};
    std::uint64_t nants{1};
    int nbit{8};
    int ndim{2};
    bool samples_signed{true};
    bool msb_first{false};
    double tsamp{0.0};
    double tstart{0.0};
    std::time_t utc_start{0};
    double fch1{0.0};
    double foff{0.0};
    double center_frequency{0.0};
    double bandwidth{0.0};
    std::string source{"Unknown"};
    std::string telescope{"Unknown"};
    std::string backend;
    std::string filename;
    std::uint64_t overlap{0};
    std::uint64_t nsamples{0};
    bool has_nsamples{false};
    std::vector<ExtraKey> extra;

    /// One time sample: every antenna, channel, polarization, and component.
    [[nodiscard]] std::uint64_t bytes_per_sample() const {
        detail::BasebandLayout layout;
        layout.nants = nants;
        layout.nchan = nchan;
        layout.npol  = npol;
        layout.nbit  = nbit;
        layout.ndim  = ndim;
        return detail::sample_stride(layout);
    }

    /// Complex or real values in one time sample (`nants * nchan * npol`).
    [[nodiscard]] std::uint64_t values_per_sample() const {
        if (nants == 0U || nchan == 0U || npol == 0U) {
            throw ValidationError("psrio: baseband geometry is empty");
        }
        return nants * nchan * npol;
    }

    [[nodiscard]] SampleType sample_type() const {
        if (nbit < 0) {
            return SampleType::kFloat32;
        }
        switch (nbit) {
        case 1:
            return SampleType::kUInt1;
        case 2:
            return SampleType::kUInt2;
        case 4:
            return SampleType::kUInt4;
        case 8:
            return samples_signed ? SampleType::kInt8 : SampleType::kUInt8;
        case 16:
            return samples_signed ? SampleType::kInt16 : SampleType::kUInt16;
        case 32:
            return SampleType::kUInt32;
        default:
            throw ValidationError("psrio: baseband bit depth is not supported");
        }
    }

    [[nodiscard]] double tobs() const noexcept {
        return static_cast<double>(nsamples) * tsamp;
    }
};

namespace concepts {

/**
 * @brief Packed baseband stream read one block of time samples at a time.
 *
 * The sample cursor counts time. `read_block` copies canonical packed bytes:
 * `[time][antenna][channel][polarization][component]`. `read_samples` unpacks
 * that same layout. A complex stream rejects a real destination.
 */
template <typename R>
concept BasebandReader = requires(R& reader,
                                  const R& reader_const,
                                  std::uint64_t count,
                                  std::int64_t delta,
                                  std::span<std::byte> bytes,
                                  std::span<float> floats,
                                  std::span<std::int8_t> i8s,
                                  std::span<std::int16_t> i16s,
                                  std::span<std::complex<float>> cfloats,
                                  std::span<std::complex<std::int8_t>> ci8s,
                                  std::span<std::complex<std::int16_t>> ci16s) {
    { reader_const.header() } -> std::same_as<const BasebandHeader&>;
    { reader_const.npol() } -> std::same_as<std::uint64_t>;
    { reader_const.nchan() } -> std::same_as<std::uint64_t>;
    { reader_const.nants() } -> std::same_as<std::uint64_t>;
    { reader_const.nbit() } -> std::same_as<int>;
    { reader_const.ndim() } -> std::same_as<int>;
    { reader_const.sample_type() } -> std::same_as<SampleType>;
    { reader_const.samples_signed() } -> std::same_as<bool>;
    { reader_const.msb_first() } -> std::same_as<bool>;
    { reader_const.bytes_per_sample() } -> std::same_as<std::uint64_t>;
    { reader_const.nsamples() } -> std::same_as<std::uint64_t>;
    { reader_const.has_nsamples() } -> std::same_as<bool>;
    { reader_const.tsamp() } -> std::same_as<double>;
    { reader_const.tstart() } -> std::same_as<double>;
    { reader_const.utc_start() } -> std::same_as<std::time_t>;
    { reader_const.fch1() } -> std::same_as<double>;
    { reader_const.foff() } -> std::same_as<double>;
    { reader_const.bandwidth() } -> std::same_as<double>;
    { reader_const.center_frequency() } -> std::same_as<double>;
    { reader_const.source_name() } -> std::convertible_to<std::string_view>;
    { reader_const.telescope() } -> std::convertible_to<std::string_view>;
    { reader_const.tell() } -> std::same_as<std::uint64_t>;
    { reader.seek(count) };
    { reader.rewind() };
    { reader.skip(delta) };
    { reader.read_block(count, bytes) } -> std::same_as<std::uint64_t>;
    { reader.read_bytes(count, bytes) } -> std::same_as<std::uint64_t>;
    {
        reader.read_voltages(count, VoltageRead{}, bytes)
    } -> std::same_as<std::uint64_t>;
    { reader.read_samples(count, floats) } -> std::same_as<std::uint64_t>;
    { reader.read_samples(count, i8s) } -> std::same_as<std::uint64_t>;
    { reader.read_samples(count, i16s) } -> std::same_as<std::uint64_t>;
    { reader.read_samples(count, cfloats) } -> std::same_as<std::uint64_t>;
    { reader.read_samples(count, ci8s) } -> std::same_as<std::uint64_t>;
    { reader.read_samples(count, ci16s) } -> std::same_as<std::uint64_t>;
};

} // namespace concepts

template <typename T> inline constexpr bool kIsBasebandComplex = false;

template <typename T>
inline constexpr bool kIsBasebandComplex<std::complex<T>> = true;

namespace detail {

inline void require_block(std::uint64_t count,
                          std::span<const std::byte> dest,
                          std::uint64_t stride) {
    if (count > 0U &&
        stride > (std::numeric_limits<std::uint64_t>::max() / count)) {
        throw ValidationError("psrio: requested sample block is too large");
    }
    if (dest.size() != count * stride) {
        throw ValidationError(
            "psrio: destination size does not match requested samples");
    }
}

template <typename T>
void require_values(std::uint64_t count,
                    std::span<const T> dest,
                    const BasebandHeader& header) {
    const auto values = header.values_per_sample();
    if (count > 0U &&
        values > (std::numeric_limits<std::uint64_t>::max() / count)) {
        throw ValidationError("psrio: requested sample block is too large");
    }
    if (dest.size() != count * values) {
        throw ValidationError(
            "psrio: destination size does not match requested samples");
    }
    if (header.ndim == 2 && !kIsBasebandComplex<T>) {
        throw ValidationError(
            "psrio: complex baseband samples need a complex destination");
    }
    if (header.ndim == 1 && kIsBasebandComplex<T>) {
        throw ValidationError(
            "psrio: real baseband samples need a real destination");
    }
}

inline void require_dual_pol_bytes(std::uint64_t nsamps,
                                   std::span<const std::byte> dest,
                                   const BasebandHeader& header,
                                   std::uint64_t stride) {
    if (header.nants != 1U || header.npol != 2U || header.ndim != 2 ||
        (header.nbit != 2 && header.nbit != 4 && header.nbit != 8)) {
        throw ValidationError(
            "psrio: voltage export needs one antenna, two polarisations, "
            "and 2-, 4-, or 8-bit complex samples");
    }
    require_block(nsamps, dest, stride);
}

} // namespace detail

/// Canonical packed voltages held in memory.
class MemoryBaseband {
public:
    MemoryBaseband(BasebandHeader header, std::vector<std::byte> samples)
        : m_header(std::move(header)),
          m_owned(std::move(samples)) {
        m_bytes = m_owned;
        init();
    }

    template <typename Span>
        requires(!std::same_as<std::remove_cvref_t<Span>,
                               std::vector<std::byte>>) &&
                    std::convertible_to<Span, std::span<const std::byte>>
    MemoryBaseband(BasebandHeader header, Span samples)
        : m_header(std::move(header)),
          m_bytes(samples) {
        init();
    }

    MemoryBaseband(const MemoryBaseband&)                = delete;
    MemoryBaseband& operator=(const MemoryBaseband&)     = delete;
    MemoryBaseband(MemoryBaseband&&) noexcept            = default;
    MemoryBaseband& operator=(MemoryBaseband&&) noexcept = default;
    ~MemoryBaseband()                                    = default;

    [[nodiscard]] const BasebandHeader& header() const noexcept {
        return m_header;
    }
    [[nodiscard]] std::uint64_t npol() const noexcept { return m_header.npol; }
    [[nodiscard]] std::uint64_t nchan() const noexcept {
        return m_header.nchan;
    }
    [[nodiscard]] std::uint64_t nants() const noexcept {
        return m_header.nants;
    }
    [[nodiscard]] int nbit() const noexcept { return m_header.nbit; }
    [[nodiscard]] int ndim() const noexcept { return m_header.ndim; }
    [[nodiscard]] SampleType sample_type() const {
        return m_header.sample_type();
    }
    [[nodiscard]] bool samples_signed() const noexcept {
        return m_header.samples_signed;
    }
    [[nodiscard]] bool msb_first() const noexcept { return m_header.msb_first; }
    void set_msb_first(bool msb) noexcept { m_header.msb_first = msb; }
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

    std::uint64_t read_block(std::uint64_t count, std::span<std::byte> dest) {
        detail::require_block(count, dest, m_stride);
        const auto take = std::min(count, available());
        if (take == 0U) {
            return 0U;
        }
        std::memcpy(dest.data(),
                    m_bytes.data() +
                        static_cast<std::size_t>(m_sample * m_stride),
                    static_cast<std::size_t>(take * m_stride));
        m_sample += take;
        return take;
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
        const auto packed =
            m_bytes.subspan(static_cast<std::size_t>(m_sample * m_stride),
                            static_cast<std::size_t>(take * m_stride));
        const auto values      = take * m_header.values_per_sample();
        auto geometry          = layout();
        geometry.channel_major = false;
        geometry.pol_major     = false;
        if constexpr (kIsBasebandComplex<T>) {
            detail::unpack_complex(
                packed, dest.first(static_cast<std::size_t>(values)), geometry);
        } else {
            detail::unpack_components(
                packed, dest.first(static_cast<std::size_t>(values)), geometry);
        }
        m_sample += take;
        return take;
    }

    template <typename T>
    [[nodiscard]] std::vector<T> read_samples(std::uint64_t count) {
        const auto values = count * m_header.values_per_sample();
        std::vector<T> out(static_cast<std::size_t>(values));
        const auto actual = read_samples(count, std::span<T>(out));
        out.resize(
            static_cast<std::size_t>(actual * m_header.values_per_sample()));
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
            native.first(static_cast<std::size_t>(got * m_stride)), got, nsamps,
            m_header.nchan, m_header.nbit,
            how.frequency_ascending && m_header.foff < 0.0,
            how.order == VoltageOrder::kFreqMajor, dest);
        return got;
    }

    [[nodiscard]] static std::size_t num_groups() noexcept { return 1; }

    std::uint64_t read_native(std::uint64_t count, std::span<std::byte> dest) {
        return read_block(count, dest);
    }

    std::uint64_t read_groups(std::uint64_t count,
                              std::span<std::span<std::byte>> dest_groups) {
        if (dest_groups.size() != 1U) {
            throw ValidationError(
                "psrio: dest_groups size must match num_groups()");
        }
        return read_voltages(count,
                             VoltageRead{
                                 .order = VoltageOrder::kFreqMajor,
                                 .frequency_ascending = true,
                             },
                             dest_groups[0]);
    }

    [[nodiscard]] std::vector<std::vector<std::byte>>
    read_groups(std::uint64_t count) {
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

    void init() {
        m_header.order        = BasebandOrder::kTimeMajor;
        m_header.has_nsamples = true;
        m_stride              = m_header.bytes_per_sample();
        if (m_stride == 0U ||
            m_bytes.size() % static_cast<std::size_t>(m_stride) != 0U) {
            throw ValidationError("psrio: memory baseband length is not a "
                                  "whole number of samples");
        }
        m_header.nsamples = m_bytes.size() / static_cast<std::size_t>(m_stride);
    }

    BasebandHeader m_header;
    std::uint64_t m_sample{0};
    std::uint64_t m_stride{0};
    mutable std::vector<std::byte> m_scratch_bytes;
    std::vector<std::byte> m_owned;
    std::span<const std::byte> m_bytes;
};

class BasebandSource {
public:
    template <concepts::BasebandReader R>
        requires(!std::same_as<std::remove_cvref_t<R>, BasebandSource>)
    explicit BasebandSource(R reader)
        : m_impl(std::make_unique<Model<R>>(std::move(reader))) {}

    BasebandSource(const BasebandSource&)                = delete;
    BasebandSource& operator=(const BasebandSource&)     = delete;
    BasebandSource(BasebandSource&&) noexcept            = default;
    BasebandSource& operator=(BasebandSource&&) noexcept = default;
    ~BasebandSource()                                    = default;

    [[nodiscard]] const BasebandHeader& header() const {
        return impl().header();
    }
    [[nodiscard]] std::uint64_t npol() const { return impl().npol(); }
    [[nodiscard]] std::uint64_t nchan() const { return impl().nchan(); }
    [[nodiscard]] std::uint64_t nants() const { return impl().nants(); }
    [[nodiscard]] int nbit() const { return impl().nbit(); }
    [[nodiscard]] int ndim() const { return impl().ndim(); }
    [[nodiscard]] SampleType sample_type() const {
        return impl().sample_type();
    }
    [[nodiscard]] bool samples_signed() const {
        return impl().samples_signed();
    }
    [[nodiscard]] bool msb_first() const { return impl().msb_first(); }
    void set_msb_first(bool msb) { impl().set_msb_first(msb); }
    [[nodiscard]] std::uint64_t bytes_per_sample() const {
        return impl().bytes_per_sample();
    }
    [[nodiscard]] std::uint64_t nsamples() const { return impl().nsamples(); }
    [[nodiscard]] bool has_nsamples() const { return impl().has_nsamples(); }
    [[nodiscard]] double tsamp() const { return impl().tsamp(); }
    [[nodiscard]] double tstart() const { return impl().tstart(); }
    [[nodiscard]] std::time_t utc_start() const { return impl().utc_start(); }
    [[nodiscard]] double fch1() const { return impl().fch1(); }
    [[nodiscard]] double foff() const { return impl().foff(); }
    [[nodiscard]] double bandwidth() const { return impl().bandwidth(); }
    [[nodiscard]] double center_frequency() const {
        return impl().center_frequency();
    }
    [[nodiscard]] std::string_view source_name() const {
        return impl().source_name();
    }
    [[nodiscard]] std::string_view telescope() const {
        return impl().telescope();
    }
    [[nodiscard]] std::uint64_t tell() const { return impl().tell(); }

    void seek(std::uint64_t sample) { impl().seek(sample); }
    void rewind() { impl().rewind(); }
    void skip(std::int64_t delta) { impl().skip(delta); }

    std::uint64_t read_block(std::uint64_t count, std::span<std::byte> dest) {
        return impl().read_block(count, dest);
    }
    std::uint64_t read_voltages(std::uint64_t nsamps,
                                VoltageRead how,
                                std::span<std::byte> dest) {
        return impl().read_voltages(nsamps, how, dest);
    }
    std::uint64_t read_bytes(std::uint64_t nbytes, std::span<std::byte> dest) {
        return impl().read_bytes(nbytes, dest);
    }
    std::uint64_t read_samples(std::uint64_t count, std::span<float> dest) {
        return impl().read_samples(count, dest);
    }
    std::uint64_t read_samples(std::uint64_t count,
                               std::span<std::int8_t> dest) {
        return impl().read_samples(count, dest);
    }
    std::uint64_t read_samples(std::uint64_t count,
                               std::span<std::int16_t> dest) {
        return impl().read_samples(count, dest);
    }
    std::uint64_t read_samples(std::uint64_t count,
                               std::span<std::complex<float>> dest) {
        return impl().read_samples(count, dest);
    }
    std::uint64_t read_samples(std::uint64_t count,
                               std::span<std::complex<std::int8_t>> dest) {
        return impl().read_samples(count, dest);
    }
    std::uint64_t read_samples(std::uint64_t count,
                               std::span<std::complex<std::int16_t>> dest) {
        return impl().read_samples(count, dest);
    }

    template <typename T = std::complex<float>>
    [[nodiscard]] std::vector<T> read_samples(std::uint64_t count) {
        const auto values = count * header().values_per_sample();
        std::vector<T> out(static_cast<std::size_t>(values));
        const auto actual = read_samples(count, std::span<T>(out));
        out.resize(
            static_cast<std::size_t>(actual * header().values_per_sample()));
        return out;
    }

    [[nodiscard]] std::vector<std::byte> read_block(std::uint64_t count) {
        std::vector<std::byte> out(
            static_cast<std::size_t>(count * bytes_per_sample()));
        const auto actual = read_block(count, std::span<std::byte>(out));
        out.resize(static_cast<std::size_t>(actual * bytes_per_sample()));
        return out;
    }

    [[nodiscard]] std::size_t num_groups() const { return impl().num_groups(); }

    std::uint64_t read_native(std::uint64_t count, std::span<std::byte> dest) {
        return impl().read_native(count, dest);
    }

    std::uint64_t read_groups(std::uint64_t count,
                              std::span<std::span<std::byte>> dest_groups) {
        return impl().read_groups(count, dest_groups);
    }

    [[nodiscard]] std::vector<std::vector<std::byte>>
    read_groups(std::uint64_t count) {
        return impl().read_groups(count);
    }

    [[nodiscard]] std::uint64_t dropped_packets() const {
        return impl().dropped_packets();
    }

    [[nodiscard]] std::uint64_t dropped_samples() const {
        return impl().dropped_samples();
    }

private:
    struct Concept {
        Concept(const Concept&)            = delete;
        Concept& operator=(const Concept&) = delete;
        Concept(Concept&&)                 = delete;
        Concept& operator=(Concept&&)      = delete;
        virtual ~Concept()                 = default;

        [[nodiscard]] virtual const BasebandHeader& header() const       = 0;
        [[nodiscard]] virtual std::uint64_t npol() const                 = 0;
        [[nodiscard]] virtual std::uint64_t nchan() const                = 0;
        [[nodiscard]] virtual std::uint64_t nants() const                = 0;
        [[nodiscard]] virtual int nbit() const                           = 0;
        [[nodiscard]] virtual int ndim() const                           = 0;
        [[nodiscard]] virtual SampleType sample_type() const             = 0;
        [[nodiscard]] virtual bool samples_signed() const                = 0;
        [[nodiscard]] virtual bool msb_first() const                     = 0;
        virtual void set_msb_first(bool msb)                             = 0;
        [[nodiscard]] virtual std::uint64_t bytes_per_sample() const     = 0;
        [[nodiscard]] virtual std::uint64_t nsamples() const             = 0;
        [[nodiscard]] virtual bool has_nsamples() const                  = 0;
        [[nodiscard]] virtual double tsamp() const                       = 0;
        [[nodiscard]] virtual double tstart() const                      = 0;
        [[nodiscard]] virtual std::time_t utc_start() const              = 0;
        [[nodiscard]] virtual double fch1() const                        = 0;
        [[nodiscard]] virtual double foff() const                        = 0;
        [[nodiscard]] virtual double bandwidth() const                   = 0;
        [[nodiscard]] virtual double center_frequency() const            = 0;
        [[nodiscard]] virtual std::string_view source_name() const       = 0;
        [[nodiscard]] virtual std::string_view telescope() const         = 0;
        [[nodiscard]] virtual std::uint64_t tell() const                 = 0;
        virtual void seek(std::uint64_t sample)                          = 0;
        virtual void rewind()                                            = 0;
        virtual void skip(std::int64_t delta)                            = 0;
        virtual std::uint64_t read_block(std::uint64_t count,
                                         std::span<std::byte> dest)      = 0;
        virtual std::uint64_t read_voltages(std::uint64_t nsamps,
                                            VoltageRead how,
                                            std::span<std::byte> dest)   = 0;
        virtual std::uint64_t read_bytes(std::uint64_t nbytes,
                                         std::span<std::byte> dest)      = 0;
        virtual std::uint64_t read_samples(std::uint64_t count,
                                           std::span<float> dest)        = 0;
        virtual std::uint64_t read_samples(std::uint64_t count,
                                           std::span<std::int8_t> dest)  = 0;
        virtual std::uint64_t read_samples(std::uint64_t count,
                                           std::span<std::int16_t> dest) = 0;
        virtual std::uint64_t
        read_samples(std::uint64_t count,
                     std::span<std::complex<float>> dest) = 0;
        virtual std::uint64_t
        read_samples(std::uint64_t count,
                     std::span<std::complex<std::int8_t>> dest) = 0;
        virtual std::uint64_t
        read_samples(std::uint64_t count,
                     std::span<std::complex<std::int16_t>> dest)     = 0;
        [[nodiscard]] virtual std::size_t num_groups() const         = 0;
        virtual std::uint64_t read_native(std::uint64_t count,
                                          std::span<std::byte> dest) = 0;
        virtual std::uint64_t
        read_groups(std::uint64_t count,
                    std::span<std::span<std::byte>> dest_groups) = 0;
        virtual std::vector<std::vector<std::byte>>
        read_groups(std::uint64_t count)                            = 0;
        [[nodiscard]] virtual std::uint64_t dropped_packets() const = 0;
        [[nodiscard]] virtual std::uint64_t dropped_samples() const = 0;

    protected:
        Concept() = default;
    };

    template <concepts::BasebandReader R> struct Model final : Concept {
        explicit Model(R reader) : m_reader(std::move(reader)) {}

        [[nodiscard]] const BasebandHeader& header() const override {
            return m_reader.header();
        }
        [[nodiscard]] std::uint64_t npol() const override {
            return m_reader.npol();
        }
        [[nodiscard]] std::uint64_t nchan() const override {
            return m_reader.nchan();
        }
        [[nodiscard]] std::uint64_t nants() const override {
            return m_reader.nants();
        }
        [[nodiscard]] int nbit() const override { return m_reader.nbit(); }
        [[nodiscard]] int ndim() const override { return m_reader.ndim(); }
        [[nodiscard]] SampleType sample_type() const override {
            return m_reader.sample_type();
        }
        [[nodiscard]] bool samples_signed() const override {
            return m_reader.samples_signed();
        }
        [[nodiscard]] bool msb_first() const override {
            return m_reader.msb_first();
        }
        void set_msb_first(bool msb) override { m_reader.set_msb_first(msb); }
        [[nodiscard]] std::uint64_t bytes_per_sample() const override {
            return m_reader.bytes_per_sample();
        }
        [[nodiscard]] std::uint64_t nsamples() const override {
            return m_reader.nsamples();
        }
        [[nodiscard]] bool has_nsamples() const override {
            return m_reader.has_nsamples();
        }
        [[nodiscard]] double tsamp() const override { return m_reader.tsamp(); }
        [[nodiscard]] double tstart() const override {
            return m_reader.tstart();
        }
        [[nodiscard]] std::time_t utc_start() const override {
            return m_reader.utc_start();
        }
        [[nodiscard]] double fch1() const override { return m_reader.fch1(); }
        [[nodiscard]] double foff() const override { return m_reader.foff(); }
        [[nodiscard]] double bandwidth() const override {
            return m_reader.bandwidth();
        }
        [[nodiscard]] double center_frequency() const override {
            return m_reader.center_frequency();
        }
        [[nodiscard]] std::string_view source_name() const override {
            return m_reader.source_name();
        }
        [[nodiscard]] std::string_view telescope() const override {
            return m_reader.telescope();
        }
        [[nodiscard]] std::uint64_t tell() const override {
            return m_reader.tell();
        }
        void seek(std::uint64_t sample) override { m_reader.seek(sample); }
        void rewind() override { m_reader.rewind(); }
        void skip(std::int64_t delta) override { m_reader.skip(delta); }
        std::uint64_t read_block(std::uint64_t count,
                                 std::span<std::byte> dest) override {
            return m_reader.read_block(count, dest);
        }
        std::uint64_t read_voltages(std::uint64_t nsamps,
                                    VoltageRead how,
                                    std::span<std::byte> dest) override {
            return m_reader.read_voltages(nsamps, how, dest);
        }
        std::uint64_t read_bytes(std::uint64_t nbytes,
                                 std::span<std::byte> dest) override {
            return m_reader.read_bytes(nbytes, dest);
        }
        std::uint64_t read_samples(std::uint64_t count,
                                   std::span<float> dest) override {
            return m_reader.read_samples(count, dest);
        }
        std::uint64_t read_samples(std::uint64_t count,
                                   std::span<std::int8_t> dest) override {
            return m_reader.read_samples(count, dest);
        }
        std::uint64_t read_samples(std::uint64_t count,
                                   std::span<std::int16_t> dest) override {
            return m_reader.read_samples(count, dest);
        }
        std::uint64_t
        read_samples(std::uint64_t count,
                     std::span<std::complex<float>> dest) override {
            return m_reader.read_samples(count, dest);
        }
        std::uint64_t
        read_samples(std::uint64_t count,
                     std::span<std::complex<std::int8_t>> dest) override {
            return m_reader.read_samples(count, dest);
        }
        std::uint64_t
        read_samples(std::uint64_t count,
                     std::span<std::complex<std::int16_t>> dest) override {
            return m_reader.read_samples(count, dest);
        }
        [[nodiscard]] std::size_t num_groups() const override {
            return m_reader.num_groups();
        }
        std::uint64_t read_native(std::uint64_t count,
                                  std::span<std::byte> dest) override {
            return m_reader.read_native(count, dest);
        }
        std::uint64_t
        read_groups(std::uint64_t count,
                    std::span<std::span<std::byte>> dest_groups) override {
            return m_reader.read_groups(count, dest_groups);
        }
        std::vector<std::vector<std::byte>>
        read_groups(std::uint64_t count) override {
            return m_reader.read_groups(count);
        }
        [[nodiscard]] std::uint64_t dropped_packets() const override {
            return m_reader.dropped_packets();
        }
        [[nodiscard]] std::uint64_t dropped_samples() const override {
            return m_reader.dropped_samples();
        }

        R m_reader;
    };

    [[nodiscard]] const Concept& impl() const {
        if (m_impl == nullptr) {
            throw ValidationError("psrio: baseband source is empty");
        }
        return *m_impl;
    }
    [[nodiscard]] Concept& impl() {
        if (m_impl == nullptr) {
            throw ValidationError("psrio: baseband source is empty");
        }
        return *m_impl;
    }

    std::unique_ptr<Concept> m_impl;
};

[[nodiscard]] inline double
lowest_channel_frequency(const BasebandHeader& header) {
    if (header.nchan == 0U) {
        return header.fch1;
    }
    const auto last =
        header.fch1 + (header.foff * static_cast<double>(header.nchan - 1U));
    return std::min(header.fch1, last);
}

[[nodiscard]] inline double
highest_channel_frequency(const BasebandHeader& header) {
    if (header.nchan == 0U) {
        return header.fch1;
    }
    const auto last =
        header.fch1 + (header.foff * static_cast<double>(header.nchan - 1U));
    return std::max(header.fch1, last);
}

/// Several baseband bands read as one stream with a single channel axis.
///
/// When foff >= 0, bands are ordered ascending in frequency.
/// When foff < 0, bands are ordered descending in frequency so channel 0 is
/// at the highest sky frequency, maintaining monotonic channel ordering.
class FrequencyStitch {
public:
    explicit FrequencyStitch(std::vector<BasebandSource> bands) {
        if (bands.empty()) {
            throw ValidationError("psrio: frequency stitch needs a band");
        }
        if (bands.front().header().foff < 0.0) {
            std::ranges::stable_sort(bands, [](const BasebandSource& left,
                                               const BasebandSource& right) {
                return highest_channel_frequency(left.header()) >
                       highest_channel_frequency(right.header());
            });
        } else {
            std::ranges::stable_sort(bands, [](const BasebandSource& left,
                                               const BasebandSource& right) {
                return lowest_channel_frequency(left.header()) <
                       lowest_channel_frequency(right.header());
            });
        }
        const auto& first      = bands.front().header();
        std::uint64_t channels = 0;
        for (const auto& band : bands) {
            check_compatible(first, band.header());
            const auto antenna_bits =
                band.header().nchan * band.header().npol *
                static_cast<std::uint64_t>(band.header().ndim) *
                static_cast<std::uint64_t>(
                    detail::component_bits(band.header().nbit));
            if (antenna_bits % 8U != 0U) {
                throw ValidationError("psrio: stitched band is not a whole "
                                      "number of bytes per antenna");
            }
            channels += band.header().nchan;
        }
        m_header         = first;
        m_header.nchan   = channels;
        m_header.order   = BasebandOrder::kTimeMajor;
        m_header.overlap = 0;
        m_header.bandwidth =
            std::abs(m_header.foff) * static_cast<double>(channels);
        m_header.center_frequency =
            m_header.fch1 +
            (m_header.foff * (static_cast<double>(channels) - 1.0) / 2.0);
        m_stride = m_header.bytes_per_sample();
        m_bands  = std::move(bands);
    }

    FrequencyStitch(const FrequencyStitch&)                = delete;
    FrequencyStitch& operator=(const FrequencyStitch&)     = delete;
    FrequencyStitch(FrequencyStitch&&) noexcept            = default;
    FrequencyStitch& operator=(FrequencyStitch&&) noexcept = default;
    ~FrequencyStitch()                                     = default;

    [[nodiscard]] const BasebandHeader& header() const noexcept {
        return m_header;
    }
    [[nodiscard]] std::uint64_t npol() const noexcept { return m_header.npol; }
    [[nodiscard]] std::uint64_t nchan() const noexcept {
        return m_header.nchan;
    }
    [[nodiscard]] std::uint64_t nants() const noexcept {
        return m_header.nants;
    }
    [[nodiscard]] int nbit() const noexcept { return m_header.nbit; }
    [[nodiscard]] int ndim() const noexcept { return m_header.ndim; }
    [[nodiscard]] SampleType sample_type() const {
        return m_header.sample_type();
    }
    [[nodiscard]] bool samples_signed() const noexcept {
        return m_header.samples_signed;
    }
    [[nodiscard]] bool msb_first() const noexcept { return m_header.msb_first; }
    void set_msb_first(bool msb) {
        m_header.msb_first = msb;
        for (auto& band : m_bands) {
            band.set_msb_first(msb);
        }
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
        for (auto& band : m_bands) {
            band.seek(sample);
        }
    }

    void rewind() { seek(0); }

    void skip(std::int64_t delta) {
        if (!m_header.has_nsamples) {
            throw ValidationError("psrio: baseband length is not known");
        }
        const auto landing =
            detail::apply_skip(m_sample, delta, m_header.nsamples);
        seek(landing);
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

    std::uint64_t read_block(std::uint64_t count, std::span<std::byte> dest) {
        detail::require_block(count, dest, m_stride);
        const auto take = std::min(count, available());
        if (take == 0U) {
            return 0U;
        }
        const auto origin = m_sample;
        try {
            std::uint64_t antenna_offset = 0;
            const auto out_antenna       = m_stride / m_header.nants;
            for (auto& band : m_bands) {
                const auto band_stride  = band.bytes_per_sample();
                const auto band_antenna = band_stride / band.nants();
                const auto needed =
                    static_cast<std::size_t>(take * band_stride);
                if (m_stitch_scratch.size() < needed) {
                    m_stitch_scratch.resize(needed);
                }
                const std::span<std::byte> scratch(m_stitch_scratch.data(),
                                                   needed);
                const auto got = band.read_block(take, scratch);
                if (got != take) {
                    throw FormatError("psrio: baseband band ended early");
                }
                for (std::uint64_t time = 0; time < take; ++time) {
                    for (std::uint64_t antenna = 0; antenna < m_header.nants;
                         ++antenna) {
                        std::memcpy(dest.data() + static_cast<std::size_t>(
                                                      (time * m_stride) +
                                                      (antenna * out_antenna) +
                                                      antenna_offset),
                                    scratch.data() +
                                        static_cast<std::size_t>(
                                            (time * band_stride) +
                                            (antenna * band_antenna)),
                                    static_cast<std::size_t>(band_antenna));
                    }
                }
                antenna_offset += band_antenna;
            }
        } catch (...) {
            for (auto& band : m_bands) {
                band.seek(origin);
            }
            m_sample = origin;
            throw;
        }
        m_sample = origin + take;
        return take;
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

    std::uint64_t read_native(std::uint64_t count, std::span<std::byte> dest) {
        return read_block(count, dest);
    }

    /// Each band reverses itself when its spacing is negative.
    /// When frequency_ascending is true, bands are arranged in ascending
    /// frequency order.
    std::uint64_t read_voltages(std::uint64_t nsamps,
                                VoltageRead how,
                                std::span<std::byte> dest) {
        detail::require_dual_pol_bytes(nsamps, dest, m_header, m_stride);
        const auto take = std::min(nsamps, available());
        if (take == 0U) {
            return 0U;
        }
        const auto origin = m_sample;
        const auto slot   = static_cast<std::uint64_t>(m_header.nbit) * 4U / 8U;
        const bool freq_major = how.order == VoltageOrder::kFreqMajor;
        try {
            std::uint64_t channel_base = 0;
            const auto process_band    = [&](BasebandSource& band) {
                const auto band_channels = band.nchan();
                const auto needed =
                    static_cast<std::size_t>(take * band.bytes_per_sample());
                if (m_stitch_scratch.size() < needed) {
                    m_stitch_scratch.resize(needed);
                }
                const std::span<std::byte> scratch(m_stitch_scratch.data(),
                                                   needed);
                const auto got = band.read_voltages(take, how, scratch);
                if (got != take) {
                    throw FormatError("psrio: baseband band ended early");
                }
                for (std::uint64_t channel = 0; channel < band_channels;
                     ++channel) {
                    const auto dest_channel = channel_base + channel;
                    if (freq_major) {
                        std::memcpy(
                            dest.data() + static_cast<std::size_t>(
                                              (dest_channel * nsamps) * slot),
                            scratch.data() + static_cast<std::size_t>(
                                                 (channel * take) * slot),
                            static_cast<std::size_t>(take * slot));
                    } else {
                        for (std::uint64_t time = 0; time < take; ++time) {
                            std::memcpy(
                                dest.data() + static_cast<std::size_t>(
                                                  ((time * m_header.nchan) +
                                                   dest_channel) *
                                                  slot),
                                scratch.data() +
                                    static_cast<std::size_t>(
                                        ((time * band_channels) + channel) *
                                        slot),
                                static_cast<std::size_t>(slot));
                        }
                    }
                }
                channel_base += band_channels;
            };

            if (how.frequency_ascending && m_header.foff < 0.0) {
                for (auto it = m_bands.rbegin(); it != m_bands.rend(); ++it) {
                    process_band(*it);
                }
            } else {
                for (auto& band : m_bands) {
                    process_band(band);
                }
            }
        } catch (...) {
            for (auto& band : m_bands) {
                band.seek(origin);
            }
            m_sample = origin;
            throw;
        }
        m_sample = origin + take;
        return take;
    }

    template <typename T>
    std::uint64_t read_samples(std::uint64_t count, std::span<T> dest) {
        detail::require_values(count, std::span<const T>(dest), m_header);
        const auto take = std::min(count, available());
        if (take == 0U) {
            return 0U;
        }
        const auto origin = m_sample;
        try {
            std::uint64_t value_offset = 0;
            const auto out_antenna     = m_header.nchan * m_header.npol;
            auto& band_scratch         = scratch_for<T>();
            for (auto& band : m_bands) {
                const auto band_values  = band.header().values_per_sample();
                const auto band_antenna = band.nchan() * band.npol();
                const auto needed =
                    static_cast<std::size_t>(take * band_values);
                if (band_scratch.size() < needed) {
                    band_scratch.resize(needed);
                }
                const std::span<T> scratch(band_scratch.data(), needed);
                const auto got = band.read_samples(take, scratch);
                if (got != take) {
                    throw FormatError("psrio: baseband band ended early");
                }
                for (std::uint64_t time = 0; time < take; ++time) {
                    for (std::uint64_t antenna = 0; antenna < m_header.nants;
                         ++antenna) {
                        const auto dest_index =
                            (((time * m_header.nants) + antenna) *
                             out_antenna) +
                            value_offset;
                        const auto source_index =
                            ((time * m_header.nants) + antenna) * band_antenna;
                        std::copy_n(
                            scratch.begin() +
                                static_cast<std::ptrdiff_t>(source_index),
                            static_cast<std::ptrdiff_t>(band_antenna),
                            dest.begin() +
                                static_cast<std::ptrdiff_t>(dest_index));
                    }
                }
                value_offset += band_antenna;
            }
        } catch (...) {
            for (auto& band : m_bands) {
                band.seek(origin);
            }
            m_sample = origin;
            throw;
        }
        m_sample = origin + take;
        return take;
    }

    template <typename T>
    [[nodiscard]] std::vector<T> read_samples(std::uint64_t count) {
        const auto values = count * m_header.values_per_sample();
        std::vector<T> out(static_cast<std::size_t>(values));
        const auto actual = read_samples(count, std::span<T>(out));
        out.resize(
            static_cast<std::size_t>(actual * m_header.values_per_sample()));
        return out;
    }

    [[nodiscard]] std::size_t num_groups() const noexcept {
        return m_bands.size();
    }

    std::uint64_t read_groups(std::uint64_t nsamps,
                              std::span<std::span<std::byte>> dest_groups) {
        if (dest_groups.size() != m_bands.size()) {
            throw ValidationError(
                "psrio: dest_groups size must match num_groups()");
        }
        const auto take = std::min(nsamps, available());
        if (take == 0U) {
            return 0U;
        }
        const auto origin = m_sample;
        try {
            const bool ascending_bands = (m_header.foff >= 0.0);
            for (std::size_t i = 0; i < m_bands.size(); ++i) {
                const std::size_t band_idx =
                    ascending_bands ? i : (m_bands.size() - 1U - i);
                auto& band             = m_bands[band_idx];
                const auto band_stride = band.bytes_per_sample();
                const auto needed =
                    static_cast<std::size_t>(take * band_stride);
                if (dest_groups[i].size() < needed) {
                    throw ValidationError(
                        "psrio: dest_group buffer is too small");
                }
                std::span<std::byte> single_group = dest_groups[i];
                const std::span<std::span<std::byte>> single_span(&single_group,
                                                                  1);
                const auto got = band.read_groups(take, single_span);
                if (got != take) {
                    throw FormatError("psrio: baseband band ended early");
                }
            }
        } catch (...) {
            for (auto& band : m_bands) {
                band.seek(origin);
            }
            m_sample = origin;
            throw;
        }
        m_sample = origin + take;
        return take;
    }

    [[nodiscard]] std::vector<std::vector<std::byte>>
    read_groups(std::uint64_t nsamps) {
        const auto n = m_bands.size();
        std::vector<std::vector<std::byte>> buffers(n);
        std::vector<std::span<std::byte>> spans(n);
        const bool ascending_bands = (m_header.foff >= 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t band_idx =
                ascending_bands ? i : (m_bands.size() - 1U - i);
            buffers[i].resize(static_cast<std::size_t>(
                nsamps * m_bands[band_idx].bytes_per_sample()));
            spans[i] = std::span<std::byte>(buffers[i]);
        }
        read_groups(nsamps, spans);
        return buffers;
    }

    [[nodiscard]] std::uint64_t dropped_packets() const {
        std::uint64_t total = 0;
        for (const auto& band : m_bands) {
            total += band.dropped_packets();
        }
        return total;
    }

    [[nodiscard]] std::uint64_t dropped_samples() const {
        std::uint64_t total = 0;
        for (const auto& band : m_bands) {
            total += band.dropped_samples();
        }
        return total;
    }

private:
    static void check_compatible(const BasebandHeader& first,
                                 const BasebandHeader& other) {
        if (first.npol != other.npol || first.nbit != other.nbit ||
            first.ndim != other.ndim || first.nants != other.nants ||
            first.samples_signed != other.samples_signed ||
            first.nsamples != other.nsamples ||
            first.has_nsamples != other.has_nsamples) {
            throw ValidationError(
                "psrio: stitched bands do not share a geometry");
        }
        if (std::abs(first.tsamp - other.tsamp) > 1e-12) {
            throw ValidationError(
                "psrio: stitched bands do not share a sample interval");
        }
        const auto tolerance =
            first.tsamp > 0.0 ? (0.5 * first.tsamp / 86400.0) : 0.0;
        if (std::abs(first.tstart - other.tstart) > tolerance) {
            throw ValidationError(
                "psrio: stitched bands do not share a start time");
        }
        if (std::abs(first.foff - other.foff) > 1e-6) {
            throw ValidationError(
                "psrio: stitched bands do not share a channel spacing");
        }
    }

    template <typename T> std::vector<T>& scratch_for() const {
        if constexpr (std::is_same_v<T, float>) {
            return m_scratch_f;
        } else if constexpr (std::is_same_v<T, std::int8_t>) {
            return m_scratch_i8;
        } else if constexpr (std::is_same_v<T, std::int16_t>) {
            return m_scratch_i16;
        } else if constexpr (std::is_same_v<T, std::complex<float>>) {
            return m_scratch_cf;
        } else if constexpr (std::is_same_v<T, std::complex<std::int8_t>>) {
            return m_scratch_ci8;
        } else if constexpr (std::is_same_v<T, std::complex<std::int16_t>>) {
            return m_scratch_ci16;
        } else {
            static_assert(sizeof(T) == 0, "psrio: unsupported sample type");
        }
    }

    BasebandHeader m_header;
    std::uint64_t m_sample{0};
    std::uint64_t m_stride{0};
    std::vector<BasebandSource> m_bands;
    mutable std::vector<std::byte> m_stitch_scratch;
    mutable std::vector<float> m_scratch_f;
    mutable std::vector<std::int8_t> m_scratch_i8;
    mutable std::vector<std::int16_t> m_scratch_i16;
    mutable std::vector<std::complex<float>> m_scratch_cf;
    mutable std::vector<std::complex<std::int8_t>> m_scratch_ci8;
    mutable std::vector<std::complex<std::int16_t>> m_scratch_ci16;
};

static_assert(concepts::BasebandReader<MemoryBaseband>);
static_assert(concepts::BasebandReader<FrequencyStitch>);

} // namespace psrio
