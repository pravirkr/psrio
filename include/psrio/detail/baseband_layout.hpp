#pragma once

#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"

#include <complex>
#include <cstdint>
#include <cstring>
#include <limits>
#include <span>
#include <utility>

namespace psrio::detail {

/// Axis description for one baseband stream.
///
/// Canonical order, and the order `copy_to_canonical` writes, is
/// `[time][antenna][channel][polarization][component]`. `channel_major` is the
/// GUPPI `1SFA` on-disk order `[antenna][channel][time][polarization]`.
/// `pol_major` is the DADA order `[time][polarization][channel]` used when
/// both axes are longer than one.
struct BasebandLayout {
    std::uint64_t nants{1};
    std::uint64_t nchan{1};
    std::uint64_t npol{1};
    int nbit{8};
    int ndim{1};
    bool samples_signed{true};
    bool offset_binary{false};
    bool msb_first{false};
    bool channel_major{false};
    bool pol_major{false};
};

[[nodiscard]] inline int component_bits(int nbit) {
    switch (nbit) {
    case 1:
    case 2:
    case 4:
    case 8:
    case 16:
    case 32:
        return nbit;
    case -32:
        return 32;
    default:
        throw ValidationError("psrio: baseband bit depth is not supported");
    }
}

[[nodiscard]] inline std::uint64_t sample_bits(const BasebandLayout& layout) {
    if (layout.nants == 0U || layout.nchan == 0U || layout.npol == 0U ||
        (layout.ndim != 1 && layout.ndim != 2)) {
        throw ValidationError("psrio: baseband geometry is empty");
    }
    const auto bits = static_cast<std::uint64_t>(component_bits(layout.nbit));
    return layout.nants * layout.nchan * layout.npol *
           static_cast<std::uint64_t>(layout.ndim) * bits;
}

[[nodiscard]] inline std::uint64_t sample_stride(const BasebandLayout& layout) {
    const auto bits = sample_bits(layout);
    if (bits % 8U != 0U) {
        throw ValidationError(
            "psrio: baseband sample is not a whole number of bytes");
    }
    return bits / 8U;
}

[[nodiscard]] inline std::uint64_t
values_per_sample(const BasebandLayout& layout) {
    if (layout.nants == 0U || layout.nchan == 0U || layout.npol == 0U) {
        throw ValidationError("psrio: baseband geometry is empty");
    }
    return layout.nants * layout.nchan * layout.npol;
}

inline void apply_offset_binary(std::span<std::byte> bytes) {
    for (std::byte& item : bytes) {
        const auto raw = std::to_integer<std::uint8_t>(item);
        if (raw != 0U) {
            item =
                static_cast<std::byte>(static_cast<std::uint8_t>(raw ^ 0x80U));
        }
    }
}

inline void copy_bits(const std::byte* source,
                      std::uint64_t source_bit,
                      std::byte* dest,
                      std::uint64_t dest_bit,
                      std::uint64_t nbits,
                      bool msb_first) {
    for (std::uint64_t bit = 0; bit < nbits; ++bit) {
        const auto src_index =
            static_cast<std::size_t>((source_bit + bit) / 8U);
        const auto src_shift =
            msb_first ? static_cast<unsigned>(7U - ((source_bit + bit) % 8U))
                      : static_cast<unsigned>((source_bit + bit) % 8U);
        const auto src_byte  = std::to_integer<unsigned>(source[src_index]);
        const auto value     = (src_byte >> src_shift) & 1U;
        const auto dst_index = static_cast<std::size_t>((dest_bit + bit) / 8U);
        const auto dst_shift =
            msb_first ? static_cast<unsigned>(7U - ((dest_bit + bit) % 8U))
                      : static_cast<unsigned>((dest_bit + bit) % 8U);
        const auto mask = static_cast<unsigned char>(1U << dst_shift);
        auto current    = std::to_integer<unsigned char>(dest[dst_index]);
        if (value != 0U) {
            current = static_cast<unsigned char>(current | mask);
        } else {
            current = static_cast<unsigned char>(
                current & static_cast<unsigned char>(~mask));
        }
        dest[dst_index] = static_cast<std::byte>(current);
    }
}

inline void copy_packed(const std::byte* source,
                        std::uint64_t source_bit,
                        std::byte* dest,
                        std::uint64_t dest_bit,
                        std::uint64_t nbits,
                        bool msb_first = false) {
    if (nbits == 0U) {
        return;
    }
    if (source_bit % 8U == 0U && dest_bit % 8U == 0U && nbits % 8U == 0U) {
        std::memcpy(dest + (dest_bit / 8U), source + (source_bit / 8U),
                    static_cast<std::size_t>(nbits / 8U));
        return;
    }
    if (nbits == 4U && source_bit % 4U == 0U && dest_bit % 4U == 0U) {
        const auto src_byte =
            std::to_integer<unsigned>(source[source_bit / 8U]);
        const auto k_src = static_cast<unsigned>((source_bit / 4U) % 2U);
        unsigned src_shift;
        if (msb_first) {
            if (k_src == 0U) {
                src_shift = 4U;
            } else {
                src_shift = 0U;
            }
        } else {
            if (k_src == 0U) {
                src_shift = 0U;
            } else {
                src_shift = 4U;
            }
        }

        const unsigned val = (src_byte >> src_shift) & 0x0FU;

        const auto dst_idx = static_cast<std::size_t>(dest_bit / 8U);
        const auto k_dst   = static_cast<unsigned>((dest_bit / 4U) % 2U);
        unsigned dst_shift;
        if (msb_first) {
            if (k_dst == 0U) {
                dst_shift = 4U;
            } else {
                dst_shift = 0U;
            }
        } else {
            if (k_dst == 0U) {
                dst_shift = 0U;
            } else {
                dst_shift = 4U;
            }
        }

        auto cur      = std::to_integer<unsigned>(dest[dst_idx]);
        cur           = (cur & ~(0x0FU << dst_shift)) | (val << dst_shift);
        dest[dst_idx] = static_cast<std::byte>(cur);
        return;
    }
    copy_bits(source, source_bit, dest, dest_bit, nbits, msb_first);
}

[[nodiscard]] inline std::uint64_t
native_group_bit(const BasebandLayout& layout,
                 std::uint64_t native_ntime,
                 std::uint64_t time,
                 std::uint64_t antenna,
                 std::uint64_t channel,
                 std::uint64_t pol) {
    const auto bits  = static_cast<std::uint64_t>(component_bits(layout.nbit));
    const auto group = static_cast<std::uint64_t>(layout.ndim) * bits;
    if (layout.channel_major) {
        const auto index =
            ((((((antenna * layout.nchan) + channel) * native_ntime) + time) *
              layout.npol) +
             pol);
        return index * group;
    }
    if (layout.pol_major) {
        const auto index =
            ((((time * layout.npol) + pol) * layout.nchan) + channel);
        return index * group;
    }
    const auto index =
        ((((((time * layout.nants) + antenna) * layout.nchan) + channel) *
          layout.npol) +
         pol);
    return index * group;
}

[[nodiscard]] inline std::uint64_t
canonical_group_bit(const BasebandLayout& layout,
                    std::uint64_t time,
                    std::uint64_t antenna,
                    std::uint64_t channel,
                    std::uint64_t pol) {
    const auto bits  = static_cast<std::uint64_t>(component_bits(layout.nbit));
    const auto group = static_cast<std::uint64_t>(layout.ndim) * bits;
    const auto index =
        ((((((time * layout.nants) + antenna) * layout.nchan) + channel) *
          layout.npol) +
         pol);
    return index * group;
}

/// Copy @p count time samples from a native block into canonical packed bytes.
///
/// @p native holds `native_ntime` on-disk time samples. @p time0 is the first
/// sample to copy. @p dest must hold `count` canonical samples.
inline void copy_to_canonical(std::span<const std::byte> native,
                              std::uint64_t native_ntime,
                              std::uint64_t time0,
                              std::uint64_t count,
                              const BasebandLayout& layout,
                              std::span<std::byte> dest) {
    const auto stride = sample_stride(layout);
    if (count == 0U) {
        return;
    }
    if (time0 > native_ntime || count > native_ntime - time0) {
        throw ValidationError("psrio: baseband copy extends past the block");
    }
    if (dest.size() != count * stride) {
        throw ValidationError(
            "psrio: canonical destination does not match the sample count");
    }
    const auto need = (time0 + count) * stride;
    if (!layout.channel_major && native.size() < need) {
        throw ValidationError(
            "psrio: baseband block is shorter than its header");
    }
    if (!layout.channel_major && !layout.pol_major) {
        std::memcpy(dest.data(), native.data() + (time0 * stride),
                    static_cast<std::size_t>(count * stride));
        if (layout.offset_binary) {
            apply_offset_binary(dest);
        }
        return;
    }

    const auto group_bits =
        static_cast<std::uint64_t>(layout.ndim) *
        static_cast<std::uint64_t>(component_bits(layout.nbit));
    std::memset(dest.data(), 0, dest.size());
    for (std::uint64_t local = 0; local < count; ++local) {
        const auto time = time0 + local;
        for (std::uint64_t antenna = 0; antenna < layout.nants; ++antenna) {
            for (std::uint64_t channel = 0; channel < layout.nchan; ++channel) {
                for (std::uint64_t pol = 0; pol < layout.npol; ++pol) {
                    const auto source_bit = native_group_bit(
                        layout, native_ntime, time, antenna, channel, pol);
                    const auto dest_bit = canonical_group_bit(
                        layout, local, antenna, channel, pol);
                    const auto source_byte = source_bit / 8U;
                    if (source_byte >= native.size() && group_bits > 0U) {
                        throw ValidationError(
                            "psrio: baseband block is shorter than its header");
                    }
                    copy_packed(native.data(), source_bit, dest.data(),
                                dest_bit, group_bits, layout.msb_first);
                }
            }
        }
    }
    if (layout.offset_binary) {
        apply_offset_binary(dest);
    }
}

template <typename T> [[nodiscard]] T narrow_component(std::int64_t value) {
    if constexpr (std::same_as<T, float>) {
        return static_cast<float>(value);
    } else {
        if (value < static_cast<std::int64_t>(std::numeric_limits<T>::min()) ||
            value > static_cast<std::int64_t>(std::numeric_limits<T>::max())) {
            throw ValidationError(
                "psrio: baseband sample does not fit in the destination");
        }
        return static_cast<T>(value);
    }
}

template <typename T>
[[nodiscard]] T load_component(const std::byte* data,
                               std::uint64_t bit,
                               int nbit,
                               bool samples_signed,
                               bool offset_binary,
                               bool msb_first = false) {
    if constexpr (std::same_as<T, std::int8_t> ||
                  std::same_as<T, std::int16_t> || std::same_as<T, float>) {
        if (nbit < 0) {
            if constexpr (!std::same_as<T, float>) {
                throw ValidationError(
                    "psrio: float baseband samples need a float destination");
            } else {
                return load_little_endian<float>(data + (bit / 8U));
            }
        }
        if (nbit == 8) {
            auto raw = std::to_integer<std::uint8_t>(data[bit / 8U]);
            if (offset_binary && raw != 0U) {
                raw = static_cast<std::uint8_t>(raw ^ 0x80U);
            }
            if (samples_signed || offset_binary) {
                return narrow_component<T>(static_cast<std::int8_t>(raw));
            }
            return narrow_component<T>(raw);
        }
        if (nbit == 16) {
            if (samples_signed) {
                return narrow_component<T>(
                    load_little_endian<std::int16_t>(data + (bit / 8U)));
            }
            return narrow_component<T>(
                load_little_endian<std::uint16_t>(data + (bit / 8U)));
        }
        if (nbit == 32) {
            if (samples_signed) {
                return narrow_component<T>(
                    load_little_endian<std::int32_t>(data + (bit / 8U)));
            }
            return narrow_component<T>(static_cast<int64_t>(
                load_little_endian<std::uint32_t>(data + (bit / 8U))));
        }

        const auto elem_idx =
            (nbit > 0) ? (bit / static_cast<std::uint64_t>(nbit)) : 0U;
        const auto per_byte = 8U / static_cast<unsigned>(nbit);
        const auto byte_idx = static_cast<std::size_t>(elem_idx / per_byte);
        const auto k        = static_cast<unsigned>(elem_idx % per_byte);
        const auto shift = msb_first
                               ? (8U - (static_cast<unsigned>(nbit) * (k + 1U)))
                               : (static_cast<unsigned>(nbit) * k);
        const auto unbit = static_cast<unsigned>(nbit);
        const auto mask  = (1U << unbit) - 1U;
        const auto byte_val = std::to_integer<unsigned>(data[byte_idx]);
        const auto code     = (byte_val >> shift) & mask;

        std::int64_t val = 0;
        const auto half = static_cast<int>(1U << (unbit - 1U));
        if (samples_signed) {
            if (std::cmp_greater_equal(code, half)) {
                val = static_cast<std::int64_t>(static_cast<int>(code) -
                                                (2 * half));
            } else {
                val = static_cast<std::int64_t>(code);
            }
        } else if (offset_binary) {
            val = static_cast<std::int64_t>(static_cast<int>(code) - half);
        } else {
            val = static_cast<std::int64_t>(code);
        }
        return narrow_component<T>(val);
    } else {
        static_assert(sizeof(T) == 0, "unsupported baseband destination");
    }
}

template <typename T>
void unpack_components(std::span<const std::byte> packed,
                       std::span<T> dest,
                       const BasebandLayout& layout) {
    const auto bits = static_cast<std::uint64_t>(component_bits(layout.nbit));
    if (dest.empty()) {
        return;
    }
    const auto need_bits = dest.size() * bits;
    if (need_bits % 8U != 0U || need_bits / 8U != packed.size()) {
        throw ValidationError(
            "psrio: baseband packed bytes do not match the destination");
    }
    if constexpr (std::same_as<T, std::int8_t>) {
        if (layout.nbit == 8 && layout.samples_signed &&
            !layout.offset_binary) {
            std::memcpy(dest.data(), packed.data(), packed.size());
            return;
        }
    }
    for (std::size_t index = 0; index < dest.size(); ++index) {
        dest[index] = load_component<T>(packed.data(), index * bits,
                                        layout.nbit, layout.samples_signed,
                                        layout.offset_binary, layout.msb_first);
    }
}

template <typename T>
void unpack_complex(std::span<const std::byte> packed,
                    std::span<std::complex<T>> dest,
                    const BasebandLayout& layout) {
    if (layout.ndim != 2) {
        throw ValidationError(
            "psrio: complex destination needs complex baseband samples");
    }
    static_assert(sizeof(std::complex<T>) == 2U * sizeof(T));
    if constexpr (std::same_as<T, std::int8_t>) {
        if (layout.nbit == 8 && layout.samples_signed &&
            !layout.offset_binary && packed.size() == dest.size() * 2U) {
            std::memcpy(dest.data(), packed.data(), packed.size());
            return;
        }
    }
    const auto bits = static_cast<std::uint64_t>(component_bits(layout.nbit));
    const auto need_bits = dest.size() * 2U * bits;
    if (need_bits % 8U != 0U || need_bits / 8U != packed.size()) {
        throw ValidationError(
            "psrio: baseband packed bytes do not match the destination");
    }
    for (std::size_t index = 0; index < dest.size(); ++index) {
        const auto real = load_component<T>(
            packed.data(), (index * 2U) * bits, layout.nbit,
            layout.samples_signed, layout.offset_binary, layout.msb_first);
        const auto imag = load_component<T>(
            packed.data(), ((index * 2U) + 1U) * bits, layout.nbit,
            layout.samples_signed, layout.offset_binary, layout.msb_first);
        dest[index] = std::complex<T>(real, imag);
    }
}

/// Rearrange a time-major dual-polarisation block.
///
/// @p native is `[time][channel][pol][real, imag]` for @p nsamps times.
/// @p dest is sized for @p dest_nsamps times. Only the first @p nsamps times
/// are written. `freq_major` stores `[channel][time][pol][real, imag]`.
/// `reverse_channels` maps native channel 0 onto the last destination channel.
inline void arrange_dual_pol(std::span<const std::byte> native,
                             std::uint64_t nsamps,
                             std::uint64_t dest_nsamps,
                             std::uint64_t nchan,
                             int nbits,
                             bool reverse_channels,
                             bool freq_major,
                             std::span<std::byte> dest) {
    const auto slot = static_cast<std::uint64_t>(nbits) * 4U / 8U;
    if (slot == 0U || nsamps > dest_nsamps) {
        throw ValidationError("psrio: voltage block geometry is invalid");
    }
    const auto native_bytes = nsamps * nchan * slot;
    if (native.size() < native_bytes) {
        throw ValidationError("psrio: voltage block is shorter than requested");
    }
    if (!reverse_channels && !freq_major) {
        std::memcpy(dest.data(), native.data(),
                    static_cast<std::size_t>(native_bytes));
        return;
    }
    for (std::uint64_t time = 0; time < nsamps; ++time) {
        for (std::uint64_t channel = 0; channel < nchan; ++channel) {
            const auto dest_channel =
                reverse_channels ? (nchan - 1U - channel) : channel;
            const auto source = ((time * nchan) + channel) * slot;
            const auto target =
                freq_major ? (((dest_channel * dest_nsamps) + time) * slot)
                           : (((time * nchan) + dest_channel) * slot);
            std::memcpy(dest.data() + static_cast<std::size_t>(target),
                        native.data() + static_cast<std::size_t>(source),
                        static_cast<std::size_t>(slot));
        }
    }
}

} // namespace psrio::detail
