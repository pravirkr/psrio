#pragma once

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <span>

namespace psrio::detail {

template <typename T>
concept ByteLike =
    std::same_as<T, std::byte> || std::same_as<T, const std::byte>;

template <typename T>
concept ByteSpan = requires(T span) {
    { span.data() } -> std::convertible_to<const std::byte*>;
    { span.size() } -> std::convertible_to<std::size_t>;
};

} // namespace psrio::detail

namespace psrio::concepts {

// Multi-channel intensity streams model psrio::concepts::BlockReader in
// psrio/block_source.hpp. That concept is separate from TimeSeriesReader.

/**
 * @brief Concept modeling a single-channel time-series streaming reader.
 *
 * Modeled by psrio::TimeSeriesReader, allowing downstream algorithms (FFTs,
 * periodic folding, single-pulse boxcars) to operate generically on any
 * dedispersed stream.
 */
template <typename R>
concept TimeSeriesReader =
    requires(R& reader, std::span<float> dest, std::uint64_t count) {
        { reader.header() };
        { reader.tell() } -> std::same_as<std::uint64_t>;
        { reader.seek(count) };
        { reader.rewind() };
        { reader.read_samples(count, dest) } -> std::same_as<std::uint64_t>;
        { reader.read(count, dest) } -> std::same_as<std::uint64_t>;
    };

} // namespace psrio::concepts
