#pragma once

#include <concepts>
#include <cstddef>

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
