#pragma once

#include "psrio/detail/endian.hpp"

#include <bit>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace psrio::test {

/// Minimal FITS writer for search-mode reader tests.
class FitsFileBytes {
public:
    void begin_hdu() { m_hdu = m_out.size(); }

    void key(std::string_view name, std::string_view value) {
        if (name.size() > 8 || value.size() > 70) {
            throw std::invalid_argument("psrio test: FITS card does not fit");
        }
        std::string card(80, ' ');
        card.replace(0, name.size(), name);
        card[8] = '=';
        card[9] = ' ';
        card.replace(10, value.size(), value);
        append(card);
    }

    void end_header() {
        std::string card(80, ' ');
        card.replace(0, 3, "END");
        append(card);
        const auto used = m_out.size() - m_hdu;
        const auto pad  = (2880U - (used % 2880U)) % 2880U;
        m_out.insert(m_out.end(), pad, std::byte{0});
    }

    void raw(std::span<const std::byte> bytes) {
        m_out.insert(m_out.end(), bytes.begin(), bytes.end());
    }

    void pad_data() {
        const auto pad = (2880U - (m_out.size() % 2880U)) % 2880U;
        m_out.insert(m_out.end(), pad, std::byte{0});
    }

    [[nodiscard]] const std::vector<std::byte>& bytes() const { return m_out; }

    void write(const std::filesystem::path& path) const {
        std::ofstream out(path, std::ios::binary);
        out.write(reinterpret_cast<const char*>(m_out.data()),
                  static_cast<std::streamsize>(m_out.size()));
        if (!out) {
            throw std::runtime_error("psrio test: failed to write FITS file");
        }
    }

private:
    void append(std::string_view card) {
        for (const char character : card) {
            m_out.push_back(std::byte{static_cast<unsigned char>(character)});
        }
    }

    std::vector<std::byte> m_out;
    std::size_t m_hdu{0};
};

inline void append_be(std::vector<std::byte>& out, std::uint8_t value) {
    out.push_back(std::byte{value});
}

inline void append_be(std::vector<std::byte>& out, std::int32_t value) {
    if (detail::host_is_little_endian()) {
        value = detail::swap_endian(value);
    }
    const auto* raw = reinterpret_cast<const std::byte*>(&value);
    out.insert(out.end(), raw, raw + sizeof(value));
}

inline void append_be(std::vector<std::byte>& out, std::int16_t value) {
    if (detail::host_is_little_endian()) {
        value = detail::swap_endian(value);
    }
    const auto* raw = reinterpret_cast<const std::byte*>(&value);
    out.insert(out.end(), raw, raw + sizeof(value));
}

inline void append_be(std::vector<std::byte>& out, float value) {
    auto bits = std::bit_cast<std::uint32_t>(value);
    if (detail::host_is_little_endian()) {
        bits = detail::swap_endian(bits);
    }
    const auto* raw = reinterpret_cast<const std::byte*>(&bits);
    out.insert(out.end(), raw, raw + sizeof(bits));
}

inline void append_be(std::vector<std::byte>& out, double value) {
    auto bits = std::bit_cast<std::uint64_t>(value);
    if (detail::host_is_little_endian()) {
        bits = detail::swap_endian(bits);
    }
    const auto* raw = reinterpret_cast<const std::byte*>(&bits);
    out.insert(out.end(), raw, raw + sizeof(bits));
}

} // namespace psrio::test
