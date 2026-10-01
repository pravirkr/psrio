#pragma once

#include "psrio/detail/endian.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <span>
#include <string_view>
#include <vector>

namespace psrio::test {

/// Little-endian SIGPROC byte builder used only by the test suite.
class SigprocBytes {
public:
    void bytes(std::span<const std::byte> extra) {
        out.insert(out.end(), extra.begin(), extra.end());
    }

    void u8(std::uint8_t value) { out.push_back(std::byte{value}); }

    void u32(std::uint32_t value) { append_integer(value); }

    void i32(std::int32_t value) { append_integer(value); }

    void f64(double value) { append_integer(value); }

    void string(std::string_view value) {
        u32(static_cast<std::uint32_t>(value.size()));
        for (const char character : value) {
            out.push_back(std::byte{static_cast<unsigned char>(character)});
        }
    }

    void key_i32(std::string_view key, std::int32_t value) {
        string(key);
        i32(value);
    }

    void key_u32(std::string_view key, std::uint32_t value) {
        string(key);
        u32(value);
    }

    void key_f64(std::string_view key, double value) {
        string(key);
        f64(value);
    }

    void key_i8(std::string_view key, std::int8_t value) {
        string(key);
        out.push_back(std::byte{static_cast<unsigned char>(value)});
    }

    void key_str(std::string_view key, std::string_view value) {
        string(key);
        string(value);
    }

    std::vector<std::byte> out;

private:
    template <typename T> void append_integer(T value) {
        std::array<unsigned char, sizeof(T)> raw{};
        std::memcpy(raw.data(), &value, raw.size());
        if (!::psrio::detail::host_is_little_endian()) {
            std::reverse(raw.begin(), raw.end());
        }
        for (const unsigned char byte : raw) {
            out.push_back(std::byte{byte});
        }
    }
};

[[nodiscard]] inline SigprocBytes minimal_filterbank(std::int32_t nchans = 1,
                                                     std::int32_t nbits  = 8,
                                                     double tsamp        = 1.0,
                                                     double fch1 = 1400.0,
                                                     double foff = -1.0) {
    SigprocBytes bytes;
    bytes.string("HEADER_START");
    bytes.key_i32("nchans", nchans);
    bytes.key_i32("nbits", nbits);
    bytes.key_f64("tsamp", tsamp);
    bytes.key_f64("fch1", fch1);
    bytes.key_f64("foff", foff);
    bytes.string("HEADER_END");
    return bytes;
}

} // namespace psrio::test
