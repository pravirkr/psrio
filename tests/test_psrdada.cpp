#include "psrio/block_source.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/formats/psrdada.hpp"

#include <catch2/catch_test_macros.hpp>

#include <ascii_header.h>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <functional>
#include <ipcbuf.h>
#include <ipcio.h>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <unistd.h>

using psrio::BlockSource;
using psrio::IoError;
using psrio::SampleType;
using psrio::ValidationError;
using psrio::formats::psrdada::RingReader;

namespace {

void set_header(char* header, const char* key, const char* fmt, auto value) {
    if (ascii_header_set(header, key, fmt, value) < 0) {
        throw std::runtime_error(std::string("ascii_header_set failed for ") +
                                 key);
    }
}

class RingPair {
public:
    RingPair(const RingPair&)            = delete;
    RingPair& operator=(const RingPair&) = delete;

    RingPair(RingPair&& other) noexcept
        : m_key(other.m_key),
          m_header(other.m_header),
          m_data(other.m_data),
          m_header_bufsz(other.m_header_bufsz),
          m_live(other.m_live) {
        other.m_live = false;
    }

    RingPair& operator=(RingPair&&) = delete;

    ~RingPair() {
        if (!m_live) {
            return;
        }
        ipcio_destroy(&m_data);
        ipcbuf_destroy(&m_header);
    }

    static RingPair create(std::uint64_t data_bytes) {
        const auto page = static_cast<std::uint64_t>(sysconf(_SC_PAGESIZE));
        const std::uint64_t header_bufsz = page;
        const std::uint64_t pages =
            data_bytes == 0U ? 4U : ((data_bytes + page - 1U) / page) + 2U;
        std::mt19937 rng{std::random_device{}()};
        std::uniform_int_distribution<int> dist(0x1000, 0x7ffe);
        for (int attempt = 0; attempt < 32; ++attempt) {
            const auto key  = static_cast<key_t>(dist(rng) & ~1);
            ipcbuf_t header = IPCBUF_INIT;
            ipcio_t data    = IPCIO_INIT;
            if (ipcbuf_create(&header, key + 1, 4, header_bufsz, 1) < 0) {
                continue;
            }
            if (ipcio_create(&data, key, pages, page, 1) < 0) {
                ipcbuf_destroy(&header);
                continue;
            }
            return RingPair(key, header, data, header_bufsz);
        }
        throw std::runtime_error("could not create a PSRDADA ring");
    }

    [[nodiscard]] key_t key() const { return m_key; }

    void write_header(const std::function<void(char*)>& fill) {
        ipcbuf_t writer = IPCBUF_INIT;
        if (ipcbuf_connect(&writer, m_key + 1) < 0) {
            throw std::runtime_error("header connect failed");
        }
        struct Disconnect {
            ipcbuf_t* buf;
            ~Disconnect() { ipcbuf_disconnect(buf); }
        } disconnect{&writer};

        if (ipcbuf_lock_write(&writer) < 0) {
            throw std::runtime_error("header lock_write failed");
        }
        char* buf = ipcbuf_get_next_write(&writer);
        if (buf == nullptr) {
            throw std::runtime_error("header get_next_write failed");
        }
        std::memset(buf, 0, static_cast<std::size_t>(m_header_bufsz));
        fill(buf);
        if (ipcbuf_mark_filled(&writer, m_header_bufsz) < 0) {
            throw std::runtime_error("header mark_filled failed");
        }
        if (ipcbuf_unlock_write(&writer) < 0) {
            throw std::runtime_error("header unlock_write failed");
        }
    }

    void write_data(std::span<const std::uint8_t> bytes) {
        ipcio_t writer = IPCIO_INIT;
        if (ipcio_connect(&writer, m_key) < 0) {
            throw std::runtime_error("data connect failed");
        }
        struct Disconnect {
            ipcio_t* buf;
            bool open = false;
            ~Disconnect() {
                if (open) {
                    ipcio_close(buf);
                }
                ipcio_disconnect(buf);
            }
        } disconnect{&writer, false};

        // 'W' raises start-of-data on the first buffer. 'w' leaves the ring
        // inactive, so a reader blocks forever on a partial buffer.
        if (ipcio_open(&writer, 'W') < 0) {
            throw std::runtime_error("data open failed");
        }
        disconnect.open = true;
        if (!bytes.empty()) {
            const auto wrote = ipcio_write(
                &writer,
                const_cast<char*>(reinterpret_cast<const char*>(bytes.data())),
                bytes.size());
            if (wrote != static_cast<ssize_t>(bytes.size())) {
                throw std::runtime_error("data write failed");
            }
        }
        if (ipcio_stop(&writer) < 0) {
            throw std::runtime_error("data stop failed");
        }
    }

private:
    RingPair(key_t key,
             ipcbuf_t header,
             ipcio_t data,
             std::uint64_t header_bufsz)
        : m_key(key),
          m_header(header),
          m_data(data),
          m_header_bufsz(header_bufsz) {}

    key_t m_key;
    ipcbuf_t m_header{};
    ipcio_t m_data{};
    std::uint64_t m_header_bufsz;
    bool m_live = true;
};

void fill_base(char* header, int nchan, int nbit) {
    set_header(header, "NCHAN", "%d", nchan);
    set_header(header, "NBIT", "%d", nbit);
    set_header(header, "NPOL", "%d", 1);
    set_header(header, "NDIM", "%d", 1);
    set_header(header, "ORDER", "%s", "TF");
    set_header(header, "CFREQ", "%f", 1400.0);
    set_header(header, "BANDWIDTH", "%f", -200.0);
    set_header(header, "TSAMP", "%f", 64.0);
    set_header(header, "UTC_START", "%s", "2020-01-02-03:04:05");
    set_header(header, "BEAM", "%d", 7);
}

std::time_t expected_utc() {
    std::tm broken{};
    if (strptime("2020-01-02-03:04:05", "%Y-%m-%d-%H:%M:%S", &broken) ==
        nullptr) {
        throw std::runtime_error("strptime failed");
    }
    return timegm(&broken);
}

std::vector<std::byte> as_bytes(const std::vector<std::uint8_t>& bytes) {
    std::vector<std::byte> out(bytes.size());
    for (std::size_t index = 0; index < bytes.size(); ++index) {
        out[index] = std::byte{bytes[index]};
    }
    return out;
}

} // namespace

static_assert(psrio::concepts::BlockReader<RingReader>);
static_assert(std::is_constructible_v<BlockSource, RingReader>);

TEST_CASE("PSRDADA header sets frequency, sampling and beam", "[psrdada]") {
    auto ring = RingPair::create(0);
    ring.write_header([](char* header) { fill_base(header, 4, 8); });
    ring.write_data({});

    RingReader source(static_cast<std::uint32_t>(ring.key()));
    CHECK(source.nchans() == 4U);
    CHECK(source.nifs() == 1U);
    CHECK(source.nbits() == 8);
    CHECK(source.sample_type() == SampleType::kUInt8);
    CHECK(source.bytes_per_sample() == 4U);
    CHECK(source.fch1() == 1475.0);
    CHECK(source.foff() == -50.0);
    CHECK(source.spectra_rate() == 1.0e6 / 64.0);
    CHECK(source.tsamp() == 64.0 / 1.0e6);
    CHECK(source.utc_start() == expected_utc());
    CHECK(source.beam() == 7);
    CHECK(source.tstart() == 0.0);
    CHECK_FALSE(source.has_nsamples());
    CHECK(source.nsamples() == 0U);
    CHECK(source.tell() == 0U);

    std::vector<std::byte> dest(4);
    CHECK(source.read_block(1, dest) == 0U);
    CHECK(source.has_nsamples());
    CHECK(source.nsamples() == 0U);
}

TEST_CASE("PSRDADA read returns the byte pattern then end of data",
          "[psrdada]") {
    constexpr int kNChan            = 4;
    constexpr std::uint64_t kNSamps = 16;
    auto ring                       = RingPair::create(kNChan * kNSamps);
    ring.write_header([](char* header) { fill_base(header, kNChan, 8); });
    std::vector<std::uint8_t> bytes(kNChan * kNSamps);
    for (std::size_t index = 0; index < bytes.size(); ++index) {
        bytes[index] = static_cast<std::uint8_t>(index);
    }
    ring.write_data(bytes);

    RingReader source(static_cast<std::uint32_t>(ring.key()));
    const auto expected = as_bytes(bytes);
    std::vector<std::byte> got(expected.size());
    CHECK(source.read_block(kNSamps, got) == kNSamps);
    CHECK(got == expected);
    CHECK(source.tell() == kNSamps);

    std::vector<std::byte> extra(expected.size(), std::byte{7});
    CHECK(source.read_block(kNSamps, extra) == 0U);
    CHECK(source.has_nsamples());
    CHECK(source.nsamples() == kNSamps);

    source.rewind();
    CHECK(source.tell() == 0U);
    std::vector<float> floats(kNSamps * static_cast<std::uint64_t>(kNChan));
    CHECK(source.read_samples(kNSamps, floats) == kNSamps);
    for (std::size_t index = 0; index < floats.size(); ++index) {
        CHECK(floats[index] == static_cast<float>(index));
    }

    BlockSource front(std::move(source));
    CHECK(front.nchans() == 4U);
    CHECK(front.sample_type() == SampleType::kUInt8);
    std::vector<std::byte> tail(4);
    CHECK(front.read_block(1, tail) == 0U);
}

TEST_CASE("PSRDADA 4-bit samples use a packed stride", "[psrdada]") {
    constexpr int kNChan            = 8;
    constexpr std::uint64_t kStride = 4;
    constexpr std::uint64_t kNSamps = 10;
    auto ring                       = RingPair::create(kStride * kNSamps);
    ring.write_header([](char* header) { fill_base(header, kNChan, 4); });
    std::vector<std::uint8_t> bytes(kStride * kNSamps);
    for (std::size_t index = 0; index < bytes.size(); ++index) {
        bytes[index] = static_cast<std::uint8_t>(0xA0 + index);
    }
    ring.write_data(bytes);

    RingReader source(static_cast<std::uint32_t>(ring.key()));
    CHECK(source.sample_type() == SampleType::kUInt4);
    CHECK(source.bytes_per_sample() == kStride);
    const auto expected = as_bytes(bytes);
    std::vector<std::byte> got(expected.size());
    CHECK(source.read_block(kNSamps, got) == kNSamps);
    CHECK(got == expected);
}

TEST_CASE("PSRDADA skip moves back inside the current page", "[psrdada]") {
    constexpr int kNChan            = 4;
    constexpr std::uint64_t kNSamps = 16;
    auto ring                       = RingPair::create(kNChan * kNSamps);
    ring.write_header([](char* header) { fill_base(header, kNChan, 8); });
    std::vector<std::uint8_t> bytes(kNChan * kNSamps);
    for (std::size_t index = 0; index < bytes.size(); ++index) {
        bytes[index] = static_cast<std::uint8_t>(index);
    }
    ring.write_data(bytes);

    RingReader source(static_cast<std::uint32_t>(ring.key()));
    CHECK_THROWS_AS(source.skip(-1), ValidationError);

    std::vector<std::byte> first(16);
    REQUIRE(source.read_block(4, first) == 4U);
    CHECK(source.tell() == 4U);
    source.skip(-2);
    CHECK(source.tell() == 2U);

    std::vector<std::byte> again(16);
    REQUIRE(source.read_block(4, again) == 4U);
    const auto expected = as_bytes(bytes);
    CHECK(std::vector<std::byte>(expected.begin() + 8, expected.begin() + 24) ==
          again);
    CHECK(source.tell() == 6U);
}

TEST_CASE("PSRDADA rejects non-intensity headers and a missing ring",
          "[psrdada]") {
    auto expect_invalid = [](const std::function<void(char*)>& fill) {
        auto ring = RingPair::create(0);
        ring.write_header(fill);
        REQUIRE_THROWS_AS(RingReader(static_cast<std::uint32_t>(ring.key())),
                          ValidationError);
    };

    expect_invalid([](char* header) {
        fill_base(header, 4, 8);
        set_header(header, "NPOL", "%d", 2);
    });
    expect_invalid([](char* header) {
        fill_base(header, 4, 8);
        set_header(header, "ORDER", "%s", "FT");
    });
    expect_invalid([](char* header) {
        fill_base(header, 4, 8);
        set_header(header, "NDIM", "%d", 2);
    });
    expect_invalid([](char* header) {
        set_header(header, "NBIT", "%d", 8);
        set_header(header, "NPOL", "%d", 1);
        set_header(header, "TSAMP", "%f", 64.0);
        set_header(header, "CFREQ", "%f", 1400.0);
        set_header(header, "BANDWIDTH", "%f", -200.0);
        set_header(header, "UTC_START", "%s", "2020-01-02-03:04:05");
    });

    std::mt19937 rng{std::random_device{}()};
    std::uniform_int_distribution<int> dist(0x1000, 0x7ffe);
    bool threw = false;
    for (int attempt = 0; attempt < 32; ++attempt) {
        const auto key = static_cast<key_t>(dist(rng) & ~1);
        ipcio_t probe  = IPCIO_INIT;
        if (ipcio_connect(&probe, key) == 0) {
            ipcio_disconnect(&probe);
            continue;
        }
        REQUIRE_THROWS_AS(RingReader(static_cast<std::uint32_t>(key)), IoError);
        threw = true;
        break;
    }
    REQUIRE(threw);
}
