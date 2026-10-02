#include "psrio/block_source.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/formats/memory.hpp"
#include "psrio/formats/sigproc/reader.hpp"
#include "sigproc_bytes.hpp"

#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <span>
#include <vector>

using psrio::formats::sigproc::FilterbankReader;

namespace {

class TempFile {
public:
    explicit TempFile(std::span<const std::byte> bytes) {
        m_path =
            std::filesystem::temp_directory_path() /
            ("psrio-block-" +
             std::to_string(
                 std::chrono::steady_clock::now().time_since_epoch().count()));
        std::ofstream stream(m_path, std::ios::binary);
        REQUIRE(stream.good());
        stream.write(reinterpret_cast<const char*>(bytes.data()),
                     static_cast<std::streamsize>(bytes.size()));
        REQUIRE(stream.good());
    }

    ~TempFile() { std::filesystem::remove(m_path); }

    TempFile(const TempFile&)            = delete;
    TempFile& operator=(const TempFile&) = delete;
    TempFile(TempFile&&)                 = delete;
    TempFile& operator=(TempFile&&)      = delete;

    [[nodiscard]] const std::filesystem::path& path() const { return m_path; }

private:
    std::filesystem::path m_path;
};

} // namespace

TEST_CASE("filterbank readers model BlockReader", "[block]") {
    static_assert(psrio::concepts::BlockReader<FilterbankReader>);
    static_assert(psrio::concepts::BlockReader<psrio::MemoryBlock>);
    static_assert(psrio::concepts::BlockReader<psrio::BlockSource>);
}

TEST_CASE("read_block short-reads packed samples and skip moves backward",
          "[block]") {
    auto bytes = psrio::test::minimal_filterbank(4, 8, 0.001, 1500.0, -0.5);
    const std::array<std::uint8_t, 12> payload = {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
    };
    for (const auto value : payload) {
        bytes.u8(value);
    }
    const TempFile file(bytes.out);
    FilterbankReader reader(file.path());

    REQUIRE(reader.nchans() == 4U);
    REQUIRE(reader.nbits() == 8);
    REQUIRE(reader.sample_type() == psrio::SampleType::kUInt8);
    REQUIRE(reader.bytes_per_sample() == 4U);
    REQUIRE(reader.nsamples() == 3U);
    REQUIRE(reader.has_nsamples());
    REQUIRE(reader.fch1() == 1500.0);
    REQUIRE(reader.foff() == -0.5);
    REQUIRE(reader.tsamp() == 0.001);
    REQUIRE(reader.spectra_rate() == 1000.0);
    REQUIRE(reader.utc_start() == psrio::astro::mjd_to_time(0.0));

    std::vector<std::byte> block(16);
    REQUIRE(reader.read_block(2, std::span<std::byte>{block}.first(8)) == 2U);
    REQUIRE(reader.tell() == 2U);
    REQUIRE(block[0] == std::byte{1});
    REQUIRE(block[7] == std::byte{8});

    reader.skip(-1);
    REQUIRE(reader.tell() == 1U);
    REQUIRE(reader.read_block(4, block) == 2U);
    REQUIRE(block[0] == std::byte{5});
    REQUIRE(block[7] == std::byte{12});
    REQUIRE(reader.read_block(1, std::span<std::byte>{block}.first(4)) == 0U);

    reader.rewind();
    std::vector<std::byte> too_far(16);
    REQUIRE_THROWS_AS(reader.read_bytes(16, too_far), psrio::ValidationError);
    REQUIRE(reader.tell() == 0U);

    psrio::BlockSource source{FilterbankReader{file.path()}};
    REQUIRE(source.nchans() == 4U);
    REQUIRE(source.read_block(3, block) == 3U);
    source.skip(-3);
    std::vector<float> unpacked(12);
    REQUIRE(source.read_samples(3, unpacked) == 3U);
    REQUIRE(unpacked[0] == 1.0F);
    REQUIRE(unpacked[11] == 12.0F);

    psrio::BlockSource moved = std::move(source);
    REQUIRE_THROWS_AS(source.nchans(), psrio::ValidationError);
    REQUIRE(moved.tell() == 3U);
}

TEST_CASE("BlockSource exposes metadata, fswap, multi-type unpacking, and "
          "allocating helpers",
          "[block]") {
    auto bytes = psrio::test::minimal_filterbank(4, 8, 0.001, 1500.0, -0.5);
    const std::array<std::uint8_t, 8> payload = {
        10, 20, 30, 40, 50, 60, 70, 80,
    };
    for (const auto value : payload) {
        bytes.u8(value);
    }
    const TempFile file(bytes.out);

    psrio::BlockSource source{FilterbankReader{file.path()}};

    // Metadata accessors
    REQUIRE(source.nchans() == 4U);
    REQUIRE(source.nifs() == 1U);
    REQUIRE(source.nbits() == 8);
    REQUIRE(source.sample_type() == psrio::SampleType::kUInt8);
    REQUIRE(source.bytes_per_sample() == 4U);
    REQUIRE(source.nsamples() == 2U);
    REQUIRE(source.has_nsamples());
    REQUIRE(source.tsamp() == 0.001);
    REQUIRE(source.tstart() == 0.0);
    REQUIRE(source.fch1() == 1500.0);
    REQUIRE(source.foff() == -0.5);
    REQUIRE(source.bandwidth() == 2.0);
    REQUIRE(source.center_frequency() == 1500.0 + (-0.5 * 3.0 / 2.0));
    REQUIRE(source.source_name() == "Unknown");
    REQUIRE(source.telescope() == "unknown");
    REQUIRE(source.raj() == 0.0);
    REQUIRE(source.dej() == 0.0);
    REQUIRE(source.beam() == 0);
    REQUIRE(source.spectra_rate() == 1000.0);
    REQUIRE_FALSE(source.fswap_enabled());

    // Allocating read_block without fswap
    const auto block0 = source.read_block(1);
    REQUIRE(block0.size() == 4U);
    REQUIRE(block0[0] == std::byte{10});
    REQUIRE(block0[1] == std::byte{20});
    REQUIRE(block0[2] == std::byte{30});
    REQUIRE(block0[3] == std::byte{40});
    REQUIRE(source.tell() == 1U);

    // Test fswap
    source.set_fswap(true);
    REQUIRE(source.fswap_enabled());
    const auto block1 = source.read_block(1);
    REQUIRE(block1.size() == 4U);
    // Channel-reversed: 50, 60, 70, 80 becomes 80, 70, 60, 50
    REQUIRE(block1[0] == std::byte{80});
    REQUIRE(block1[1] == std::byte{70});
    REQUIRE(block1[2] == std::byte{60});
    REQUIRE(block1[3] == std::byte{50});

    // Rewind and read uint8_t samples
    source.rewind();
    REQUIRE(source.tell() == 0U);
    std::vector<std::uint8_t> u8_out(8);
    REQUIRE(source.read_samples(2, std::span<std::uint8_t>{u8_out}) == 2U);
    REQUIRE(u8_out[0] == 10);
    REQUIRE(u8_out[3] == 40);
    REQUIRE(u8_out[4] == 50);
    REQUIRE(u8_out[7] == 80);

    // 16-bit output on 8-bit source throws ValidationError
    source.rewind();
    std::vector<std::uint16_t> u16_out(8);
    REQUIRE_THROWS_AS(source.read_samples(2, std::span<std::uint16_t>{u16_out}),
                      psrio::ValidationError);

    // Allocating read_samples<float> and read_samples<uint8_t>
    source.rewind();
    const auto float_samples = source.read_samples<float>(2);
    REQUIRE(float_samples.size() == 8U);
    REQUIRE(float_samples[0] == 10.0F);
    REQUIRE(float_samples[7] == 80.0F);

    source.rewind();
    const auto u8_samples = source.read_samples<std::uint8_t>(2);
    REQUIRE(u8_samples.size() == 8U);
    REQUIRE(u8_samples[0] == 10);
    REQUIRE(u8_samples[7] == 80);

    // read_bytes and allocating read_bytes
    source.rewind();
    const auto bytes_out = source.read_bytes(4);
    REQUIRE(bytes_out.size() == 4U);
    REQUIRE(source.tell() == 1U);
}

TEST_CASE("skip rejects a landing index outside the file", "[block]") {
    auto bytes = psrio::test::minimal_filterbank(1, 8);
    bytes.u8(1);
    bytes.u8(2);
    const TempFile file(bytes.out);
    FilterbankReader reader(file.path());
    REQUIRE_THROWS_AS(reader.skip(-1), psrio::ValidationError);
    REQUIRE(reader.tell() == 0U);
    REQUIRE_THROWS_AS(reader.skip(3), psrio::ValidationError);
    reader.skip(2);
    REQUIRE(reader.tell() == 2U);
    std::vector<std::byte> one(1);
    REQUIRE(reader.read_block(1, one) == 0U);
}
