#include "psrio/block_source.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/formats/memory.hpp"
#include "psrio/formats/sigproc/reader.hpp"
#include "sigproc_bytes.hpp"

#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <span>
#include <vector>

namespace {

class TempFile {
public:
    explicit TempFile(std::span<const std::byte> bytes) {
        m_path = std::filesystem::temp_directory_path() /
                 ("psrio-mem-" +
                  std::to_string(
                      std::chrono::steady_clock::now().time_since_epoch().count()));
        std::ofstream stream(m_path, std::ios::binary);
        REQUIRE(stream.good());
        stream.write(reinterpret_cast<const char*>(bytes.data()),
                     static_cast<std::streamsize>(bytes.size()));
        REQUIRE(stream.good());
    }

    ~TempFile() { std::filesystem::remove(m_path); }

    [[nodiscard]] const std::filesystem::path& path() const { return m_path; }

private:
    std::filesystem::path m_path;
};

} // namespace

TEST_CASE("MemoryBlock rejects an incomplete stride", "[memory]") {
    psrio::MemoryInfo empty;
    REQUIRE_THROWS_AS(psrio::MemoryBlock(empty, {}), psrio::ValidationError);

    psrio::MemoryInfo odd_bits;
    odd_bits.nchans = 3;
    odd_bits.nbits  = 1;
    REQUIRE_THROWS_AS(psrio::MemoryBlock(odd_bits, std::vector<std::byte>(1)),
                      psrio::ValidationError);

    psrio::MemoryInfo partial;
    partial.nchans = 2;
    partial.nbits  = 8;
    REQUIRE_THROWS_AS(psrio::MemoryBlock(partial, std::vector<std::byte>(3)),
                      psrio::ValidationError);
}

TEST_CASE("MemoryBlock matches FilterbankReader packed bytes and floats",
          "[memory]") {
    auto bytes = psrio::test::minimal_filterbank(4, 8, 6.4e-5, 1400.0, -1.0);
    const std::uint8_t payload[] = {9, 8, 7, 6, 5, 4, 3, 2};
    for (const auto value : payload) {
        bytes.u8(value);
    }
    const auto header_bytes = bytes.out.size() - sizeof(payload);
    std::vector<std::byte> packed(bytes.out.begin() +
                                      static_cast<std::ptrdiff_t>(header_bytes),
                                  bytes.out.end());

    const TempFile file(bytes.out);
    psrio::formats::sigproc::FilterbankReader file_reader(file.path());
    psrio::MemoryInfo info;
    info.nchans = file_reader.nchans();
    info.nifs   = file_reader.nifs();
    info.nbits  = file_reader.nbits();
    info.tsamp  = file_reader.tsamp();
    info.tstart = file_reader.tstart();
    info.fch1   = file_reader.fch1();
    info.foff   = file_reader.foff();
    info.beam   = file_reader.beam();
    psrio::MemoryBlock memory(info, packed);

    std::vector<std::byte> from_file(8);
    std::vector<std::byte> from_memory(8);
    REQUIRE(file_reader.read_block(2, from_file) == 2U);
    REQUIRE(memory.read_block(2, from_memory) == 2U);
    REQUIRE(from_file == from_memory);

    file_reader.skip(-2);
    memory.skip(-2);
    std::vector<float> file_floats(8);
    std::vector<float> memory_floats(8);
    REQUIRE(file_reader.read_samples(2, std::span<float>{file_floats}) == 2U);
    REQUIRE(memory.read_samples(2, std::span<float>{memory_floats}) == 2U);
    REQUIRE(file_floats == memory_floats);
    REQUIRE(memory.fch1() == 1400.0);
    REQUIRE(memory.foff() == -1.0);
    REQUIRE(memory.spectra_rate() == 1.0 / 6.4e-5);

    psrio::BlockSource source{std::move(memory)};
    REQUIRE(source.tell() == 2U);
    REQUIRE(source.read_block(1, std::span<std::byte>{from_memory}.first(4)) ==
            0U);
}

TEST_CASE("MemoryBlock supports non-owning spans, zero-copy views, fswap, and metadata",
          "[memory]") {
    const std::uint8_t raw[] = {1, 2, 3, 4, 5, 6, 7, 8};
    std::span<const std::byte> byte_view{
        reinterpret_cast<const std::byte*>(raw), sizeof(raw)};

    psrio::MemoryInfo info;
    info.nchans      = 4;
    info.nifs        = 1;
    info.nbits       = 8;
    info.tsamp       = 0.002;
    info.fch1        = 1420.0;
    info.foff        = -0.25;
    info.source    = "B1937+21";
    info.telescope = "Arecibo";
    info.raj         = 193700.0;
    info.dej         = 210000.0;

    // Non-owning span constructor
    psrio::MemoryBlock mem(info, byte_view);
    REQUIRE(mem.nchans() == 4U);
    REQUIRE(mem.nsamples() == 2U);
    REQUIRE(mem.source_name() == "B1937+21");
    REQUIRE(mem.telescope() == "Arecibo");
    REQUIRE(mem.raj() == 193700.0);
    REQUIRE(mem.dej() == 210000.0);
    REQUIRE(mem.bandwidth() == 1.0);
    REQUIRE(mem.center_frequency() == 1420.0 + (-0.25 * 3.0 / 2.0));
    REQUIRE_FALSE(mem.fswap_enabled());

    // Zero-copy view_block
    const auto v0 = mem.view_block(1);
    REQUIRE(v0.data() == byte_view.data());
    REQUIRE(v0.size() == 4U);
    REQUIRE(mem.tell() == 1U);

    // Zero-copy view_bytes
    const auto v1 = mem.view_bytes(4);
    REQUIRE(v1.data() == byte_view.data() + 4);
    REQUIRE(v1.size() == 4U);
    REQUIRE(mem.tell() == 2U);

    // Rewind and fswap
    mem.rewind();
    mem.set_fswap(true);
    REQUIRE(mem.fswap_enabled());
    std::vector<std::byte> reversed(4);
    REQUIRE(mem.read_block(1, reversed) == 1U);
    REQUIRE(reversed[0] == std::byte{4});
    REQUIRE(reversed[1] == std::byte{3});
    REQUIRE(reversed[2] == std::byte{2});
    REQUIRE(reversed[3] == std::byte{1});

    // Multi-type unpacking
    mem.rewind();
    std::vector<std::uint8_t> u8_out(8);
    REQUIRE(mem.read_samples(2, std::span<std::uint8_t>{u8_out}) == 2U);
    REQUIRE(u8_out[0] == 1);
    REQUIRE(u8_out[7] == 8);

    // Allocating read_samples<float>
    mem.rewind();
    const auto float_out = mem.read_samples<float>(2);
    REQUIRE(float_out.size() == 8U);
    REQUIRE(float_out[0] == 1.0F);
    REQUIRE(float_out[7] == 8.0F);
}
