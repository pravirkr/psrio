#include "psrio/common/types.hpp"
#include "psrio/detail/endian.hpp"
#include "psrio/detail/exceptions.hpp"
#include "psrio/formats/sigproc.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <limits>
#include <span>
#include <sstream>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {

using psrio::BitOrder;
using psrio::formats::sigproc::ExtraKey;
using psrio::formats::sigproc::FilterbankHeader;
using psrio::formats::sigproc::FilterbankReader;
using psrio::formats::sigproc::FilterbankWriter;

class TempPath {
public:
    TempPath() {
        static std::atomic<int> counter{0};
        m_path = std::filesystem::temp_directory_path() /
                 ("psrio_writer_" + std::to_string(::getpid()) + "_" +
                  std::to_string(counter.fetch_add(1)) + ".fil");
    }
    ~TempPath() { std::filesystem::remove(m_path); }

    TempPath(const TempPath&)            = delete;
    TempPath& operator=(const TempPath&) = delete;
    TempPath(TempPath&&)                 = delete;
    TempPath& operator=(TempPath&&)      = delete;

    [[nodiscard]] const std::filesystem::path& path() const noexcept {
        return m_path;
    }

private:
    std::filesystem::path m_path;
};

std::vector<std::byte> read_file(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    std::vector<char> raw((std::istreambuf_iterator<char>(in)),
                          std::istreambuf_iterator<char>());
    std::vector<std::byte> bytes(raw.size());
    for (std::size_t i = 0; i < raw.size(); ++i) {
        bytes[i] = static_cast<std::byte>(raw[i]);
    }
    return bytes;
}

FilterbankHeader
basic_header(std::int32_t nchans, std::int32_t nbits, std::int32_t nifs = 1) {
    FilterbankHeader header;
    header.nchans = nchans;
    header.nbits  = nbits;
    header.nifs   = nifs;
    header.tsamp  = 64e-6;
    header.tstart = 60000.5;
    header.fch1   = 1500.0;
    header.foff   = -1.0;
    return header;
}

ExtraKey extra_key(std::string name, std::vector<std::byte> raw) {
    ExtraKey key;
    key.name = std::move(name);
    key.raw  = std::move(raw);
    return key;
}

std::vector<std::byte> u32_bytes(std::uint32_t value) {
    std::vector<std::byte> raw(sizeof(value));
    const auto stored = psrio::detail::to_little_endian(value);
    std::memcpy(raw.data(), &stored, sizeof(value));
    return raw;
}

/// Deterministic levels in [0, 2^nbits) for nbits <= 16.
float level_at(std::size_t index, int nbits) {
    if (nbits == 32) {
        return (static_cast<float>(index) * 0.37F) - 11.5F;
    }
    const auto modulus = std::size_t{1} << static_cast<unsigned>(nbits);
    return static_cast<float>(((index * 7U) + (index / 5U)) % modulus);
}

} // namespace

TEST_CASE("FilterbankHeader::to_bytes round-trips every field",
          "[sigproc][writer]") {
    FilterbankHeader header  = basic_header(4, 8, 2);
    header.frequency_table   = {1500.0, 1499.5, 1498.25, 1497.0};
    header.data_type         = 1;
    header.telescope_id      = 4;
    header.machine_id        = 10;
    header.source_name       = "J0437-4715";
    header.rawdatafile       = "raw.fil";
    header.src_raj           = 43715.3;
    header.src_dej           = -471509.2;
    header.az_start          = 12.5;
    header.za_start          = 33.25;
    header.barycentric       = 1;
    header.pulsarcentric     = 0;
    header.ibeam             = 3;
    header.nbeams            = 13;
    header.refdm             = 2.64;
    header.period            = 0.00575;
    header.nbins             = 128;
    header.npuls             = 7;
    header.signed_data       = 1;
    header.declared_nsamples = 99;
    header.extra.push_back(extra_key("custom_int", u32_bytes(42)));
    std::vector<std::byte> text = u32_bytes(3);
    for (const char c : std::string("abc")) {
        text.push_back(static_cast<std::byte>(c));
    }
    header.extra.push_back(extra_key("custom_str", text));

    const auto bytes  = header.to_bytes();
    const auto parsed = FilterbankHeader::parse(bytes);

    CHECK(parsed.header_bytes == bytes.size());
    CHECK(parsed.nchans == 4);
    CHECK(parsed.nbits == 8);
    CHECK(parsed.nifs == 2);
    CHECK(parsed.tsamp == header.tsamp);
    CHECK(parsed.tstart == header.tstart);
    CHECK(parsed.fch1 == header.fch1);
    CHECK(parsed.foff == header.foff);
    CHECK(parsed.frequency_table == header.frequency_table);
    CHECK(parsed.data_type == header.data_type);
    CHECK(parsed.telescope_id == header.telescope_id);
    CHECK(parsed.machine_id == header.machine_id);
    CHECK(parsed.source_name == header.source_name);
    CHECK(parsed.rawdatafile == header.rawdatafile);
    CHECK(parsed.src_raj == header.src_raj);
    CHECK(parsed.src_dej == header.src_dej);
    CHECK(parsed.az_start == header.az_start);
    CHECK(parsed.za_start == header.za_start);
    CHECK(parsed.barycentric == header.barycentric);
    CHECK(parsed.pulsarcentric == header.pulsarcentric);
    CHECK(parsed.ibeam == header.ibeam);
    CHECK(parsed.nbeams == header.nbeams);
    CHECK(parsed.refdm == header.refdm);
    CHECK(parsed.period == header.period);
    CHECK(parsed.nbins == header.nbins);
    CHECK(parsed.npuls == header.npuls);
    CHECK(parsed.signed_data == header.signed_data);
    CHECK(parsed.declared_nsamples == header.declared_nsamples);
    REQUIRE(parsed.extra.size() == 2U);
    CHECK(parsed.extra[0].name == "custom_int");
    CHECK(parsed.extra[0].raw == header.extra[0].raw);
    CHECK(parsed.extra[1].name == "custom_str");
    CHECK(parsed.extra[1].raw == header.extra[1].raw);

    CHECK(parsed.to_bytes() == bytes);

    std::ostringstream stream;
    header.write_to(stream);
    const auto streamed = stream.str();
    REQUIRE(streamed.size() == bytes.size());
    CHECK(std::memcmp(streamed.data(), bytes.data(), bytes.size()) == 0);
}

TEST_CASE("FilterbankHeader::to_bytes defaults and omissions",
          "[sigproc][writer]") {
    FilterbankHeader header  = basic_header(2, 8);
    header.source_name       = "";
    header.declared_nsamples = 0;

    const auto parsed = FilterbankHeader::parse(header.to_bytes());
    CHECK(parsed.data_type == 1);
    CHECK_FALSE(parsed.source_name.has_value());
    CHECK_FALSE(parsed.declared_nsamples.has_value());
    CHECK_FALSE(parsed.telescope_id.has_value());
    CHECK(parsed.extra.empty());
}

TEST_CASE("FilterbankHeader::to_bytes rejects unwritable headers",
          "[sigproc][writer]") {
    FilterbankHeader header = basic_header(4, 8);

    SECTION("bad nbits") {
        header.nbits = 3;
        CHECK_THROWS_AS(header.to_bytes(), psrio::ValidationError);
    }
    SECTION("no channels") {
        header.nchans = 0;
        CHECK_THROWS_AS(header.to_bytes(), psrio::ValidationError);
    }
    SECTION("frequency table length") {
        header.frequency_table = {1.0, 2.0};
        CHECK_THROWS_AS(header.to_bytes(), psrio::ValidationError);
    }
    SECTION("over-long string") {
        header.source_name = std::string(5000, 'x');
        CHECK_THROWS_AS(header.to_bytes(), psrio::ValidationError);
    }
    SECTION("more extra keys than schema keys to separate them") {
        for (int i = 0; i < 32; ++i) {
            header.extra.push_back(
                extra_key("k" + std::to_string(i), u32_bytes(1)));
        }
        CHECK_THROWS_AS(header.to_bytes(), psrio::ValidationError);
    }
}

TEST_CASE("telescope and machine ids from names", "[sigproc][writer]") {
    using psrio::formats::sigproc::machine_id_from_name;
    using psrio::formats::sigproc::telescope_id_from_name;
    using psrio::formats::sigproc::telescope_name;
    const auto parkes = telescope_id_from_name("parkes");
    REQUIRE(parkes.has_value());
    CHECK(telescope_name(*parkes) == "Parkes");
    CHECK_FALSE(telescope_id_from_name("Nowhere").has_value());
    CHECK_FALSE(machine_id_from_name("").has_value());
}

TEST_CASE("FilterbankWriter round-trips every bit depth and order",
          "[sigproc][writer]") {
    const int nbits       = GENERATE(1, 2, 4, 8, 16, 32);
    const BitOrder order  = GENERATE(BitOrder::kLsbFirst, BitOrder::kMsbFirst);
    constexpr int kNchans = 8;
    constexpr int kNifs   = 2;
    constexpr std::uint64_t kNsamps = 37;
    const auto values_per_sample    = static_cast<std::size_t>(kNchans) * kNifs;

    std::vector<float> levels(kNsamps * values_per_sample);
    for (std::size_t i = 0; i < levels.size(); ++i) {
        levels[i] = level_at(i, nbits);
    }

    const TempPath first;
    {
        FilterbankWriter writer(first.path(),
                                basic_header(kNchans, nbits, kNifs), order);
        const std::span<const float> all(levels);
        writer.write_samples(10, all.first(10 * values_per_sample));
        writer.write_samples(kNsamps - 10, all.subspan(10 * values_per_sample));
        CHECK(writer.nsamples_written() == kNsamps);
        CHECK(writer.header().data_bytes ==
              kNsamps * writer.bytes_per_sample());
        writer.close();
    }

    FilterbankReader reader(first.path(), order);
    REQUIRE(reader.nsamples() == kNsamps);
    CHECK(reader.nifs() == kNifs);
    CHECK(reader.read_samples<float>(kNsamps) == levels);

    SECTION("write_block copies the payload unchanged") {
        reader.rewind();
        const auto payload =
            reader.view_bytes(kNsamps * reader.bytes_per_sample());
        const TempPath second;
        {
            FilterbankWriter writer(second.path(), reader.header(), order);
            writer.write_block(5, payload);
            writer.write_bytes(payload.subspan(5 * reader.bytes_per_sample()));
        }
        CHECK(read_file(second.path()) == read_file(first.path()));
    }

    if (nbits <= 8) {
        SECTION("uint8_t input") {
            std::vector<std::uint8_t> bytes(levels.size());
            for (std::size_t i = 0; i < levels.size(); ++i) {
                bytes[i] = static_cast<std::uint8_t>(levels[i]);
            }
            const TempPath second;
            {
                FilterbankWriter writer(second.path(),
                                        basic_header(kNchans, nbits, kNifs));
                writer.write_samples(
                    kNsamps, std::span<const std::uint8_t>(bytes), order);
            }
            CHECK(read_file(second.path()) == read_file(first.path()));
        }
    }
    if (nbits == 16) {
        SECTION("uint16_t input") {
            std::vector<std::uint16_t> words(levels.size());
            for (std::size_t i = 0; i < levels.size(); ++i) {
                words[i] = static_cast<std::uint16_t>(levels[i]);
            }
            const TempPath second;
            {
                FilterbankWriter writer(second.path(),
                                        basic_header(kNchans, nbits, kNifs));
                writer.write(kNsamps, std::span<const std::uint16_t>(words));
            }
            CHECK(read_file(second.path()) == read_file(first.path()));
        }
    }
}

TEST_CASE("FilterbankWriter never writes nsamples", "[sigproc][writer]") {
    FilterbankHeader header  = basic_header(4, 8);
    header.declared_nsamples = 1000;
    const TempPath file;
    {
        FilterbankWriter writer(file.path(), header);
        CHECK_FALSE(writer.header().declared_nsamples.has_value());
        const std::vector<std::uint8_t> data(12, 7);
        writer.write_samples(3, std::span<const std::uint8_t>(data));
    }
    const FilterbankReader reader(file.path());
    CHECK_FALSE(reader.header().declared_nsamples.has_value());
    CHECK(reader.nsamples() == 3U);
}

TEST_CASE("FilterbankWriter quantises float input", "[sigproc][writer]") {
    const float nan = std::numeric_limits<float>::quiet_NaN();

    SECTION("unsigned 8-bit rounds and clamps") {
        const std::vector<float> input{-3.0F, 0.4F, 0.6F, 254.6F, 300.0F, nan};
        const TempPath file;
        {
            FilterbankWriter writer(file.path(), basic_header(6, 8));
            writer.write_samples(1, std::span<const float>(input));
        }
        FilterbankReader reader(file.path());
        CHECK(reader.read_samples<float>(1) ==
              std::vector<float>{0.0F, 0.0F, 1.0F, 255.0F, 255.0F, 0.0F});
    }
    SECTION("signed 8-bit") {
        FilterbankHeader header = basic_header(4, 8);
        header.signed_data      = 1;
        const std::vector<float> input{-200.0F, -1.4F, 3.6F, 200.0F};
        const TempPath file;
        {
            FilterbankWriter writer(file.path(), header);
            writer.write_samples(1, std::span<const float>(input));
        }
        FilterbankReader reader(file.path());
        CHECK(reader.read_samples<float>(1) ==
              std::vector<float>{-128.0F, -1.0F, 4.0F, 127.0F});
    }
    SECTION("2-bit clamps to three") {
        const std::vector<float> input{-1.0F, 1.2F, 2.5F, 9.0F};
        const TempPath file;
        {
            FilterbankWriter writer(file.path(), basic_header(4, 2));
            writer.write_samples(1, std::span<const float>(input));
        }
        FilterbankReader reader(file.path());
        CHECK(reader.read_samples<float>(1) ==
              std::vector<float>{0.0F, 1.0F, 3.0F, 3.0F});
    }
}

TEST_CASE("FilterbankWriter rejects bad writes", "[sigproc][writer]") {
    const TempPath file;

    SECTION("invalid header") {
        FilterbankHeader header = basic_header(4, 8);
        header.tsamp.reset();
        CHECK_THROWS_AS(FilterbankWriter(file.path(), header),
                        psrio::ValidationError);
    }
    SECTION("signed data needs 8 bits") {
        FilterbankHeader header = basic_header(4, 16);
        header.signed_data      = 1;
        CHECK_THROWS_AS(FilterbankWriter(file.path(), header),
                        psrio::ValidationError);
    }
    SECTION("unwritable path") {
        CHECK_THROWS_AS(FilterbankWriter(file.path() / "missing" / "x.fil",
                                         basic_header(4, 8)),
                        psrio::IoError);
    }

    FilterbankWriter writer(file.path(), basic_header(4, 16));
    const std::vector<std::byte> bytes(10);
    const std::vector<float> floats(7);
    const std::vector<std::uint8_t> small(8);

    SECTION("partial sample of bytes") {
        CHECK_THROWS_AS(writer.write_bytes(bytes), psrio::ValidationError);
    }
    SECTION("block shorter than count") {
        CHECK_THROWS_AS(writer.write_block(2, bytes), psrio::ValidationError);
    }
    SECTION("too few values") {
        CHECK_THROWS_AS(writer.write_samples(2, std::span<const float>(floats)),
                        psrio::ValidationError);
    }
    SECTION("uint8_t needs 8 bits or fewer") {
        CHECK_THROWS_AS(
            writer.write_samples(1, std::span<const std::uint8_t>(small)),
            psrio::ValidationError);
    }
    SECTION("closed writer") {
        writer.close();
        writer.close();
        CHECK_FALSE(writer.is_open());
        CHECK_THROWS_AS(writer.write_bytes(std::span<const std::byte>{}),
                        psrio::IoError);
    }
    CHECK(writer.nsamples_written() == 0U);
}
