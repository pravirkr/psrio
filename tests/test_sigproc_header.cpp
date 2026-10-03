#include "psrio/detail/exceptions.hpp"
#include "psrio/formats/sigproc/header.hpp"
#include "psrio/psrio.hpp"
#include "sigproc_bytes.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>

using Catch::Matchers::ContainsSubstring;
using psrio::formats::sigproc::FilterbankHeader;
using psrio::test::minimal_filterbank;
using psrio::test::SigprocBytes;

namespace {

void require_format(std::span<const std::byte> bytes, const char* fragment) {
    REQUIRE_THROWS_WITH(FilterbankHeader::parse(bytes),
                        ContainsSubstring(fragment));
}

} // namespace

TEST_CASE("library version is exposed", "[sigproc]") {
    REQUIRE(psrio::version() == "0.3.0");
}

TEST_CASE("minimal filterbank header parses and validates",
          "[sigproc][header]") {
    const auto bytes  = minimal_filterbank();
    const auto header = FilterbankHeader::parse(bytes.out);
    header.validate();
    REQUIRE(header.nchans == 1);
    REQUIRE(header.nifs == 1);
    REQUIRE(header.nbits == 8);
    REQUIRE(header.tsamp == 1.0);
    REQUIRE(header.fch1 == 1400.0);
    REQUIRE(header.foff == -1.0);
    REQUIRE(header.bytes_per_sample() == 1U);
    REQUIRE(header.nsamples() == 0U);
    REQUIRE(header.warnings().empty());
    REQUIRE_FALSE(header.telescope_id.has_value());
    REQUIRE(header.telescope_name() == "unknown");
    REQUIRE_FALSE(header.samples_are_signed());
}

TEST_CASE("optional sigproc keys round-trip", "[sigproc][header]") {
    SigprocBytes bytes;
    bytes.string("HEADER_START");
    bytes.key_i32("nchans", 4);
    bytes.key_i32("nifs", 2);
    bytes.key_i32("nbits", 8);
    bytes.key_f64("tsamp", 6.4e-5);
    bytes.key_f64("tstart", 59000.5);
    bytes.key_f64("fch1", 1500.0);
    bytes.key_f64("foff", -0.5);
    bytes.string("source_name");
    bytes.string("HEADER_END");
    bytes.string("rawdatafile");
    bytes.string("raw.dat");
    bytes.key_i32("telescope_id", 6);
    bytes.key_i32("machine_id", 10);
    bytes.key_i32("data_type", 1);
    bytes.key_i32("barycentric", 1);
    bytes.key_i32("pulsarcentric", 0);
    bytes.key_f64("az_start", 12.5);
    bytes.key_f64("za_start", 40.0);
    bytes.key_f64("src_raj", 123456.78);
    bytes.key_f64("src_dej", -451234.5);
    bytes.key_i32("ibeam", 3);
    bytes.key_i32("nbeams", 64);
    bytes.key_f64("refdm", 26.5);
    bytes.key_f64("period", 0.125);
    bytes.key_i32("nbins", 128);
    bytes.key_i32("npuls", 7);
    bytes.key_i8("signed", -1);
    bytes.key_u32("nsamples", 0);
    bytes.string("HEADER_END");
    bytes.u8(1);
    bytes.u8(2);
    bytes.u8(3);
    bytes.u8(4);
    bytes.u8(5);
    bytes.u8(6);
    bytes.u8(7);
    bytes.u8(8);
    bytes.u8(9);

    const auto header = FilterbankHeader::parse(bytes.out, "roundtrip.fil");
    header.validate();
    REQUIRE(header.source_name == "HEADER_END");
    REQUIRE(header.rawdatafile == "raw.dat");
    REQUIRE(header.telescope_id == 6);
    REQUIRE(header.telescope_name() == "GBT");
    REQUIRE(header.machine_id == 10);
    REQUIRE(header.machine_name() == "BPSR");
    REQUIRE(header.data_type == 1);
    REQUIRE(header.data_type_name() == "filterbank");
    REQUIRE(header.barycentric == 1);
    REQUIRE(header.pulsarcentric == 0);
    REQUIRE(header.az_start == 12.5);
    REQUIRE(header.za_start == 40.0);
    REQUIRE(header.src_raj == 123456.78);
    REQUIRE(header.src_dej == -451234.5);
    REQUIRE(header.ibeam == 3);
    REQUIRE(header.nbeams == 64);
    REQUIRE(header.refdm == 26.5);
    REQUIRE(header.period == 0.125);
    REQUIRE(header.nbins == 128);
    REQUIRE(header.npuls == 7);
    REQUIRE(header.tstart == 59000.5);
    REQUIRE(header.samples_are_signed());
    REQUIRE(header.nifs == 2);
    REQUIRE(header.bytes_per_sample() == 8U);
    REQUIRE(header.declared_nsamples == 0U);
    REQUIRE(header.samples_in_file() == 1U);
    REQUIRE(header.nsamples() == 1U);
    REQUIRE(header.trailing_bytes() == 1U);
    REQUIRE(header.warnings().empty());

    REQUIRE(header.bandwidth() == 2.0);
    REQUIRE(header.ftop() == 1500.25);
    REQUIRE(header.fbottom() == 1498.25);
    REQUIRE(header.fcenter() == 1499.25);
    std::array<double, 4> frequencies{};
    header.channel_frequencies(frequencies);
    REQUIRE(frequencies ==
            std::array<double, 4>{1500.0, 1499.5, 1499.0, 1498.5});
}

TEST_CASE("duplicate header keys keep the last value", "[sigproc][header]") {
    SigprocBytes bytes;
    bytes.string("HEADER_START");
    bytes.key_i32("nchans", 4);
    bytes.key_i32("nchans", 8);
    bytes.key_i32("nbits", 8);
    bytes.key_f64("tsamp", 1.0);
    bytes.key_f64("fch1", 1400.0);
    bytes.key_f64("foff", -1.0);
    bytes.string("HEADER_END");
    const auto header = FilterbankHeader::parse(bytes.out);
    REQUIRE(header.nchans == 8);
}

TEST_CASE("frequency table supplies channel centres", "[sigproc][header]") {
    SigprocBytes bytes;
    bytes.string("HEADER_START");
    bytes.key_i32("nchans", 2);
    bytes.key_i32("nbits", 8);
    bytes.key_f64("tsamp", 1.0);
    bytes.string("FREQUENCY_START");
    bytes.key_f64("fchannel", 1400.0);
    bytes.key_f64("fchannel", 1410.5);
    bytes.string("FREQUENCY_END");
    bytes.string("HEADER_END");
    const auto header = FilterbankHeader::parse(bytes.out);
    header.validate();
    REQUIRE_FALSE(header.fch1.has_value());
    std::array<double, 2> frequencies{};
    header.channel_frequencies(frequencies);
    REQUIRE(frequencies.at(0) == 1400.0);
    REQUIRE(frequencies.at(1) == 1410.5);
    REQUIRE_THROWS_AS(header.bandwidth(), psrio::ValidationError);
}

TEST_CASE("frequency table wins over fch1 and foff", "[sigproc][header]") {
    SigprocBytes bytes;
    bytes.string("HEADER_START");
    bytes.key_i32("nchans", 1);
    bytes.key_i32("nbits", 8);
    bytes.key_f64("tsamp", 1.0);
    bytes.key_f64("fch1", 100.0);
    bytes.key_f64("foff", -1.0);
    bytes.key_f64("fchannel", 1420.25);
    bytes.string("HEADER_END");
    const auto header = FilterbankHeader::parse(bytes.out);
    std::array<double, 1> frequencies{};
    header.channel_frequencies(frequencies);
    REQUIRE(frequencies.at(0) == 1420.25);
    REQUIRE(header.fch1 == 100.0);
}

TEST_CASE("frequency table length must match nchans", "[sigproc][header]") {
    SigprocBytes bytes;
    bytes.string("HEADER_START");
    bytes.key_i32("nchans", 2);
    bytes.key_i32("nbits", 8);
    bytes.key_f64("tsamp", 1.0);
    bytes.key_f64("fchannel", 1400.0);
    bytes.string("HEADER_END");
    const auto header = FilterbankHeader::parse(bytes.out);
    REQUIRE_THROWS_AS(header.validate(), psrio::ValidationError);
}

TEST_CASE("omitted nifs defaults to one", "[sigproc][header]") {
    const auto header = FilterbankHeader::parse(minimal_filterbank(4, 8).out);
    REQUIRE(header.nifs == 1);
    REQUIRE(header.bytes_per_sample() == 4U);
}

TEST_CASE("nifs widens the sample stride", "[sigproc][header]") {
    SigprocBytes file;
    file.string("HEADER_START");
    file.key_i32("nchans", 4);
    file.key_i32("nifs", 2);
    file.key_i32("nbits", 8);
    file.key_f64("tsamp", 1.0);
    file.key_f64("fch1", 1400.0);
    file.key_f64("foff", -1.0);
    file.string("HEADER_END");
    const auto header = FilterbankHeader::parse(file.out);
    REQUIRE(header.bytes_per_sample() == 8U);
}

TEST_CASE("declared nsamples uses the smaller count and warns",
          "[sigproc][header]") {
    const auto with_samples = [](std::uint32_t declared, std::size_t payload) {
        SigprocBytes file;
        file.string("HEADER_START");
        file.key_i32("nchans", 1);
        file.key_i32("nbits", 8);
        file.key_f64("tsamp", 1.0);
        file.key_f64("fch1", 1400.0);
        file.key_f64("foff", -1.0);
        if (declared != 0xFFFFFFFFU) {
            file.key_u32("nsamples", declared);
        }
        file.string("HEADER_END");
        file.out.insert(file.out.end(), payload, std::byte{7});
        return FilterbankHeader::parse(file.out);
    };

    const auto equal = with_samples(4, 4);
    equal.validate();
    REQUIRE(equal.nsamples() == 4U);
    REQUIRE(equal.warnings().empty());

    const auto smaller = with_samples(3, 10);
    smaller.validate();
    REQUIRE(smaller.nsamples() == 3U);
    REQUIRE(smaller.samples_in_file() == 10U);
    REQUIRE(smaller.trailing_bytes() == 7U);
    REQUIRE(smaller.warnings().size() == 1U);
    REQUIRE_THAT(smaller.warnings().front(), ContainsSubstring("disagrees"));

    const auto larger = with_samples(20, 10);
    larger.validate();
    REQUIRE(larger.nsamples() == 10U);
    REQUIRE(larger.trailing_bytes() == 0U);
    REQUIRE(larger.warnings().size() == 1U);

    const auto unset = with_samples(0, 10);
    REQUIRE(unset.declared_nsamples == 0U);
    REQUIRE(unset.nsamples() == 10U);
    REQUIRE(unset.warnings().empty());

    const auto missing = with_samples(0xFFFFFFFFU, 10);
    REQUIRE_FALSE(missing.declared_nsamples.has_value());
    REQUIRE(missing.nsamples() == 10U);
    REQUIRE(missing.warnings().empty());
}

TEST_CASE("malformed sigproc headers throw format errors",
          "[sigproc][header]") {
    require_format({}, "empty");

    SigprocBytes foo;
    foo.string("FOO");
    require_format(foo.out, "HEADER_START");

    SigprocBytes truncated_key;
    truncated_key.string("HEADER_START");
    truncated_key.u8(1);
    truncated_key.u8(2);
    require_format(truncated_key.out, "HEADER_END");

    SigprocBytes truncated_value;
    truncated_value.string("HEADER_START");
    truncated_value.string("nchans");
    truncated_value.u8(1);
    truncated_value.u8(2);
    require_format(truncated_value.out, "truncated value");

    SigprocBytes truncated_string;
    truncated_string.string("HEADER_START");
    truncated_string.string("source_name");
    truncated_string.u32(5);
    truncated_string.u8('a');
    require_format(truncated_string.out, "truncated string");

    SigprocBytes zero_length;
    zero_length.string("HEADER_START");
    zero_length.u32(0);
    require_format(zero_length.out, "implausible string length");

    SigprocBytes too_long;
    too_long.u32(4097);
    require_format(too_long.out, "implausible string length");

    SigprocBytes swapped;
    swapped.u8(0x00);
    swapped.u8(0x00);
    swapped.u8(0x00);
    swapped.u8(0x0C);
    require_format(swapped.out, "big-endian");
}

TEST_CASE("unknown keys are skipped when the next token is recognisable",
          "[sigproc][header]") {
    SigprocBytes bytes;
    bytes.string("HEADER_START");
    bytes.key_i32("nchans", 1);
    bytes.string("custom");
    bytes.i32(42);
    bytes.key_i32("nbits", 8);
    bytes.string("scale");
    bytes.f64(0.0);
    bytes.key_f64("tsamp", 1.0);
    bytes.string("flag");
    bytes.u8(0x7F);
    bytes.key_f64("fch1", 1400.0);
    bytes.key_f64("foff", -1.0);
    bytes.string("HEADER_END");

    const auto header = FilterbankHeader::parse(bytes.out);
    header.validate();
    REQUIRE(header.extra.size() == 3U);
    REQUIRE(header.extra.at(0).name == "custom");
    REQUIRE(header.extra.at(0).raw.size() == 4U);
    REQUIRE(header.extra.at(1).name == "scale");
    REQUIRE(header.extra.at(1).raw.size() == 8U);
    REQUIRE(header.extra.at(2).name == "flag");
    REQUIRE(header.extra.at(2).raw.size() == 1U);
    REQUIRE(std::to_integer<unsigned char>(header.extra.at(2).raw.at(0)) ==
            0x7F);
    REQUIRE(header.nbits == 8);
}

TEST_CASE("an unknown key that cannot be skipped is an error",
          "[sigproc][header]") {
    SigprocBytes bytes;
    bytes.string("HEADER_START");
    bytes.key_i32("nchans", 1);
    bytes.key_i32("nbits", 8);
    bytes.string("weird");
    for (int index = 0; index < 16; ++index) {
        bytes.u8(0xFF);
    }
    require_format(bytes.out, "unknown header key");
}

TEST_CASE("invalid filterbank geometry fails validation", "[sigproc][header]") {
    const auto expect = [](SigprocBytes bytes) {
        const auto header = FilterbankHeader::parse(bytes.out);
        REQUIRE_THROWS_AS(header.validate(), psrio::ValidationError);
    };

    SigprocBytes zero_channels;
    zero_channels.string("HEADER_START");
    zero_channels.key_i32("nchans", 0);
    zero_channels.key_i32("nbits", 8);
    zero_channels.key_f64("tsamp", 1.0);
    zero_channels.key_f64("fch1", 1400.0);
    zero_channels.key_f64("foff", -1.0);
    zero_channels.string("HEADER_END");
    expect(zero_channels);

    for (const int nbits : {0, 3, 12, 64}) {
        expect(minimal_filterbank(1, nbits));
    }

    expect(minimal_filterbank(3, 1));

    SigprocBytes no_time;
    no_time.string("HEADER_START");
    no_time.key_i32("nchans", 1);
    no_time.key_i32("nbits", 8);
    no_time.key_f64("fch1", 1400.0);
    no_time.key_f64("foff", -1.0);
    no_time.string("HEADER_END");
    expect(no_time);

    SigprocBytes no_frequency;
    no_frequency.string("HEADER_START");
    no_frequency.key_i32("nchans", 1);
    no_frequency.key_i32("nbits", 8);
    no_frequency.key_f64("tsamp", 1.0);
    no_frequency.string("HEADER_END");
    expect(no_frequency);

    SigprocBytes series;
    series.string("HEADER_START");
    series.key_i32("data_type", 2);
    series.key_i32("nchans", 1);
    series.key_i32("nbits", 8);
    series.key_f64("tsamp", 1.0);
    series.key_f64("fch1", 1400.0);
    series.key_f64("foff", -1.0);
    series.string("HEADER_END");
    REQUIRE_THROWS_WITH(FilterbankHeader::parse(series.out).validate(),
                        ContainsSubstring("separate reader"));
}

TEST_CASE("unlisted telescope ids are unknown", "[sigproc][header]") {
    REQUIRE(psrio::formats::sigproc::telescope_name(999) == "unknown");
    REQUIRE(psrio::formats::sigproc::data_type_name(2) == "time series");
    REQUIRE(psrio::formats::sigproc::telescope_name(21) == "FAST");
    REQUIRE(psrio::formats::sigproc::telescope_name(65) == "KAT-7");
    REQUIRE(psrio::formats::sigproc::machine_name(83) == "ROACH");
}

TEST_CASE("unknown keys with string values are skipped and preserved",
          "[sigproc][header]") {
    SigprocBytes bytes;
    bytes.string("HEADER_START");
    bytes.key_i32("nchans", 2);
    bytes.key_i32("nbits", 8);
    bytes.key_f64("tsamp", 1.0);
    bytes.key_f64("fch1", 1400.0);
    bytes.key_f64("foff", -1.0);

    // Custom unknown string key
    bytes.string("observer_name");
    bytes.string("Antigravity");

    bytes.string("HEADER_END");

    const auto header = FilterbankHeader::parse(bytes.out);
    header.validate();
    REQUIRE(header.extra.size() == 1U);
    REQUIRE(header.extra.front().name == "observer_name");
    // Raw value bytes should be 4-byte length + string characters
    const auto raw = header.extra.front().raw;
    REQUIRE(raw.size() == 4U + 11U);
}

TEST_CASE("header astronomical and physical convenience helpers",
          "[sigproc][header]") {
    SigprocBytes bytes;
    bytes.string("HEADER_START");
    bytes.key_i32("nchans", 4);
    bytes.key_i32("nbits", 8);
    bytes.key_f64("tsamp", 0.001); // 1 ms
    bytes.key_f64("tstart", 59000.25);
    bytes.key_f64("fch1", 1500.0);
    bytes.key_f64("foff", -10.0);
    bytes.key_f64("src_raj", 183045.5);
    bytes.key_f64("src_dej", -043015.0);
    bytes.string("HEADER_END");
    for (int i = 0; i < 40; ++i) { // 10 samples * 4 chans
        bytes.u8(1);
    }

    const auto header = FilterbankHeader::parse(bytes.out);
    header.validate();
    REQUIRE(header.nsamples() == 10U);

    // RA/Dec formatting & conversions
    REQUIRE(header.ra_string() == "18:30:45.5000");
    REQUIRE(header.dec_string() == "-04:30:15.0000");
    REQUIRE(header.ra_hours().has_value());
    REQUIRE(header.ra_degrees().has_value());
    REQUIRE(header.ra_radians().has_value());
    REQUIRE(header.dec_degrees().has_value());
    REQUIRE(header.dec_radians().has_value());

    // Observation duration & MJD
    REQUIRE(header.observation_duration() == 0.010);
    REQUIRE(header.duration_string() == "0.0 seconds");
    REQUIRE(header.gregorian_date() == "2020-05-31");
    REQUIRE(header.mjd_after_samples(1000) == 59000.25 + (1.0 / 86400.0));

    // Dispersion delays across the 4 channels (1500, 1490, 1480, 1470 MHz)
    std::array<double, 4> delays{};
    header.dispersion_delays(10.0, delays);
    // Highest freq channel (1500 MHz) has delay 0.0 relative to highest
    REQUIRE(delays[0] == 0.0);
    REQUIRE(delays[1] > 0.0);
    REQUIRE(delays[2] > delays[1]);
    REQUIRE(delays[3] > delays[2]);
}
