#include "psrio/astro.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <numbers>
#include <stdexcept>
#include <string>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("sexagesimal parser extracts fields accurately", "[astro]") {
    const auto coord = psrio::astro::parse_sexagesimal("12:34:56.789");
    REQUIRE_FALSE(coord.negative);
    REQUIRE(coord.major == 12);
    REQUIRE(coord.minutes == 34);
    REQUIRE_THAT(coord.seconds, WithinAbs(56.789, 1e-6));

    const auto neg = psrio::astro::parse_sexagesimal(" -00:30:15.5 ");
    REQUIRE(neg.negative);
    REQUIRE(neg.major == 0);
    REQUIRE(neg.minutes == 30);
    REQUIRE_THAT(neg.seconds, WithinAbs(15.5, 1e-6));

    const auto plus = psrio::astro::parse_sexagesimal("+45:00:00.0");
    REQUIRE_FALSE(plus.negative);
    REQUIRE(plus.major == 45);
    REQUIRE(plus.minutes == 0);
    REQUIRE_THAT(plus.seconds, WithinAbs(0.0, 1e-6));
}

TEST_CASE("sexagesimal parser rejects malformed inputs", "[astro]") {
    REQUIRE_THROWS_AS(psrio::astro::parse_sexagesimal(""),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(psrio::astro::parse_sexagesimal("12:34"),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(psrio::astro::parse_sexagesimal("12:60:00.0"),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(psrio::astro::parse_sexagesimal("12:00:60.0"),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(psrio::astro::parse_sexagesimal("12:34:56.78extra"),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(psrio::astro::parse_sexagesimal("12:34:nan"),
                      std::invalid_argument);
}

TEST_CASE("coordinate string conversion and formatting", "[astro]") {
    REQUIRE(psrio::astro::ra_to_string(123456.78) == "12:34:56.7800");
    REQUIRE(psrio::astro::dec_to_string(-451234.5) == "-45:12:34.5000");
    REQUIRE(psrio::astro::dec_to_string(451234.5) == "+45:12:34.5000");

    // Edge case: Negative declination between 0 and -1 degree
    REQUIRE(psrio::astro::dec_to_string(-001234.5) == "-00:12:34.5000");
}

TEST_CASE("packed coordinate decimal and radian conversions", "[astro]") {
    // 12h 00m 00s = 180 degrees = pi radians
    REQUIRE_THAT(psrio::astro::ra_to_hours(120000.0), WithinAbs(12.0, 1e-9));
    REQUIRE_THAT(psrio::astro::ra_to_degrees(120000.0), WithinAbs(180.0, 1e-9));
    REQUIRE_THAT(psrio::astro::ra_to_radians(120000.0),
                 WithinAbs(std::numbers::pi_v<double>, 1e-9));

    // Dec: +30 00 00 = +30 degrees
    REQUIRE_THAT(psrio::astro::dec_to_degrees(300000.0), WithinAbs(30.0, 1e-9));
    REQUIRE_THAT(psrio::astro::dec_to_radians(300000.0),
                 WithinAbs(std::numbers::pi_v<double> / 6.0, 1e-9));

    // Dec: -30 00 00 = -30 degrees
    REQUIRE_THAT(psrio::astro::dec_to_degrees(-300000.0),
                 WithinAbs(-30.0, 1e-9));
    REQUIRE_THAT(psrio::astro::dec_to_radians(-300000.0),
                 WithinAbs(-std::numbers::pi_v<double> / 6.0, 1e-9));

    // Dec edge case: -00 30 00 = -0.5 degrees
    REQUIRE_THAT(psrio::astro::dec_to_degrees(-003000.0),
                 WithinAbs(-0.5, 1e-9));
}

TEST_CASE("string to radians conversion", "[astro]") {
    REQUIRE_THAT(psrio::astro::ra_to_rad("12:00:00.0"),
                 WithinAbs(std::numbers::pi_v<double>, 1e-9));
    REQUIRE_THAT(psrio::astro::dec_to_rad("+90:00:00.0"),
                 WithinAbs(std::numbers::pi_v<double> / 2.0, 1e-9));
    REQUIRE_THAT(psrio::astro::dec_to_rad("-90:00:00.0"),
                 WithinAbs(-std::numbers::pi_v<double> / 2.0, 1e-9));

    REQUIRE_THROWS_AS(psrio::astro::ra_to_rad("24:00:00.0"),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(psrio::astro::ra_to_rad("-01:00:00.0"),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(psrio::astro::dec_to_rad("+90:00:01.0"),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(psrio::astro::dec_to_rad("-91:00:00.0"),
                      std::invalid_argument);
}

TEST_CASE("deg_to_dms round-trip", "[astro]") {
    const double angle = 45.5; // 45 deg 30 min 00 sec
    REQUIRE_THAT(psrio::astro::deg_to_dms(angle), WithinAbs(453000.0, 1e-6));
    REQUIRE_THAT(psrio::astro::deg_to_dms(-angle), WithinAbs(-453000.0, 1e-6));
}

TEST_CASE("MJD to Gregorian calendar date", "[astro]") {
    REQUIRE(psrio::astro::mjd_to_time(40587.0) == 0);
    REQUIRE(psrio::astro::mjd_to_time(40587.5) == 43200);
    REQUIRE(psrio::astro::mjd_to_time(58000.0) == 1504483200);
    REQUIRE(psrio::astro::mjd_to_time(58000.0 + (1.5 / 86400.0)) == 1504483201);

    REQUIRE(psrio::astro::mjd_to_gregorian(50000) == "1995-10-10");
    REQUIRE(psrio::astro::mjd_to_gregorian(59000) == "2020-05-31");
    REQUIRE(psrio::astro::mjd_to_gregorian(59000.75) == "2020-05-31");
}

TEST_CASE("duration formatting", "[astro]") {
    REQUIRE(psrio::astro::format_duration(45.2) == "45.2 seconds");
    REQUIRE(psrio::astro::format_duration(150.0) == "2.5 minutes");
    REQUIRE(psrio::astro::format_duration(7200.0) == "2.0 hours");
    REQUIRE(psrio::astro::format_duration(172800.0) == "2.0 days");
}

TEST_CASE("dispersion delay calculation", "[astro]") {
    constexpr double kDm = 100.0;
    // Low freq: 1000 MHz, High freq: 2000 MHz
    // dt = 4.148808e3 * 100 * (1/1000^2 - 1/2000^2)
    //    = 414880.8 * (1e-6 - 0.25e-6) = 414880.8 * 0.75e-6 = 0.3111606 s
    const double delay = psrio::astro::dispersion_delay(kDm, 1000.0, 2000.0);
    REQUIRE_THAT(delay, WithinRel(0.3111606, 1e-5));

    // Infinite reference frequency
    const double delay_inf = psrio::astro::dispersion_delay(kDm, 1000.0, 0.0);
    REQUIRE_THAT(delay_inf, WithinRel(0.4148808, 1e-5));
}
