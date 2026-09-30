#include <catch2/catch_test_macros.hpp>

#include <psrio/psrio.hpp>

TEST_CASE("library version is exposed", "[sigproc]")
{
    REQUIRE(psrio::version() == "0.1.0");
}
