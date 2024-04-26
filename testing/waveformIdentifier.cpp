#include <string>
#include "massociate/waveformIdentifier.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

TEST_CASE("MAssociate", "[WaveformIdentifier]")
{
    MAssociate::WaveformIdentifier waveid;
    std::string network{"UU"};
    std::string station{"CAPU"};
    std::string channel{"EHZ"};
    std::string location{"01"};

    REQUIRE(waveid.isEmpty());

    waveid.setNetwork(network);
    REQUIRE(!waveid.isEmpty());
    REQUIRE(waveid.getNetwork() == network);
    waveid.setStation(station);
    REQUIRE(waveid.getStation() == station);
    waveid.setChannel(channel);
    REQUIRE(waveid.getChannel() == channel);
    waveid.setLocationCode(location);
    REQUIRE(waveid.getLocationCode() == location);

    MAssociate::WaveformIdentifier waveidCopy{waveid};
    REQUIRE(waveidCopy.getNetwork() == network);
    REQUIRE(waveidCopy.getStation() == station);
    REQUIRE(waveidCopy.getChannel() == channel);
    REQUIRE(waveidCopy.getLocationCode() == location);

    REQUIRE(waveidCopy == waveid);

    waveid.clear();
    REQUIRE(waveidCopy != waveid);
    REQUIRE(waveid.isEmpty());

    MAssociate::WaveformIdentifier waveid2(network, station, channel, location);
    REQUIRE(waveid2 == waveidCopy);
}
