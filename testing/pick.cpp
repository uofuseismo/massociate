#include <cmath>
#include "massociate/pick.hpp"
#include "massociate/waveformIdentifier.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

TEST_CASE("MAssociate::Pick", "[pick]")
{
    MAssociate::Pick pick;
    constexpr std::chrono::microseconds time{1000000};
    constexpr double standardError{0.2};
    constexpr double weight{1/standardError};
    constexpr double firstMotionWeight{0.5};
    constexpr uint64_t id{66823};
    std::string phase{"P"};
    auto firstMotion = MAssociate::Arrival::FirstMotion::Up;
    MAssociate::WaveformIdentifier sncl;
    std::string network{"UU"};
    std::string station{"CAPU"};
    std::string channel{"EHZ"};
    std::string location{"01"};
    sncl.setNetwork(network);
    sncl.setStation(station);
    sncl.setChannel(channel);
    sncl.setLocationCode(location);

    pick.setWaveformIdentifier(sncl);
    pick.setPhaseHint(MAssociate::Pick::PhaseHint::P);
    pick.setTime(time);
    pick.setFirstMotion(firstMotion);
    pick.setStandardError(standardError);
    pick.setIdentifier(id);
    pick.setFirstMotionWeight(firstMotionWeight);

    MAssociate::Pick copyPick{pick};
    REQUIRE(copyPick.getWaveformIdentifier() == sncl);
    REQUIRE(copyPick.getPhaseHint() == phase);
    REQUIRE(copyPick.getTime() == time);
    REQUIRE(copyPick.getFirstMotion() == firstMotion);
    REQUIRE(std::abs(copyPick.getStandardError()
                   - standardError) < 1.e-14);
    REQUIRE(std::abs(copyPick.getWeight() - weight) < 1.e-14);
    REQUIRE(copyPick.getIdentifier() == id);
    REQUIRE(std::abs(copyPick.getFirstMotionWeight()
                   - firstMotionWeight) < 1.e-14);
}

