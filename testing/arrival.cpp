#include <cmath>
#include "massociate/arrival.hpp"
#include "massociate/pick.hpp"
#include "massociate/waveformIdentifier.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

TEST_CASE("MAssociate::Arrival", "[arrival]")
{
    MAssociate::Pick pick;
    constexpr std::chrono::microseconds time{1000000};
    constexpr double standardError{0.2};
    constexpr double travelTime{4};
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

    // Create an arrival from a pick
    MAssociate::Arrival copyArrival{pick};
    copyArrival.setTravelTime(travelTime);

    REQUIRE(copyArrival.getWaveformIdentifier() == sncl);
    REQUIRE(copyArrival.getPhase() == phase);
    REQUIRE(copyArrival.getTime() == time);
    REQUIRE(copyArrival.getFirstMotion() == firstMotion);
    REQUIRE(std::abs(copyArrival.getWeight() - weight) < 1.e-14);
    REQUIRE(std::abs(copyArrival.getStandardError() - standardError) < 1.e-14);
    REQUIRE(copyArrival.getIdentifier() == id);
    REQUIRE(std::abs(copyArrival.getTravelTime() - travelTime) < 1.e-14);
    REQUIRE(std::abs(copyArrival.getFirstMotionWeight() - firstMotionWeight) < 1.e-14);
}

