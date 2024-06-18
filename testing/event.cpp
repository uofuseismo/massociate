#include <vector>
#include <chrono>
#include <cmath>
#include "massociate/event.hpp"
#include "massociate/waveformIdentifier.hpp"
#include "massociate/arrival.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

TEST_CASE("MAssociate::Event", "[event]")
{
    uint64_t evid{5923};
    double latitude{40.7608};
    double longitude{-111.8910};
    double depth{5000};
    double originTime{50};
    //double x{2};
    //double y{3};
    //double z{4};
    MAssociate::Event::Type type{MAssociate::Event::Type::Trigger};

    MAssociate::Arrival arrival;
    MAssociate::WaveformIdentifier waveid;
    MAssociate::Event event;
 
    REQUIRE_NOTHROW(event.setIdentifier(evid));
    REQUIRE(event.getIdentifier() == evid);

    REQUIRE_NOTHROW(event.setLatitude(latitude));
    REQUIRE(std::abs(event.getLatitude() - latitude) < 1.e-10);

    REQUIRE_NOTHROW(event.setLongitude(longitude));
    REQUIRE(std::abs(event.getLongitude() - (longitude + 360)) < 1.e-10);
    event.setLongitude(event.getLongitude());
    REQUIRE(std::abs(event.getLongitude() - (longitude + 360)) < 1.e-10);

    REQUIRE_NOTHROW(event.setDepth(depth));
    REQUIRE(std::abs(event.getDepth() - depth) < 1.e-10);

    REQUIRE_NOTHROW(event.setOriginTime(originTime));
    REQUIRE(std::abs(event.getOriginTime().count()*1.e-6 - originTime) < 1.e-10);

    //event.setXPosition(x);
    //REQUIRE(std::abs(event.getXPosition() - x) < 1.e-14);
    //event.setYPosition(y);
    //REQUIRE(std::abs(event.getYPosition() - y) < 1.e-14);
    //event.setZPosition(z);
    //REQUIRE(std::abs(event.getZPosition() - z) < 1.e-14);

    event.setType(type);
    REQUIRE(event.getType() == type);

    // Create some arrivals - check this is copy c'tor
    std::vector<MAssociate::Arrival> arrivals;
    waveid.setNetwork("UU");
    waveid.setStation("LKWY");
    waveid.setChannel("ENZ");
    waveid.setLocationCode("00");
    arrival.setIdentifier(0);
    arrival.setTime(originTime + 4);
    arrival.setWaveformIdentifier(waveid);
    arrival.setPhase("P");
    arrivals.push_back(arrival);
 
    waveid.setStation("VEC");
    arrival.setTime(originTime + 8);
    arrival.setWaveformIdentifier(waveid);
    arrival.setPhase("S");
    arrival.setIdentifier(1);
    arrivals.push_back(arrival);

    waveid.setStation("COY");
    arrival.setTime(originTime + 6);
    arrival.setWaveformIdentifier(waveid);
    arrival.setPhase("P");
    arrival.setIdentifier(2);
    arrivals.push_back(arrival);
    for (const auto &arrival : arrivals)
    {
        REQUIRE_NOTHROW(event.addArrival(arrival));
    }
    REQUIRE(event.getNumberOfArrivals() == 3);
    REQUIRE(event.getNumberOfPArrivals() == 2);
    // Make some arrivals we can't add
    auto badArrivals = arrivals;
    badArrivals[0].setTime(originTime - 1); // Arrival before origin time
    badArrivals[1].setPhase("P"); // P at S
    badArrivals[2].setPhase("S"); // S at P

    // Copy constructor
    MAssociate::Event eCopy(event);
    REQUIRE(eCopy.getIdentifier() == evid);
    REQUIRE(std::abs(eCopy.getLatitude() - latitude) < 1.e-10);
    REQUIRE(std::abs(eCopy.getLongitude()- (longitude + 360)) < 1.e-10);
    REQUIRE(std::abs(eCopy.getDepth() - depth) < 1.e-10);
    REQUIRE(std::abs(eCopy.getOriginTime().count()*1.e-6 - originTime) < 1.e-10); 
    //REQUIRE(std::abs(eCopy.getXPosition() - x) < 1.e-14);
    //REQUIRE(std::abs(eCopy.getYPosition() - y) < 1.e-14);
    //REQUIRE(std::abs(eCopy.getZPosition() - z) < 1.e-14);
    std::array<int, 3> permutation{0, 2, 1}; // Arrivals are sorted in increasing time
    auto arrivalsBack = eCopy.getArrivals(); 
    for (int ia = 0; ia < static_cast<int> (arrivals.size()); ++ia)
    {
        REQUIRE(arrivalsBack[ia].getWaveformIdentifier() ==
                arrivals[permutation[ia]].getWaveformIdentifier());
        REQUIRE(arrivalsBack[ia].getPhase() ==
                arrivals[permutation[ia]].getPhase());
        REQUIRE(arrivalsBack[ia].getTime() ==
                arrivals[permutation[ia]].getTime());
        REQUIRE(eCopy.canAddArrival(arrivals[ia], true) >= 0);
        REQUIRE(eCopy.canAddArrival(arrivals[ia], false) <= 0);
        REQUIRE(eCopy.canAddArrival(badArrivals[ia], true) < 0);
    }
    event.clearArrivals();
    REQUIRE(event.getNumberOfArrivals() == 0);
    event.clear();
    REQUIRE(event.getType() == MAssociate::Event::Type::Event);
}

