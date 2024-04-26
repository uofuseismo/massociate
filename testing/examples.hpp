#include <string>
#ifndef NDEBUG
#include <cassert>
#endif
#include <uLocator/station.hpp>
#include <uLocator/travelTimeCalculatorMap.hpp>
#include <uLocator/uussRayTracer.hpp>
#include <uLocator/position/utahRegion.hpp>
#include <uLocator/position/ynpRegion.hpp>
#include <uLocator/position/wgs84.hpp>
#include "massociate/waveformIdentifier.hpp"
#include "massociate/arrival.hpp"
#include "massociate/pick.hpp"

namespace
{

MAssociate::Pick arrivalToPick(const MAssociate::Arrival &arrival)
{
    MAssociate::Pick pick;

    pick.setIdentifier(arrival.getIdentifier());
    pick.setWaveformIdentifier(arrival.getWaveformIdentifier());
    pick.setTime(arrival.getTime());
    pick.setStandardError(arrival.getStandardError());
    pick.setPhaseHint(arrival.getPhase());
    return pick;
}

std::pair<MAssociate::Arrival, ULocator::Station>
    toArrival(const int64_t identifier,
              const std::string &network,
                              const std::string &station,
                              const std::string &phase,
                              const double time,
                              const double standardError2,
                              const double stationLatitude,
                              const double stationLongitude,
                              const double stationElevation,
                              const bool isUtah = true,
                              const int utmZone = 12) 
{
    MAssociate::WaveformIdentifier waveformIdentifier;
    waveformIdentifier.setNetwork(network);
    waveformIdentifier.setStation(station);
    MAssociate::Arrival arrival;
    arrival.setIdentifier(identifier);
    arrival.setWaveformIdentifier(waveformIdentifier);
    arrival.setPhase(phase);
    arrival.setTime(time);
    arrival.setStandardError(standardError2);

    ULocator::Station uStation;
    uStation.setNetwork(network);
    uStation.setName(station);
    uStation.setElevation(stationElevation);
    ULocator::Position::WGS84 stationLocation{stationLatitude, stationLongitude, utmZone};
    if (isUtah)
    {   
        uStation.setGeographicPosition(stationLocation,
                                       ULocator::Position::UtahRegion {});
    }   
    else
    {   
        uStation.setGeographicPosition(stationLocation,
                                       ULocator::Position::YNPRegion {});
    }   
    return std::pair {arrival, uStation};
}

struct CreateTestCase60557072
{
    CreateTestCase60557072(bool addNoisePicks = false)
    {   
        auto [a1,  s1]  = ::toArrival(1,  "UU", "KNB",  "P", 1703612966.086906, 0.066805, 37.0166,-112.822,  1715.0, true, 12);
        auto [a2,  s2]  = ::toArrival(2,  "AE", "U15A", "P", 1703612967.977675, 0.102402, 36.428, -112.2915, 2489.0, true, 12);
        auto [a3,  s3]  = ::toArrival(3,  "UU", "LCMT", "P", 1703612971.257604, 0.065439, 37.0118,-113.2439, 1411.0, true, 12);
        auto [a4,  s4]  = ::toArrival(4,  "UU", "KNB",  "S", 1703612971.659966, 0.131456, 37.0166,-112.822,  1715.0, true, 12);
        auto [a5,  s5]  = ::toArrival(5,  "UU", "ZNPU", "P", 1703612972.431031, 0.072258, 37.3561,-113.1254, 1953.0, true, 12);
        auto [a6,  s6]  = ::toArrival(6,  "AE", "U15A", "S", 1703612974.938671, 0.185014, 36.428, -112.2915, 2489.0, true, 12);
        auto [a7,  s7]  = ::toArrival(7,  "UU", "SZCU", "P", 1703612975.187562, 0.312345, 37.5954,-113.0875, 2026.0, true, 12);
        auto [a8,  s8]  = ::toArrival(8,  "UU", "LCMT", "S", 1703612980.490127, 0.196639, 37.0118,-113.2439, 1411.0, true, 12);
        auto [a9,  s9]  = ::toArrival(9,  "UU", "ZNPU", "S", 1703612982.488189, 0.30128,  37.3561,-113.1254, 1953.0, true, 12);
        auto [a10, s10] = ::toArrival(10, "UU", "BHU",  "S", 1703612985.001897, 0.200406, 37.5939,-112.8621, 3250.0, true, 12);
        auto [a11, s11] = ::toArrival(11, "UU", "SZCU", "S", 1703612986.842699, 0.192562, 37.5954,-113.0875, 2026.0, true, 12);
        auto [a12, s12] = ::toArrival(12, "US", "WUAZ", "P", 1703612986.500783, 0.26241,  35.5169,-111.3739, 1592.0, true, 12); // Big residual so I moved it
        arrivals = std::vector {a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12};
        stations = std::vector {s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12};
        if (addNoisePicks)
        {
            auto [a13, s13] = ::toArrival(13, "UU", "KNB",  "P", 1703612936.086906, 0.066805, 37.0166,-112.822,  1715.0, true, 12);
            auto [a14, s14] = ::toArrival(14, "AE", "U15A", "S", 1703612901.938671, 0.185014, 36.428, -112.2915, 2489.0, true, 12);
            auto [a15, s15] = ::toArrival(15, "US", "WUAZ", "P", 1703612917.500783, 0.26241,  35.5169,-111.3739, 1592.0, true, 12);
            arrivals.push_back(a13);
            arrivals.push_back(a14);
            arrivals.push_back(a15);
            stations.push_back(s13);
            stations.push_back(s14);
            stations.push_back(s15);
        }
        for (int i = 0; i < static_cast<int> (arrivals.size()); ++i)
        {
            auto [xSource, ySource] = region.geographicToLocalCoordinates(eventLatitude, eventLongitude);
            constexpr bool applyCorrection{true};
            if (arrivals.at(i).getPhase() == "P")
            {
                if (!travelTimeCalculatorMap->contains(stations[i], "P"))
                {
                    std::unique_ptr<ULocator::ITravelTimeCalculator> pCalculator
                       = std::make_unique<ULocator::UUSSRayTracer>
                         (stations.at(i),
                          ULocator::UUSSRayTracer::Phase::P,
                          ULocator::UUSSRayTracer::Region::Utah);
                    travelTimes.push_back(pCalculator->evaluate(
                        0.0, xSource, ySource, eventDepth, applyCorrection));
                    travelTimeCalculatorMap->insert(stations[i], "P",
                                                    std::move(pCalculator));
                }   
            }
            else
            {
                if (!travelTimeCalculatorMap->contains(stations[i], "S"))
                {
                    std::unique_ptr<ULocator::ITravelTimeCalculator> sCalculator
                       = std::make_unique<ULocator::UUSSRayTracer>
                         (stations.at(i),
                          ULocator::UUSSRayTracer::Phase::S,
                          ULocator::UUSSRayTracer::Region::Utah);
                    travelTimes.push_back(sCalculator->evaluate(
                        0.0, xSource, ySource, eventDepth, applyCorrection));
                    travelTimeCalculatorMap->insert(stations[i], "S",
                                                    std::move(sCalculator));
                }   
            }   
        }
#ifndef NDEBUG
        assert(arrivals.size() == travelTimes.size());
#endif
    }
    std::vector<ULocator::Station> stations;
    std::vector<MAssociate::Arrival> arrivals;
    std::vector<double> travelTimes;
    std::unique_ptr<ULocator::TravelTimeCalculatorMap>
        travelTimeCalculatorMap{
            std::make_unique<ULocator::TravelTimeCalculatorMap> ()};
    ULocator::Position::UtahRegion region;  
    int64_t eventIdentifier{60557072};
    double eventLatitude{36.8758333};
    double eventLongitude{-112.4343333};
    double eventDepth{20200.0};
    double eventTime{1703612958.6499996};
};

}
