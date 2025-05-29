#include <iostream>
#include <iomanip>
#include <set>
#include <string>
#include <vector>
#include <map>
#include <uLocator/station.hpp>
#include <uLocator/optimizers/originTime.hpp>
#include <uLocator/travelTimeCalculatorMap.hpp>
#include <uLocator/travelTimeCalculator.hpp>
#include <uLocator/position/knownLocalLocation.hpp>
#include <uLocator/position/geographicRegion.hpp>
#ifdef WITH_LOCAL_UMPS_LOGGING
#include "logging/standardOut.hpp"
#else
#include <umps/logging/standardOut.hpp>
#endif
#include "massociate/migrator.hpp"
#include "massociate/arrival.hpp"
#include "massociate/waveformIdentifier.hpp"
#include "analyticSignals.hpp"
#include "knn.hpp"

using namespace MAssociate; 

namespace
{
double truncatedGaussianNormalization(const double maxStandardDeviations)
{
    double PhiLower = 0.5*(1 + std::erf(-maxStandardDeviations*M_SQRT1_2));
    double PhiUpper = 0.5*(1 + std::erf( maxStandardDeviations*M_SQRT1_2));
    double PhiNormalization = 1./(PhiUpper - PhiLower);
    const double gaussianNormalization = 1./std::sqrt(2*M_PI);
    return gaussianNormalization*PhiNormalization;
}
double truncatedGaussianContribution(
    const double observation,
    const double estimate,
    const double standardDeviation,
    const double normalization, // 1/[ sqrt(2*pi)) x (Phi_u - Phi_l) ]
    const double maxStandardDeviations,
    const bool acceptAnyway = false)
{
    auto residual = observation - estimate;
    double likelihood = 0; 
    if (acceptAnyway ||
        std::abs(residual) < maxStandardDeviations*standardDeviation)
    {
        double xi = residual/standardDeviation;
        likelihood = (normalization/standardDeviation)*std::exp(-0.5*(xi*xi));
    }
    return likelihood;
}
}

class IMigrator::IMigratorImpl
{
public:
    IMigratorImpl(std::shared_ptr<UMPS::Logging::ILog> logger = nullptr) :
        mLogger(logger)
    {
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> ();
        }
    }
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    std::unique_ptr<ULocator::TravelTimeCalculatorMap>
        mTravelTimeCalculatorMap{nullptr};
    std::unique_ptr<ULocator::Position::IGeographicRegion>
        mGeographicRegion{nullptr};
    std::vector<Arrival> mArrivals;
    std::vector<std::unique_ptr<ULocator::Position::IKnownLocalLocation>> mKnownLocations;
    std::vector<std::pair<Arrival, double>> mContributingArrivals;
    std::vector<std::pair<ULocator::Station, std::string>> mUniqueStationPhases;
    std::vector<double> mObservations;
    std::vector<double> mStandardDeviations;
    std::vector<int> mArrivalToUniqueStationPhaseMap;
    std::set<std::set<int>> mFeasibleArrivalSets;
    std::vector<std::pair<double, double>> mDifferentialArrivalPairs;
    std::vector<std::pair<double, double>>
         mDifferentialArrivalPairsStandardDeviations;
    std::vector<std::pair<int, int>> mDifferentialArrivalIndicesPairs;
    std::vector<std::pair<int, int>> mDifferentialTablePairs;
    std::map<std::string, std::vector<std::string>> mStationNeighbors;
    double mMaximumStandardDeviation{3};
    double mMaximumEpicentralDistance{1.e30};
    IMigrator::SignalType
        mSignalType{IMigrator::SignalType::DoubleDifference};
    IMigrator::PickSignal
        mPickSignal{IMigrator::PickSignal::Boxcar};
    bool mHaveContributingArrivals{false};
    bool mSaveArrivalContributionList{false};
};

/// Constructor
IMigrator::IMigrator() :
    pImpl(std::make_unique<IMigratorImpl> (nullptr))
{
}

/// Constructor
IMigrator::IMigrator(std::shared_ptr<UMPS::Logging::ILog> &logger) :
    pImpl(std::make_unique<IMigratorImpl> (logger))
{
}

/// Move constructor
IMigrator::IMigrator(IMigrator &&migrator) noexcept
{
    *this = std::move(migrator);
}

/// Move assignment
IMigrator& IMigrator::operator=(IMigrator &&migrator) noexcept
{
    if (&migrator == this){return *this;}
    pImpl = std::move(migrator.pImpl);
    return *this;
}

/// Geographic region
void IMigrator::setGeographicRegion(
    const ULocator::Position::IGeographicRegion &region)
{
    pImpl->mGeographicRegion = region.clone();
}

std::unique_ptr<ULocator::Position::IGeographicRegion>
IMigrator::getGeographicRegion() const
{
    if (!haveGeographicRegion())
    {
        throw std::runtime_error("Geographic region not set");
    }
    return pImpl->mGeographicRegion->clone();
}

bool IMigrator::haveGeographicRegion() const noexcept
{
    return pImpl->mGeographicRegion != nullptr;
}

/// Signal type
void IMigrator::setSignalType(const IMigrator::SignalType signalType) noexcept
{
    pImpl->mSignalType = signalType;
}

IMigrator::SignalType IMigrator::getSignalType() const noexcept
{
    return pImpl->mSignalType;
}

/// Signal to migrate
void IMigrator::setPickSignalToMigrate(const PickSignal pickSignal) noexcept
{
    pImpl->mPickSignal = pickSignal;
}

IMigrator::PickSignal IMigrator::getPickSignalToMigrate() const noexcept
{
    return pImpl->mPickSignal;
}

/// Travel time calculator map
void IMigrator::setTravelTimeCalculatorMap(
    std::unique_ptr<ULocator::TravelTimeCalculatorMap> &&map,
    int numberOfNeighbors)
{
    if (map == nullptr)
    {   
        throw std::invalid_argument("Travel time calculator map is NULL");
    }   
    pImpl->mStationNeighbors.clear();
    pImpl->mTravelTimeCalculatorMap = std::move(map);
    // Figure out each station's neighbors
    auto stationPhasePairs = pImpl->mTravelTimeCalculatorMap->getStationPhasePairs();
    std::vector<std::string> stationNames;
    std::vector<double> xs;
    std::vector<double> ys;
    for (const auto &stationPhasePair : stationPhasePairs)
    {
        auto stationName = stationPhasePair.first.getNetwork() + "."
                         + stationPhasePair.first.getName();
        auto [x, y] = stationPhasePair.first.getLocalCoordinates();
        if (std::find(stationNames.begin(),
                      stationNames.end(),
                      stationName) == stationNames.end())
        {
            xs.push_back(x);
            ys.push_back(y); 
            stationNames.push_back(stationName);
        }
    }
    if (numberOfNeighbors > 0)
    {
        pImpl->mStationNeighbors = ::cluster(xs, ys, stationNames, numberOfNeighbors);
    }
}

const ULocator::TravelTimeCalculatorMap 
    *IMigrator::getTravelTimeCalculatorMap() const
{
    if (!haveTravelTimeCalculatorMap())
    {   
        throw std::runtime_error("Travel time calculator map not set");
    }   
    return pImpl->mTravelTimeCalculatorMap.get();
}

std::unique_ptr<ULocator::TravelTimeCalculatorMap>
IMigrator::releaseTravelTimeCalculatorMap()
{
    if (!haveTravelTimeCalculatorMap())
    {
        throw std::runtime_error("Travel time calculator map not set");
    }
    std::unique_ptr<ULocator::TravelTimeCalculatorMap>
        result{pImpl->mTravelTimeCalculatorMap.release()};
    pImpl->mTravelTimeCalculatorMap = nullptr;
    return result;
}

bool IMigrator::haveTravelTimeCalculatorMap() const noexcept
{
    return (pImpl->mTravelTimeCalculatorMap != nullptr);
}

/// Insert a travel time calculator
void IMigrator::addTravelTimeCalculator(
    const ULocator::Station &station,
    const std::string &phase,
    std::unique_ptr<const ULocator::ITravelTimeCalculator> &&calculator)
{
    if (calculator == nullptr)
    {
        throw std::invalid_argument("Calculator is NULL");
    }
    if (!station.haveNetwork())
    {
        throw std::invalid_argument("Network not set");
    }
    if (!station.haveName())
    {
        throw std::invalid_argument("Station not set");
    }
    if (phase.empty()){throw std::invalid_argument("Phase not defined");}
    if (!haveTravelTimeCalculatorMap())
    {
        pImpl->mLogger->debug("No map exists - creating new one...");
        auto travelTimeCalculatorMap
            = std::make_unique<ULocator::TravelTimeCalculatorMap> ();
    }
#ifndef NDEBUG
    assert(haveTravelTimeCalculatorMap());
#endif
    auto name = station.getNetwork() + "."
              + station.getName() + "."
              + phase;
    if (!pImpl->mTravelTimeCalculatorMap->contains(station, phase))
    {
        pImpl->mLogger->debug("Inserting " + name + "...");
        pImpl->mTravelTimeCalculatorMap->insert(station, phase,
                                                std::move(calculator));
    }
    else
    {
        pImpl->mLogger->error("Could not add " + name);
        throw std::runtime_error("Calculator already exists for " + name);
    }
}


/// Sets the arrivals
void IMigrator::setArrivals(const std::vector<Arrival> &arrivals)
{
    if (arrivals.empty())
    {
        throw std::invalid_argument("No arrivals!");
    }   
    if (!haveTravelTimeCalculatorMap())
    {
        throw std::runtime_error("Travel time calculator map not set");
    }
    pImpl->mArrivals.clear();
    pImpl->mUniqueStationPhases.clear();
    pImpl->mObservations.clear();
    pImpl->mStandardDeviations.clear();
    pImpl->mArrivalToUniqueStationPhaseMap.clear();
    pImpl->mContributingArrivals.clear();
    pImpl->mHaveContributingArrivals = false;

    pImpl->mArrivals.reserve(arrivals.size());
    pImpl->mUniqueStationPhases.reserve(arrivals.size());
    pImpl->mObservations.reserve(arrivals.size());
    pImpl->mStandardDeviations.reserve(arrivals.size());
    pImpl->mArrivalToUniqueStationPhaseMap.reserve(arrivals.size());

    pImpl->mDifferentialArrivalPairs.clear();
    pImpl->mDifferentialArrivalPairsStandardDeviations.clear();
    pImpl->mDifferentialArrivalIndicesPairs.clear();
    pImpl->mDifferentialTablePairs.clear();

    pImpl->mFeasibleArrivalSets.clear();

    std::vector<std::string> namePhases;
    namePhases.reserve(arrivals.size());
    for (const auto &arrival : arrivals)
    {
        if (!arrival.haveWaveformIdentifier())
        {
            pImpl->mLogger->warn(
                "Arrival's waveform identifier not set.  Skipping...");
            continue;
        }
        auto waveformIdentifier = arrival.getWaveformIdentifier();
        auto name = waveformIdentifier.getNetwork() + "."
                  + waveformIdentifier.getStation();
        if (!arrival.havePhase())
        {
            pImpl->mLogger->warn("Arrival's phase not set on " + name
                               + ".  Skipping...");
            continue;
        }
        auto phase = arrival.getPhase();
        auto namePhase = name + "." + phase;
        if (!arrival.haveTime())
        {
            pImpl->mLogger->warn("Arrival time not set on "
                               + namePhase + ".  Skipping...");
            continue;
        }
 /*
        if (!arrival.haveStandardDeviation())
        {
            pImpl->mLogger->warn("Standard error not set on "
                              + namePhase + ".  Skipping...");
            continue;
        }
*/
        ULocator::Station station;
        station.setNetwork(waveformIdentifier.getNetwork());
        station.setName(waveformIdentifier.getStation());
        if (!pImpl->mTravelTimeCalculatorMap->contains(station, phase))
        {
            pImpl->mLogger->warn("No travel time calculator for "
                               + namePhase
                               + ".  Skipping...");
            continue;
        }
        bool duplicate{false};
        int uniqueStationPhaseIndex{-1};
        auto indexPtr = std::find(namePhases.begin(), namePhases.end(),
                                  namePhase);
        if (indexPtr == namePhases.end())
        {
            duplicate = false; 
            uniqueStationPhaseIndex
                = static_cast<int> (namePhases.size());
        }
        else
        {
            duplicate = true;
            uniqueStationPhaseIndex 
               = static_cast<int> (std::distance(namePhases.begin(), indexPtr));
        }
        // Add everything
        pImpl->mArrivals.push_back(arrival);
        pImpl->mArrivalToUniqueStationPhaseMap.push_back(
            uniqueStationPhaseIndex);
        double arrivalTime = arrival.getTime().count()*1.e-6;
        pImpl->mObservations.push_back(arrivalTime);
        //reductionTime = std::min(arrivalTime, reductionTime);
        pImpl->mStandardDeviations.push_back(arrival.getStandardError());
        // For computational reasons keep a unique list 
        if (!duplicate)
        {
            pImpl->mUniqueStationPhases.push_back(std::pair {station, phase});
            namePhases.push_back(namePhase);
        }
    }
    if (pImpl->mArrivals.empty())
    {
        throw std::invalid_argument("No arrivals to set");
    }
    // Make the differential tables
    auto nArrivalsSet = static_cast<int> (pImpl->mArrivals.size());
    if (getSignalType() == IMigrator::SignalType::DoubleDifference)
    {
        auto spaceEstimate = std::max(1, (nArrivalsSet*(nArrivalsSet - 1))/2);
        pImpl->mDifferentialArrivalPairs.reserve(spaceEstimate);
        pImpl->mDifferentialArrivalPairsStandardDeviations.reserve(
            spaceEstimate);
        pImpl->mDifferentialArrivalIndicesPairs.reserve(spaceEstimate);
        pImpl->mDifferentialTablePairs.reserve(spaceEstimate);
        for (int i = 0; i < nArrivalsSet; ++i)
        {
            auto iUniqueStationPhase
                = pImpl->mArrivalToUniqueStationPhaseMap.at(i); 
            auto name1 = pImpl->mUniqueStationPhases.at(iUniqueStationPhase).first.getNetwork()
                       + "."
                       + pImpl->mUniqueStationPhases.at(iUniqueStationPhase).first.getName();
            auto phase1
                = pImpl->mUniqueStationPhases.at(iUniqueStationPhase).second;
            auto arrivalTime1 = pImpl->mObservations.at(i);
            auto std1 = pImpl->mStandardDeviations.at(i);
            const auto &neighborList = pImpl->mStationNeighbors.find(name1);
            // Do I even have enough neighbors?
            int nNeighboringArrivals{0};
            if (neighborList != pImpl->mStationNeighbors.end())
            {
                for (int j = i + 1; j < nArrivalsSet; ++j)
                {
                    auto jUniqueStationPhase
                        = pImpl->mArrivalToUniqueStationPhaseMap.at(j);
                    auto name2 = pImpl->mUniqueStationPhases.at(jUniqueStationPhase).first.getNetwork()
                               + "."
                               + pImpl->mUniqueStationPhases.at(jUniqueStationPhase).first.getName();
                    if (std::find(neighborList->second.begin(),
                                  neighborList->second.end(),
                                  name2) != neighborList->second.end())
                    {
                        nNeighboringArrivals = nNeighboringArrivals + 1;
                    }
                }
            }
            // Make the pairs
            for (int j = i + 1; j < nArrivalsSet; ++j)
            {
                auto jUniqueStationPhase
                    = pImpl->mArrivalToUniqueStationPhaseMap[j];
                auto name2 = pImpl->mUniqueStationPhases[jUniqueStationPhase].first.getNetwork()
                           + "."
                           + pImpl->mUniqueStationPhases[jUniqueStationPhase].first.getName();
                auto phase2
                    = pImpl->mUniqueStationPhases[jUniqueStationPhase].second;
                auto arrivalTime2 = pImpl->mObservations.at(j);
                auto std2 = pImpl->mStandardDeviations.at(j);
                if (phase2 != "P" && phase2 != "S"){continue;}
                if (name1 == name2)
                {
                    // Can't have two phases of same type
                    if (phase1 == phase2){continue;}
                    // Can't have S waves arrive before P waves
                    if (phase1 == "P" && phase2 == "S" && 
                        arrivalTime1 > arrivalTime2)
                    {
                        continue;
                    }
                    if (phase1 == "S" && phase2 == "P" &&
                        arrivalTime1 < arrivalTime2) 
                    {
                        continue;
                    }
                }
                // Is this one of my neighbors?
                if (neighborList != pImpl->mStationNeighbors.end())
                {
                    bool isNeighbor{false};
                    if (nNeighboringArrivals >= 1)
                    {
                        if (std::find(neighborList->second.begin(),
                                      neighborList->second.end(),
                                      name2) != neighborList->second.end())
                        {
                            isNeighbor = true;
                        }
                    }
                    if (!isNeighbor){continue;}
                }
                // Add it
                auto iTable = pImpl->mArrivalToUniqueStationPhaseMap.at(i);
                auto jTable = pImpl->mArrivalToUniqueStationPhaseMap.at(j);
#ifndef NDEBUG
                assert(iTable != jTable);
#endif
                pImpl->mDifferentialArrivalPairs.push_back(
                    std::pair {arrivalTime1, arrivalTime2});
                double std1Use = std1;
                double std2Use = std2;
                if (getPickSignalToMigrate() == IMigrator::PickSignal::Boxcar)
                {
                    // Standard errors come from Gaussian distributions.
                    // Instead, we want boxcar widths which are uniform 
                    // distributions.  This is really just solving for x where
                    // sigma = sqrt( (2x)^2 / 12) = 2x/sqrt(12) = x/sqrt(3)
                    // so the new sigma is sqrt(3)*sigma
                    constexpr double sqrt3{1.7320508075688772};
                    std1Use = sqrt3*std1;
                    std2Use = sqrt3*std2;
                }
                pImpl->mDifferentialArrivalPairsStandardDeviations.push_back(
                    std::pair {std1Use, std2Use});
                pImpl->mDifferentialArrivalIndicesPairs.push_back(
                    std::pair {i, j});
                pImpl->mDifferentialTablePairs.push_back(
                    std::pair {iTable, jTable});
            }
        }
        pImpl->mLogger->debug("Set " + std::to_string(nArrivalsSet)
                            + " phase arrivals which has "
                            + std::to_string(pImpl->mDifferentialTablePairs.size())
                            + " phase pairs");
    }
    else
    {
        // Loosely speaking - we have to compute the sets of feasible
        // arrival collections.  For example, if station UU.FORK has 
        // two P arrivals then we would have 2 sets.  Add in another
        // station, say UU.FORU that has 2 S arrivals then we now have
        // four feasible sets (each of the previous two feasible sets
        // has a feasible S arrival set for one FORU arrival and the 
        // the other FORU arrival).  My guess is the number of
        // feasible sets is the product of duplicate arrivals 
        // on each station.
        pImpl->mLogger->debug("Set " + std::to_string(nArrivalsSet)
                            + " phase arrivals");
        namePhases.clear();
        namePhases.reserve(nArrivalsSet);
        for (const auto &arrival : pImpl->mArrivals)
        {
            const auto &waveformIdentifier = arrival.getWaveformIdentifier();
            auto name = waveformIdentifier.getNetwork() + "." 
                      + waveformIdentifier.getStation() + "."
                      + arrival.getPhase();
            namePhases.push_back(std::move(name));
        }
        // Make the realizable combinations
        std::set<std::set<int>> arrivalSets;
        for (int iArrival = 0; iArrival < nArrivalsSet; ++iArrival)
        {
            std::set<int> initialArrivalSet{ iArrival };
            std::set<std::string> namePhasesInSet{ namePhases[iArrival] };
            // Now try to build the set
            for (int jArrival = 0; jArrival < nArrivalsSet; ++jArrival)
            {
                if (iArrival == jArrival){continue;} // Already have it
                if (!namePhasesInSet.contains(namePhases[jArrival]))
                {
                    initialArrivalSet.insert(jArrival);
                    namePhasesInSet.insert(namePhases[jArrival]);
                }
            }
//std::cout << "candidate set set: " << initialArrivalSet.size() << std::endl;
            // Does each arrival in the set have a neighbor?
            std::set<int> arrivalSet;
            for (const auto &candidateArrivalIndex : initialArrivalSet)
            {
                auto arrival = pImpl->mArrivals[candidateArrivalIndex];
                auto name = arrival.getWaveformIdentifier().getNetwork()
                          + "."
                          + arrival.getWaveformIdentifier().getStation();
                const auto &neighborList = pImpl->mStationNeighbors.find(name);
                // Look through the spatial neighbors and try to find a
                // corresponding pick
                if (neighborList != pImpl->mStationNeighbors.end())
                {
//for (const auto &neighbor : neighborList->second){std::cout<< neighbor << std::endl;}
//std::cout << std::endl;
                    bool hasNeighbor{false};
                    for (const auto &otherArrivalIndex : initialArrivalSet)
                    {
                        // This is me
                        if (candidateArrivalIndex == otherArrivalIndex)
                        {
                            continue;
                        }
                        // 
                        arrival = pImpl->mArrivals[otherArrivalIndex]; 
                        auto otherName
                            = arrival.getWaveformIdentifier().getNetwork()
                            + "."
                            + arrival.getWaveformIdentifier().getStation();
                        if (std::find(neighborList->second.begin(),
                                      neighborList->second.end(),
                                      otherName) != neighborList->second.end())
                        {
                            hasNeighbor = true;
                            break;
                        }
                    }
                    if (hasNeighbor)
                    {
                        arrivalSet.insert(candidateArrivalIndex);
                    }
                    else
                    {
                    }
                }
                else
                {
                    arrivalSet.insert(candidateArrivalIndex);
                }
            }
                          
            // Does the set already exist?
            if (!arrivalSets.contains(arrivalSet))
            {
//                std::cout << "--------------" << std::endl;
//                for (auto &idx : arrivalSet){std::cout << pImpl->mArrivals.at(idx).getIdentifier() << std::endl;}
//                std::cout << "--------------" << std::endl;
                if (arrivalSet.size() >= 4) // TODO hardwired
                {
                    pImpl->mLogger->debug(
                        "Adding candidate arrival set with "
                      + std::to_string(arrivalSet.size()) + " picks");
/*
                    std::cout << "--------------" << std::endl;
                    for (auto &idx : arrivalSet)
                    {
                        std::cout << pImpl->mArrivals.at(idx).getWaveformIdentifier().getNetwork() + "."
                                  << pImpl->mArrivals.at(idx).getWaveformIdentifier().getStation() + "."
                                  << pImpl->mArrivals.at(idx).getPhase() << std::endl;
                    }
                    std::cout << "--------------" << std::endl;
*/
                    arrivalSets.insert(std::move(arrivalSet));
                }
            }
        }
        pImpl->mFeasibleArrivalSets = std::move(arrivalSets);
    }
}

std::vector<Arrival> IMigrator::getArrivals() const
{
    if (!haveArrivals()){throw std::runtime_error("Arrivals not set");}
    return pImpl->mArrivals;
}

const std::vector<Arrival> &IMigrator::getArrivalsReference() const
{
    if (!haveArrivals()){throw std::runtime_error("Arrivals not set");}
    return *&pImpl->mArrivals;
}

int IMigrator::getNumberOfArrivals() const noexcept
{
    return static_cast<int> (pImpl->mArrivals.size());
}

bool IMigrator::haveArrivals() const noexcept
{
    return !pImpl->mArrivals.empty();
}

/// Save contributing picks?
void IMigrator::enableSaveArrivalContributionList() noexcept
{
     pImpl->mSaveArrivalContributionList = true;
}

void IMigrator::disableSaveArrivalContributionList() noexcept
{
    pImpl->mContributingArrivals.clear();
    pImpl->mHaveContributingArrivals = false;
    pImpl->mSaveArrivalContributionList = false;
}

bool IMigrator::saveArrivalContributionList() const noexcept
{
    return pImpl->mSaveArrivalContributionList;
}

void IMigrator::setDefaultSearchLocations(
    const std::vector<std::unique_ptr<ULocator::Position::IKnownLocalLocation>>
        &locations)
{
    pImpl->mKnownLocations.clear();
    for (const auto &location : locations)
    {
        if (location)
        {
            pImpl->mKnownLocations.push_back(location->clone());
        }
    }
}

std::unique_ptr<ULocator::Position::IKnownLocalLocation>
IMigrator::getKnownSearchLocation(size_t index) const
{
    return pImpl->mKnownLocations.at(index)->clone();
}

/// TODO: Introduce travel time database
/// Evaluate at the known locations
std::vector<std::pair<double, double>> IMigrator::evaluateAtKnownLocations() const
{
    auto nLocations = static_cast<int> (pImpl->mKnownLocations.size());
    std::vector<std::pair<double, double>> values(nLocations, std::pair {0.0, 0.0});
    if (getSignalType() == IMigrator::SignalType::DoubleDifference)
    {
        for (int i = 0; i < nLocations; ++i) 
        {
            if (pImpl->mKnownLocations[i])
            {
                constexpr double time{0}; // Doesn't matter
                try
                {
                    auto x = pImpl->mKnownLocations[i]->x();
                    auto y = pImpl->mKnownLocations[i]->y();
                    auto z = pImpl->mKnownLocations[i]->z();
                    constexpr double time{0};
                    double value = evaluate(x, y, z, time);
                    values[i] = std::pair{time, value};
                }
                catch (const std::exception &e)
                {    
                    pImpl->mLogger->warn("Failed to evaluate migration for "
                                        +       std::to_string(i)
                                        + ".  Failed with: "
                                        + std::string {e.what()});
                   values[i] = std::pair{0, 0}; // Contribute nothing
                } 
            } // End check on locations exists
        } // Loop on locations
    }
    else
    {
        ULocator::Optimizers::OriginTime originTimeCalculator;
        originTimeCalculator.setNorm(ULocator::Optimizers::IOptimizer::Norm::L1, 1);
        for (int i = 0; i < nLocations; ++i)
        {
            if (pImpl->mKnownLocations[i])
            {
                auto x = pImpl->mKnownLocations[i]->x();
                auto y = pImpl->mKnownLocations[i]->y();
                auto z = pImpl->mKnownLocations[i]->z();
                // Evaluate
                std::vector<double> uniqueTravelTimes;
                constexpr double originTime{0};
                constexpr bool applyCorrection{true};
                try
                {
                    pImpl->mTravelTimeCalculatorMap->evaluate(
                        pImpl->mUniqueStationPhases,
                        originTime,
                        x, y, z,
                        &uniqueTravelTimes,
                        applyCorrection);
                }
                catch (const std::exception &e)
                {
                    pImpl->mLogger->warn("Failed to compute travel time for location "
                                        +       std::to_string(i)
                                        + ".  Failed with: "
                                        + std::string {e.what()});
                    values[i] = std::pair{0, 0}; // Contribute nothing
                    continue;
                }
                // Solve for the optimal origin time in each feasible set
                double bestValue{0};
                for (const auto &feasibleArrivalSet :
                     pImpl->mFeasibleArrivalSets)
                {
                    auto nArrivals = static_cast<int> (feasibleArrivalSet.size());
                    std::vector<double> observations(nArrivals, 0);
                    std::vector<double> estimates(nArrivals, 0);
                    std::vector<double> weights(nArrivals, 0);
                    size_t index{0}; 
                    for (const auto &iArrival : feasibleArrivalSet)
                    {
                        // Get the arrival time
                        auto uniqueStationPhaseIndex
                            = pImpl->mArrivalToUniqueStationPhaseMap[iArrival];
                        estimates[index] = uniqueTravelTimes[uniqueStationPhaseIndex];
                        observations[index] = pImpl->mObservations[iArrival];
                        weights[index] = 1./pImpl->mStandardDeviations[iArrival];
                        index = index + 1;
                    }
                    try
                    {
                        originTimeCalculator.setArrivalTimes(observations,
                                                             weights);
                        originTimeCalculator.setTravelTimes(estimates);
                        originTimeCalculator.optimize();
                        auto originTimeForSet = originTimeCalculator.getTime();
                        auto valueForSet = evaluate(x, y, z, originTimeForSet);
                        if (valueForSet > bestValue)
                        {
                            bestValue = valueForSet;
                            values[i] = std::pair{originTimeForSet, valueForSet};
                        }
                    }
                    catch (const std::exception &e)
                    {
                        pImpl->mLogger->warn("Failed to evaluate migration for "
                                            +       std::to_string(i)
                                            + ".  Failed with: "
                                            + std::string {e.what()});
                    }
                } // Loop on candidate arrivals
            }
        } // Loop on locations
    } // End check on double difference / absolute
    return values;
}

/// Evaluate at a given point and depth
double IMigrator::evaluate(
    const double x, const double y, const double z, double time) const
{
    if (!haveArrivals()){throw std::runtime_error("Arrivals not set");}
    if (!haveTravelTimeCalculatorMap())
    {
        throw std::runtime_error("Travel time calculator not set");
    }
    auto maxDistance = getMaximumEpicentralDistance();
    auto saveContributingArrivals = saveArrivalContributionList();
    std::vector<double> cumulativeContributions;
    std::vector<bool> arrivalContributed;
    pImpl->mHaveContributingArrivals = false;
    if (saveContributingArrivals)
    {
        cumulativeContributions.resize(pImpl->mArrivals.size(), 0);
        arrivalContributed.resize(pImpl->mArrivals.size(), false);
    }
    auto pickSignal = getPickSignalToMigrate();
    std::vector<double> uniqueTravelTimes;
    double originTime{0};
    constexpr bool applyCorrection{true};
    if (getSignalType() == IMigrator::SignalType::Absolute){originTime = time;}
    pImpl->mTravelTimeCalculatorMap->evaluate(pImpl->mUniqueStationPhases,
                                              originTime,
                                              x, y, z,
                                              &uniqueTravelTimes,
                                              applyCorrection);
    auto uniqueDistances
        = pImpl->mTravelTimeCalculatorMap->computeDistance(
             pImpl->mUniqueStationPhases, x, y);
    double stack{0};
    std::set<int> bestSet;
    if (getSignalType() == IMigrator::SignalType::DoubleDifference)
    {
        if (pickSignal != IMigrator::PickSignal::Boxcar &&
            pickSignal != IMigrator::PickSignal::Gaussian)
        {
            throw std::runtime_error("Unhandled pick signal for double difference");
        }
#ifndef NDEBUG
        assert(pImpl->mDifferentialTablePairs.size() ==
               pImpl->mDifferentialArrivalPairs.size());
        assert(pImpl->mDifferentialTablePairs.size() ==
               pImpl->mDifferentialArrivalPairsStandardDeviations.size());
#endif
        auto nPairs
            = static_cast<int> (pImpl->mDifferentialTablePairs.size());
        for (int iPair = 0; iPair < nPairs; ++iPair)
        {
            auto iTable = pImpl->mDifferentialTablePairs[iPair].first;
            auto jTable = pImpl->mDifferentialTablePairs[iPair].second;
            auto distance1 = uniqueDistances[iTable];
            auto distance2 = uniqueDistances[jTable];
            if (distance1 < maxDistance && distance2 < maxDistance)
            {
                const auto &observedArrivalPair
                    = pImpl->mDifferentialArrivalPairs[iPair];
                auto tObserved1 = observedArrivalPair.first;
                auto tObserved2 = observedArrivalPair.second;
                const auto &observedArrivalPairWeight 
                    = pImpl->mDifferentialArrivalPairsStandardDeviations[iPair];
                auto std1 = observedArrivalPairWeight.first;
                auto std2 = observedArrivalPairWeight.second;
                auto tEstimated1 = uniqueTravelTimes[iTable];
                auto tEstimated2 = uniqueTravelTimes[jTable];
                auto dtEstimated = tEstimated2 - tEstimated1;
                double stackContributionIJ{0};
                if (pickSignal == IMigrator::PickSignal::Boxcar)
                {
                    stackContributionIJ
                        = ::analyticBoxcarCorrelation<double>(
                              tObserved1, tObserved2,
                              tEstimated1, tEstimated2,
                              std1, std2); 
                }
                else if (pickSignal == IMigrator::PickSignal::Gaussian)
                {
                    auto amplitudeNormalization
                        = ::analyticGaussianCorrelationAmplitudeNormalization(
                             std1, std2);
                    auto exponentNormalization
                        = ::analyticGaussianCorrelationExponentNormalization(
                             std1, std2);
                    stackContributionIJ
                          = ::analyticGaussianCorrelation(dtEstimated,
                                                          tObserved1,
                                                          tObserved2,
                                                          amplitudeNormalization,
                                                          exponentNormalization);
                }
                else
                {
                    throw std::runtime_error("Unhandled pick signal");
                }
                // Save the contribution?
                if (saveContributingArrivals)
                {
                    double maxContribution{0};
                    if (pickSignal == IMigrator::PickSignal::Boxcar)
                    {
                        maxContribution
                           = ::analyticBoxcarCorrelation<double>(0, 0, 0, 0,
                                                                 std1, std2);
                    }
                    else if (pickSignal == IMigrator::PickSignal::Gaussian)
                    {
                        auto amplitudeNormalization
                          = ::analyticGaussianCorrelationAmplitudeNormalization(
                              std1, std2);
                        auto exponentNormalization
                          = ::analyticGaussianCorrelationExponentNormalization(
                              std1, std2);
                        maxContribution
                          = ::analyticGaussianCorrelation<double>(0, 0, 0,
                                                          amplitudeNormalization,
                                                          exponentNormalization);

                    }
                    if (stackContributionIJ > maxContribution*1.e-12)
                    { 
                        auto arrivalIndexPair
                            = pImpl->mDifferentialArrivalIndicesPairs[iPair];
                        auto iArrival = arrivalIndexPair.first;
                        auto jArrival = arrivalIndexPair.second;
//std::cout << "add: " << iTable << " " << jTable << " " << iArrival << " " << jArrival << " " << stackContributionIJ << " " << maxContribution*1.e-12 << std::endl;
                        cumulativeContributions[iArrival]
                            = cumulativeContributions[iArrival]
                            + stackContributionIJ;
                        cumulativeContributions[jArrival]
                            = cumulativeContributions[jArrival]
                            + stackContributionIJ;
                        arrivalContributed[iArrival] = true;
                        arrivalContributed[jArrival] = true;
                    }
                }   
                // Finally update the stack
                stack = stack + stackContributionIJ;
            } // End check on distance
        } // Loop on pairs
    }
    else
    {
        // For each feasible arrival set tabulate 
        double normalization = ::truncatedGaussianNormalization(
            pImpl->mMaximumStandardDeviation);
        double bestStackContribution{-1};
        std::vector<bool> contributed;
        if (saveContributingArrivals)
        {
            contributed.resize(pImpl->mObservations.size());
        }
        for (const auto &feasibleArrivalSet : pImpl->mFeasibleArrivalSets)
        {
            if (saveContributingArrivals)
            {
                std::fill(contributed.begin(), contributed.end(), false);
            }
            double stackContribution{0};
            for (const auto iArrival : feasibleArrivalSet)
            {
                auto uniqueStationPhaseIndex
                    = pImpl->mArrivalToUniqueStationPhaseMap[iArrival];
                double estimate = uniqueTravelTimes[uniqueStationPhaseIndex];
                double observation = pImpl->mObservations[iArrival];
                double standardDeviation = pImpl->mStandardDeviations[iArrival];
                constexpr bool acceptAnyway = true;
                double likelihood = ::truncatedGaussianContribution(observation,
                                                       estimate,
                                                       standardDeviation,
                                                       normalization,
                                                       pImpl->mMaximumStandardDeviation,
                                                       acceptAnyway);
                if (saveContributingArrivals && likelihood > 0)
                {
                    contributed.at(iArrival) = true;
                }
                stackContribution = stackContribution + likelihood;
            }
            if (stackContribution > bestStackContribution)
            {
                if (saveContributingArrivals)
                {
                    bestSet = feasibleArrivalSet;
                    if (saveContributingArrivals){arrivalContributed = contributed;}
                }
                bestStackContribution = stackContribution;
            }
        }
        stack = 0;
        if (bestStackContribution >-1){stack = bestStackContribution;}
    }
    // If we are saving the contributing arrivals then note the travel times
    if (saveContributingArrivals)
    {
        pImpl->mContributingArrivals.clear();
        std::set<int> contributionList;
        if (getSignalType() == IMigrator::SignalType::DoubleDifference)
        {
            auto nPairs
                = static_cast<int> (pImpl->mDifferentialTablePairs.size());
            for (int iPair = 0; iPair < nPairs; ++iPair)
            {
                auto iTable = pImpl->mDifferentialTablePairs[iPair].first;
                auto jTable = pImpl->mDifferentialTablePairs[iPair].second;
                auto arrivalIndexPair
                    = pImpl->mDifferentialArrivalIndicesPairs[iPair];
                auto iArrival = arrivalIndexPair.first;
                auto jArrival = arrivalIndexPair.second;
                auto tEstimated1 = uniqueTravelTimes[iTable];
                auto tEstimated2 = uniqueTravelTimes[jTable];
                if (arrivalContributed[iArrival])
                {
                    if (!contributionList.contains(iArrival))
                    {
                        auto arrival = pImpl->mArrivals[iArrival];
                        arrival.setTravelTime(tEstimated1);
                        pImpl->mContributingArrivals.push_back(
                            std::pair {std::move(arrival),
                                       cumulativeContributions[iArrival]});
                        contributionList.insert(iArrival);
                    }
                }
                if (arrivalContributed[jArrival])
                {
                    if (!contributionList.contains(jArrival))
                    {
                        auto arrival = pImpl->mArrivals[jArrival];
                        arrival.setTravelTime(tEstimated2);
                        pImpl->mContributingArrivals.push_back(
                            std::pair {std::move(arrival),
                                       cumulativeContributions[jArrival]});
                        contributionList.insert(jArrival);
                    }
                }
            }
        }
        else
        {
            double normalization = ::truncatedGaussianNormalization(
                pImpl->mMaximumStandardDeviation);
            for (const auto iArrival : bestSet)
            {
                auto uniqueStationPhaseIndex
                    = pImpl->mArrivalToUniqueStationPhaseMap[iArrival];
                double estimate = uniqueTravelTimes[uniqueStationPhaseIndex];
                double observation = pImpl->mObservations[iArrival];
                double standardDeviation = pImpl->mStandardDeviations[iArrival];
                constexpr bool acceptAnyway = false;
                auto likelihood = ::truncatedGaussianContribution(observation,
                                                     estimate,
                                                     standardDeviation,
                                                     normalization,
                                                     pImpl->mMaximumStandardDeviation,
                                                     acceptAnyway);

                if (likelihood > 0)
                {
//std::cout << "keep it! " << likelihood << std::endl;
                    auto arrival = pImpl->mArrivals[iArrival];
                    arrival.setTravelTime(estimate - originTime);
                    pImpl->mContributingArrivals.push_back(
                        std::pair {std::move(arrival), likelihood});
                }
            }
        }
        // Sort by arrival time
        std::sort(pImpl->mContributingArrivals.begin(),
                  pImpl->mContributingArrivals.end(),
                  [](const std::pair<Arrival, double> &lhs,
                     const std::pair<Arrival, double> &rhs)
                  {
                     return lhs.first.getTime() < rhs.first.getTime();
                  });
        pImpl->mHaveContributingArrivals = true;
    }
    return stack;
}


std::vector<std::pair<Arrival, double>> 
    IMigrator::getContributingArrivals() const
{
    if (!haveArrivals()){throw std::runtime_error("Arrivals not set");}
    if (!pImpl->mHaveContributingArrivals)
    {
        throw std::runtime_error("Contributing arrival list not tabulated");
    }
    return pImpl->mContributingArrivals;
} 

/// Maximum source-receiver distance
void IMigrator::setMaximumEpicentralDistance(double distance)
{
    if (distance <= 0)
    {
        throw std::invalid_argument("Maximum distance must be positive");
    }
    pImpl->mMaximumEpicentralDistance = distance;
} 

double IMigrator::getMaximumEpicentralDistance() const noexcept
{
    return pImpl->mMaximumEpicentralDistance;
}

/// Travel time calculator exists?
bool IMigrator::haveTravelTimeCalculator(
    const std::string &network,
    const std::string &station,
    const std::string &phase) const noexcept
{
    if (!haveTravelTimeCalculatorMap()){return false;}
    try
    {
        ULocator::Station uStation;
        uStation.setNetwork(network);
        uStation.setName(station);
        return pImpl->mTravelTimeCalculatorMap->contains(uStation, phase);
    }
    catch (...)
    {
    }
    return false;
}

/// Destructor
IMigrator::~IMigrator() = default;
