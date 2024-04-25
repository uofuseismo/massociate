#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <uLocator/station.hpp>
#include <uLocator/travelTimeCalculatorMap.hpp>
#include <uLocator/travelTimeCalculator.hpp>
#include <uLocator/position/knownLocalLocation.hpp>
#include <uLocator/position/geographicRegion.hpp>
#include <umps/logging/standardOut.hpp>
#include "massociate/migrator.hpp"
#include "massociate/arrival.hpp"
#include "massociate/waveformIdentifier.hpp"
#include "analyticSignals.hpp"

using namespace MAssociate; 

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
    std::vector<std::pair<double, double>> mDifferentialArrivalPairs;
    std::vector<std::pair<double, double>>
         mDifferentialArrivalPairsStandardDeviations;
    std::vector<std::pair<int, int>> mDifferentialArrivalIndicesPairs;
    std::vector<std::pair<int, int>> mDifferentialTablePairs;
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
    std::unique_ptr<ULocator::TravelTimeCalculatorMap> &&map)
{
    if (map == nullptr)
    {   
        throw std::invalid_argument("Travel time calculator map is NULL");
    }   
    pImpl->mTravelTimeCalculatorMap = std::move(map);
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
    pImpl->mDifferentialArrivalPairs.clear();
    pImpl->mDifferentialArrivalPairsStandardDeviations.clear();
    pImpl->mDifferentialArrivalIndicesPairs.clear();
    pImpl->mDifferentialTablePairs.clear();
    pImpl->mContributingArrivals.clear();
    pImpl->mHaveContributingArrivals = false;

    pImpl->mArrivals.reserve(arrivals.size());
    pImpl->mUniqueStationPhases.reserve(arrivals.size());
    pImpl->mObservations.reserve(arrivals.size());
    pImpl->mStandardDeviations.reserve(arrivals.size());
    pImpl->mArrivalToUniqueStationPhaseMap.reserve(arrivals.size());

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
        pImpl->mObservations.push_back(arrival.getTime().count()*1.e-6);
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
            for (int j = i + 1; j < nArrivalsSet; ++j)
            {
                auto jUniqueStationPhase
                    = pImpl->mArrivalToUniqueStationPhaseMap.at(j);
                auto name2 = pImpl->mUniqueStationPhases.at(jUniqueStationPhase).first.getNetwork()
                           + "."
                           + pImpl->mUniqueStationPhases.at(jUniqueStationPhase).first.getName();
                auto phase2
                    = pImpl->mUniqueStationPhases.at(jUniqueStationPhase).second;
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
                    std1Use = std::sqrt(3)*std1;
                    std2Use = std::sqrt(3)*std2;
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
        throw std::runtime_error("Option not handled");
        pImpl->mLogger->debug("Set " + std::to_string(nArrivalsSet)
                            + " phase arrivals");
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
std::vector<double> IMigrator::evaluateAtKnownLocations() const
{
    auto nLocations = static_cast<int> (pImpl->mKnownLocations.size());
    std::vector<double> values(nLocations, 0.0);
    for (int i = 0; i < nLocations; ++i)
    {
        if (pImpl->mKnownLocations[i])
        {
            try
            {
                auto x = pImpl->mKnownLocations[i]->x();
                auto y = pImpl->mKnownLocations[i]->y();
                auto z = pImpl->mKnownLocations[i]->z();
                values[i] = evaluate(x, y, z);
            }
            catch (const std::exception &e)
            {
                pImpl->mLogger->warn("Failed to evaluate migration for "
                                    +       std::to_string(i)
                                    + ".  Failed with: "
                                    + std::string {e.what()});
               values[i] = 0; // Contribute nothing
            }
        }
    }
    return values;
}

/// Evaluate at a given point and depth
double IMigrator::evaluate(
    const double x, const double y, const double z) const
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
    constexpr double originTime{0};
    constexpr bool applyCorrection{true};
    pImpl->mTravelTimeCalculatorMap->evaluate(pImpl->mUniqueStationPhases,
                                              originTime,
                                              x, y, z,
                                              &uniqueTravelTimes,
                                              applyCorrection);
    auto uniqueDistances
        = pImpl->mTravelTimeCalculatorMap->computeDistance(
             pImpl->mUniqueStationPhases, x, y);
    double stack{0};
    if (getSignalType() == IMigrator::SignalType::DoubleDifference)
    {
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
        throw std::runtime_error("Absolute not yet handled");
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
