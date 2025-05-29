#include <iostream>
#include <iomanip>
#include <string>
#include <limits>
#include <cmath>
#include <map>
#include <set>
#ifdef WITH_LOCAL_UMPS_LOGGING
#include "logging/standardOut.hpp"
#else
#include <umps/logging/standardOut.hpp>
#endif
#include <uLocator/optimizers/originTime.hpp>
//#include <uLocator/station.hpp>
#include "massociate/associator.hpp"
#include "massociate/optimizer.hpp"
#include "massociate/migrator.hpp"
#include "massociate/arrival.hpp"
#include "massociate/clusterer.hpp"
#include "massociate/waveformIdentifier.hpp"
#include "massociate/event.hpp"
#include "massociate/pick.hpp"

using namespace MAssociate;

namespace
{

/*
auto pickComparitor = [](const Pick &lhs, const Pick &rhs)
{
    return lhs.getIdentifier() < rhs.getIdentifier();
};
*/

std::pair<int, std::vector<int>>
    makeCausalClusters(const int nClusters,
                       const int minimumClusterSize,
                       const std::vector<int> &labels,
                       const std::vector<std::pair<Arrival, double>> &contributions)
{
    auto nArrivals = static_cast<int> (contributions.size());
    std::vector<int> newLabels(nArrivals, -1);
    int nNewClusters{0};
    // Nothing to do
    if (nClusters == 0 || nArrivals == 0)
    {
        return std::pair{nNewClusters, newLabels};
    }   
    std::vector<bool> checked(nArrivals, false);
    std::vector<std::pair<int, double>> ranks(nClusters);
    std::vector<int> temporaryLabels(labels.size(), -1);
    // Fix each cluster:
    //   (1) Ensure the cluster is causal
    //   (2) Note the cluster's cumulative contribution sum 
    for (int ic = 0; ic < nClusters; ++ic)
    {
        std::fill(checked.begin(), checked.end(), false);
        for (int i = 0; i < nArrivals; ++i)
        {
            // If the pick is in the cluster and has not been checked...
            if (labels[i] == ic && !checked[i])
            {
                auto waveid1 = contributions[i].first.getWaveformIdentifier();
                auto phase1 = contributions[i].first.getPhase();
                int bestIndex = i; 
                auto maxContribution = contributions[i].second;
                for (int j = i + 1; j < nArrivals; ++j) 
                {
                    // Check this label is in the cluster
                    if (labels[j] == ic && !checked[j])
                    {
                        // If the SNLC/phase match then run the tie-breaker
                        auto waveid2
                            = contributions[j].first.getWaveformIdentifier();
                        auto phase2 = contributions[j].first.getPhase();
                        if (phase1 == phase2 && waveid1 == waveid2)
                        {
                            if (contributions[j].second > maxContribution)
                            {
                                bestIndex = j; 
                                maxContribution = contributions[j].second;
                            }
                            else
                            {
                                checked[j] = true;
                            }
                        } // Check check on phase and waveid match
                    } // End check on another pick in this cluster
                } // Loop on `upper triangle' of pick matrix
                // Save the pick that contributed most to the migration 
                temporaryLabels[bestIndex] = ic;
                checked[i] = true;
            }
        } // Loop on arrivals 
        // Ensure that P arrivals come before S arrivals
        for (int i = 0; i < nArrivals; ++i) 
        {
            if (temporaryLabels[i] == ic)
            {
                auto waveid1 = contributions[i].first.getWaveformIdentifier();
                auto phase1 = contributions[i].first.getPhase();
                auto time1 = contributions[i].first.getTime();
                for (int j = i + 1; j < nArrivals; ++j)
                {
                    if (temporaryLabels[j] == ic)
                    {
                        auto waveid2
                            = contributions[j].first.getWaveformIdentifier();
                        auto phase2 = contributions[j].first.getPhase();
                        auto time2 = contributions[j].first.getTime();
                        if (waveid1 == waveid2)
                        {
                            // Make sure P precedes S.  Otherwise, select the 
                            // pick that contributed more
                            if (phase1 == "P" && phase2 == "S")
                            {
                                if (time1 >= time2) // P arrives after S
                                {
//std::cout << "p after s - removing temporaryLabels" << std::endl;
                                    if (contributions[i].second >
                                        contributions[j].second)
                                    {
                                        temporaryLabels[j] =-1; // Disable S
                                    }
                                    else
                                    {
                                        temporaryLabels[i] =-1; // Disable P
                                    }
                                }
                            }
                            else if (phase1 == "S" && phase2 == "P")
                            {
//std::cout << "p after s - removing tempLabels 2" << std::endl;
                                if (time1 <= time2) // S arrives before P
                                {
                                    if (contributions[i].second >
                                        contributions[j].second)
                                    {
                                        temporaryLabels[j] =-1; // Disable P
                                    }
                                    else
                                    {
                                        temporaryLabels[i] =-1; // Disable S
                                    }
                                }
                            }
                        } // End check on waveid match
                    } // End check on temporary label j in this cluster
                } // Loop on j arrivals
            } // End check on temporary label i in this cluster
        } // Loop on arrivals
        // Ensure this cluster has enough arrivals to nucleate
        auto nSubCluster = std::count(temporaryLabels.begin(),
                                      temporaryLabels.end(), ic); 
        double contributionSum = 0;
        if (nSubCluster < minimumClusterSize)
        {
            nSubCluster = 0;
            contributionSum =-1;
        }
        else
        {
            contributionSum = 0;
            for (int i = 0; i < nArrivals; ++i)
            {
                if (temporaryLabels[i] == ic)
                {
                    contributionSum = contributionSum + contributions[i].second;
                }
           }
        }
        ranks[ic] = std::pair(ic, contributionSum);
    } // Loop on clusters
    // Now sort the clusters in descending order 
    std::sort(ranks.begin(), ranks.end(),
              [](const std::pair<int, double> &a,
                 const std::pair<int, double> &b)
              {
                 return a.second > b.second;
              });
    // Reprioritize the clusters so that:
    //  (1) The cluster who contributes the most is processed first 
    //  (2) If a causal cluster is too small to nucleate then release all
    //      arrivals
    nNewClusters = 0; 
    for (int ic = 0; ic < nClusters; ++ic)
    {    
        int newLabel = ranks[ic].first;
        if (ranks[ic].second < 0) 
        {
            newLabel =-1; 
        }
        else
        {
            nNewClusters = nNewClusters + 1; 
        }
        for (int i = 0;  i< nArrivals; ++i) 
        {
            if (temporaryLabels[i] == ic)
            {
                newLabels.at(i) = newLabel;
            }
        }
    }
    return std::pair {nNewClusters, newLabels};
}

/*
void purgePicks(const std::vector<std::pair<Arrival, double>> &clusterContributions,
                std::set<Pick, decltype(::pickComparitor)> &unassociatedPicks)
{
    for (const auto &arrivalsToRemove : clusterContributions)
    {
        auto identifier = arrivalsToRemove.first.getIdentifier();
        for (auto &pick : unassociatedPicks)
        {
            if (pick.getIdentifier() == identifier)
            {
                unassociatedPicks.erase(pick);
                break;
            }
        }
    }
}
*/

int countPArrivals(const std::vector<std::pair<Arrival, double>> &clusterContributions)
{
    int nP = 0;
    for (const auto &contribution : clusterContributions)
    {
        if (contribution.first.getPhase() == "P"){nP = nP + 1;}
    }
    return nP;
}

}

class Associator::AssociatorImpl
{
public:
    AssociatorImpl(std::shared_ptr<UMPS::Logging::ILog> logger = nullptr) :
        mLogger(logger)
    {
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> ();
        }
    }
    [[nodiscard]] double rankStation(const std::string &channel,
                                     const std::string &locationCode,
                                     const double uncertainty) const
    {
        int channelRank{100};
        if (channel.size() == 3)
        {
            auto channel2 = channel;
            channel2.pop_back();
            if (mChannelCodePreference.contains(channel2))
            {
                channelRank = mChannelCodePreference.find(channel2)->second;
            }
        }
        int locationCodeRank{0};
        if (mLocationCodePreference.contains(locationCode))
        {
            locationCodeRank
                = mLocationCodePreference.find(locationCode)->second;
        }
        return channelRank
             + locationCodeRank
             + std::min(10 - std::numeric_limits<double>::epsilon()*100, uncertainty);
    }
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    std::unique_ptr<IOptimizer> mOptimizer{nullptr};
    std::unique_ptr<IClusterer> mClusterer{nullptr};
    std::vector<Pick> mPicks;
    std::vector<Pick> mUnassociatedPicks;
    std::vector<Event> mEvents;
    std::map<std::string, int> mChannelCodePreference
    {
        std::pair {"GN", 1000}, // x
        std::pair {"DN", 1100}, // x
        std::pair {"HH", 1200}, // x
        std::pair {"BH", 1300}, // x
        std::pair {"EH", 1400}, // x
        std::pair {"HN", 1500}, // x
        std::pair {"EN", 1600} 
    };
    std::map<std::string, int> mLocationCodePreference
    {
        std::pair {"01", 10}, // UUSS typically uses 01
        std::pair {"02", 20}, 
        std::pair {"00", 30},
        std::pair {"--", 30},
        std::pair {"",   30}
    };
    double mPickCollocationTolerance{0.5};
    int mMinimumNumberOfArrivalsToNucleate{5};
    int mMinimumNumberOfPArrivalsToNucleate{3};
};

/// Constructor
Associator::Associator() :
    pImpl(std::make_unique<AssociatorImpl> ())
{
}

/// Constructor
Associator::Associator(std::shared_ptr<UMPS::Logging::ILog> &logger) :
    pImpl(std::make_unique<AssociatorImpl> (logger))
{
}

/// Destructor
Associator::~Associator() = default;

void Associator::initialize(int nArrivals, int nPArrivals)
{
    if (nPArrivals > nArrivals)
    {
        throw std::runtime_error(
        "P arrivals to nucleate exceeds number of arrivals to nucleate");
    }
    if (nArrivals < 3)
    {
        throw std::invalid_argument("Number of arrivals must be at least 3");
    }
    if (nPArrivals < 0)
    {
        throw std::invalid_argument("Number of P arrivals must be nonnegative");
    }
    pImpl->mMinimumNumberOfArrivalsToNucleate = nArrivals;
    pImpl->mMinimumNumberOfPArrivalsToNucleate = nPArrivals;
}

/// Number of arrivals required to nucleate an event
int Associator::getMinimumNumberOfArrivalsToNucleate() const noexcept
{
    return pImpl->mMinimumNumberOfArrivalsToNucleate;
}

/// Number of P arrivals required to nucleate an event
int Associator::getMinimumNumberOfPArrivalsToNucleate() const noexcept
{
    return pImpl->mMinimumNumberOfPArrivalsToNucleate;
}

/// Sets the optimizer
void Associator::setOptimizer(std::unique_ptr<IOptimizer> &&optimizer) 
{
    if (!optimizer->haveMigrator())
    {
        throw std::invalid_argument("Migrator not set on optimizer");
    }
    pImpl->mOptimizer = std::move(optimizer);
}

std::unique_ptr<IOptimizer> Associator::releaseOptimizer()
{
    if (!haveOptimizer()){std::runtime_error("Optimizer not set");}
    auto result = std::move(pImpl->mOptimizer);
    pImpl->mOptimizer = nullptr;
    return result;
}

bool Associator::haveOptimizer() const noexcept
{
    return pImpl->mOptimizer != nullptr;
}

/// Sets the clusterer
void Associator::setClusterer(std::unique_ptr<IClusterer> &&clusterer)
{
    if (clusterer == nullptr)
    {
        throw std::invalid_argument("Clusterer is null");
    }
/*
    if (!clusterer->haveMigrator())
    {   
        throw std::invalid_argument("Migrator not set on optimizer");
    }   
*/
    pImpl->mClusterer = std::move(clusterer);
}

std::unique_ptr<IClusterer> Associator::releaseClusterer()
{
    if (!haveClusterer()){std::runtime_error("Clusterer not set");}
    auto result = std::move(pImpl->mClusterer);
    pImpl->mClusterer = nullptr;
    return result;
}

bool Associator::haveClusterer() const noexcept
{
    return pImpl->mClusterer != nullptr;
}

bool Associator::haveTravelTimeCalculator(const Pick &pick) const
{
    if (!haveOptimizer()){return false;}
    if (!pick.haveWaveformIdentifier())
    {
        throw std::invalid_argument("Waveform identifier not set");
    }
    auto waveformIdentifier = pick.getWaveformIdentifier();
    auto network = waveformIdentifier.getNetwork();
    auto station = waveformIdentifier.getStation();
    std::string phaseHint{"P"}; // Default to P
    if (pick.havePhaseHint()){phaseHint = pick.getPhaseHint();}
    const auto migratorHandle = pImpl->mOptimizer->getMigratorHandle();
    return
        migratorHandle->haveTravelTimeCalculator(network, station, phaseHint);
}

void Associator::addTravelTimeCalculator(
    const ULocator::Station &station,
    const std::string &phase,
    std::unique_ptr<const ULocator::ITravelTimeCalculator> &&calculator)
{
    if (!haveOptimizer())
    {
        throw std::invalid_argument("Optimizer not set");
    }
    auto migratorHandle = pImpl->mOptimizer->getMigratorHandle(); 
    // Will throw
    migratorHandle->addTravelTimeCalculator(
        station, phase, std::move(calculator));
}

void Associator::setPicks(const std::vector<Pick> &picks)
{
    auto temporaryPicks = picks;
    setPicks(std::move(temporaryPicks));
}

void Associator::setPicks(std::vector<Pick> &&picks)
{
    if (!haveOptimizer()){throw std::runtime_error("Optimizer not set");}
    auto nInitialPicks = picks.size(); 
    pImpl->mPicks.clear();
    pImpl->mEvents.clear();
    pImpl->mUnassociatedPicks.clear();
    if (picks.empty()){return;}
    // First - validate the picks that can be added
    //const auto migratorHandle = pImpl->mOptimizer->getMigratorHandle();
    std::vector<Pick> pickList;
    pickList.reserve(picks.size());
    for (auto &pick : picks)
    {
        if (pick.haveWaveformIdentifier() && pick.haveTime())
        {
            auto waveformIdentifier = pick.getWaveformIdentifier();
            auto network = waveformIdentifier.getNetwork();
            auto station = waveformIdentifier.getStation();
            // These two things are required
            if (network.empty())
            {
                pImpl->mLogger->warn("Network code is empty; skipping"); 
                continue;
            }
            if (station.empty())
            {
                pImpl->mLogger->warn("Station code is empty; skpping");
                continue;
            }
            if (!pick.havePhaseHint())
            {
                auto waveformIdentifier = pick.getWaveformIdentifier();
                auto channel = waveformIdentifier.getChannel();
                if (channel.back() == 'Z' || channel.back() == 'P')
                {
                    pick.setPhaseHint(Pick::PhaseHint::P);
                }
                else if (channel.back() == 'N' ||
                         channel.back() == 'E' ||
                         channel.back() == '1' ||
                         channel.back() == '2' ||
                         channel.back() == 'S')
                {
                     pick.setPhaseHint(Pick::PhaseHint::S);
                }
            } 
            // We tried - at this stage just make it a P arrival
            if (!pick.havePhaseHint())
            {
                pick.setPhaseHint(Pick::PhaseHint::P);
            }
            try
            {
                if (!haveTravelTimeCalculator(pick))
                {
                    auto phase = pick.getPhaseHint();
                    pImpl->mLogger->warn("No calculator for "
                                       + network + "." + station + "." + phase
                                       + "; skipping...");
                    continue;
                 }
            }
            catch (const std::exception &e)
            {
                pImpl->mLogger->warn(e.what());
                continue;
            }
/*
            // Can we find this station in our map?
            try
            {
                auto phase = pick.getPhaseHint();
                if (!migratorHandle->haveTravelTimeCalculator(
                    network, station, phase))
                {
                    pImpl->mLogger->warn("No calculator for "
                                       + network + "." + station + "." + phase
                                       + "; skipping...");
                    continue;
                }
            }
            catch (const std::exception &e)
            {
                pImpl->mLogger->warn(e.what());
                continue;
            }
*/
            // Okay - add it
            pickList.push_back(std::move(pick));
        }
    }
    if (pickList.empty()){return;}
    // Next compute the station ranks (lower is better) and identifiers
    auto nCandidatePicks = static_cast<int> (pickList.size());
    std::vector<double> stationRanks(pickList.size(), 1.e30);
    for (int i = 0; i < nCandidatePicks; ++i)
    {
        auto waveformIdentifier = pickList[i].getWaveformIdentifier();
        auto uncertainty = pickList[i].getStandardError();
        auto channel = waveformIdentifier.getChannel();
        auto locationCode = waveformIdentifier.getLocationCode();
        stationRanks[i]
            = pImpl->rankStation(channel, locationCode, uncertainty);
    }
    // Now truncate duplicate picks on different channels
    std::vector<bool> keepPick(pickList.size(), false);
    for (int i = 0; i < nCandidatePicks; ++i)
    {
        auto waveformIdentifier1 = pickList[i].getWaveformIdentifier();
        auto pick1 = pickList[i].getTime().count()*1.e-6;
        auto phase1 = pickList[i].getPhaseHint();
        auto network1 = waveformIdentifier1.getNetwork();
        auto station1 = waveformIdentifier1.getStation();
        auto channel1 = waveformIdentifier1.getChannel();
        auto locationCode1 = waveformIdentifier1.getLocationCode();
        auto rank1 = stationRanks[i];
        bool haveBetterMatch{false};
        for (int j = 0; j < nCandidatePicks; ++j)
        {
            if (i == j){continue;}
            auto waveformIdentifier2 = pickList[j].getWaveformIdentifier();
            auto pick2 = pickList[j].getTime().count()*1.e-6;
            auto phase2 = pickList[j].getPhaseHint();
            auto network2 = waveformIdentifier2.getNetwork();
            auto station2 = waveformIdentifier2.getStation();
            auto channel2 = waveformIdentifier2.getChannel();
            auto locationCode2 = waveformIdentifier2.getLocationCode();
            // Potentially a duplicate
            if (network1 == network2 &&
                station1 == station2 &&
                phase1   == phase2 &&
                std::abs(pick2 - pick1) < pImpl->mPickCollocationTolerance)
            {
                // Same network/station/phase but different channel
                if (channel1 != channel2 || locationCode1 != locationCode2)
                {
                    // Lower score wins
                    if (stationRanks[j] < rank1)
                    {
                        haveBetterMatch = true;
                        break;
                    }
                }
            }
        } // Loop on other picks
        // A better match does not exist - let's associate this
        auto pickName = network1 + "." + station1 + "."
                      + channel1 + "." + locationCode1 + "."
                      + phase1;
        if (!haveBetterMatch) 
        {
            pImpl->mLogger->debug("Will retain " + pickName
                                + " with time " + std::to_string(pick1));
            keepPick[i] = true;
        }
        else
        {
            pImpl->mLogger->debug("Will not retain "
                                + pickName 
                                + " because a better pick exists");
        }
    }
    // Now let's keep the best of the matching picks
    pImpl->mPicks.reserve(pickList.size());
    uint64_t maxIdentifier{0};
    for (int i = 0; i < nCandidatePicks; ++i)
    {
        if (keepPick[i])
        {
            maxIdentifier = std::max(maxIdentifier,
                                     pickList[i].getIdentifier());
            pImpl->mPicks.push_back(std::move(pickList[i]));
        }
    }
    // Now clean up the pick identifiers.  Every time we come across
    // a duplicate identifier the duplicate identifier is set to some
    // gaurenteed larger value
    for (size_t i = 0; i < pImpl->mPicks.size(); ++i)
    {
        for (size_t j = i + 1; j < pImpl->mPicks.size(); ++j)
        {
            if (pImpl->mPicks[i].getIdentifier() ==
                pImpl->mPicks[j].getIdentifier())
            {
                maxIdentifier = maxIdentifier + 1; // Make sure this goes first
                pImpl->mPicks[j].setIdentifier(maxIdentifier);
            }
        }
    } 
    std::sort(pImpl->mPicks.begin(),
              pImpl->mPicks.end(),
              [](const Pick &lhs, const Pick &rhs)
              {
                 return lhs.getTime() < rhs.getTime();
              });
    pImpl->mUnassociatedPicks = pImpl->mPicks;
    pImpl->mLogger->debug("Set " + std::to_string(pImpl->mPicks.size())
                        + " out of " + std::to_string(nInitialPicks)
                        + " initial picks on associator"); 
}

void Associator::associate(const std::chrono::seconds &minimumOriginTime,
                           const std::chrono::seconds &maximumOriginTime)
{
    if (!haveClusterer()){throw std::runtime_error("Clusterer not set");}
    if (!haveOptimizer()){throw std::runtime_error("Optimizer not set");}
    // Clear the events and unassociated picks
    pImpl->mEvents.clear();
    pImpl->mUnassociatedPicks = pImpl->mPicks;
    //if (pImpl->mUnassociatedPicks.empty()){return;} // No picks
    // Insufficient number of picks
    if (static_cast<int> (pImpl->mUnassociatedPicks.size()) < 
        getMinimumNumberOfArrivalsToNucleate())
    {
        return;
    }
    // Create the set of unassociated picks
    auto pickComparitor = [](const Pick &lhs, const Pick &rhs)
    {
        return lhs.getIdentifier() < rhs.getIdentifier();
    };
    std::set<Pick, decltype(pickComparitor)> unassociatedPicks;
    for (const auto &pick : pImpl->mPicks)
    {
        if (!unassociatedPicks.contains(pick))
        {
            if (pick.getTime().count()*1.e-6 >= minimumOriginTime.count())
            {
                unassociatedPicks.insert(pick);
            }
        }
    }
    auto unassociatedPicksPreviousIteration = unassociatedPicks;

    ULocator::Optimizers::OriginTime originTimeCalculator;
    originTimeCalculator.setNorm(ULocator::Optimizers::IOptimizer::Norm::L1, 1);
    // Iterative loop with a bound.  We can only make as many clusters
    // as there are picks.
    auto nIterations = static_cast<int> (unassociatedPicks.size());
    for (int k = 0; k < nIterations; ++k)
    {
        // Did we make any progress last time?
        if (k > 0 &&
            unassociatedPicksPreviousIteration.size() ==
            unassociatedPicks.size())
        {
            break;
        }
        // Do we even have enough remaining picks to continue?
        if (static_cast<int> (unassociatedPicks.size()) < 
            getMinimumNumberOfArrivalsToNucleate())
        {
            pImpl->mLogger->debug("Insufficient number picks to associate");
            break;
        }
        auto nPArrivals
            = std::count_if(unassociatedPicks.begin(), unassociatedPicks.end(),
                            [](const Pick &p)
                            {
                                return p.getPhaseHint() == "P";
                            });
        if (nPArrivals < getMinimumNumberOfPArrivalsToNucleate())
        {
            pImpl->mLogger->debug(
               "Insufficient number of P picks to associate");
            break;
        }
        // Set the arrivals
        pImpl->mLogger->debug("Beginning iteration " + std::to_string(k)
                            + " with " + std::to_string(unassociatedPicks.size())
                            + " unassociated picks");
        std::vector<Arrival> arrivals;
        arrivals.reserve(unassociatedPicks.size());
        for (const auto &pick : unassociatedPicks)
        {
            arrivals.push_back(Arrival {pick});
        }
        try
        {
            pImpl->mOptimizer->setArrivals(arrivals);
        }
        catch (const std::exception &e)
        {
            pImpl->mLogger->error("Failed to set arrivals"
                                + std::string {e.what()});
            break;
        }
        // Optimize the migration image
        pImpl->mLogger->debug("Performing initial optimization for iteration: "
                            + std::to_string(k));
        try
        {
            pImpl->mOptimizer->optimize();
        }
        catch (const std::exception &e)
        {
            pImpl->mLogger->warn("Initial optimization failed with: "
                               + std::string {e.what()});
            break;
        }
        if (!pImpl->mOptimizer->haveOptimum())
        {
            pImpl->mLogger->debug("No optimum found; quitting outer loop");
            break;
        }
        // For the arrivals contributing to the maximum cluster on origin time
        auto contributions = pImpl->mOptimizer->getContributingArrivals();
        std::vector<double> originTimes, weights;
        originTimes.reserve(contributions.size());
        weights.reserve(contributions.size());
        for (const auto &contribution : contributions)
        {
            originTimes.push_back(contribution.first.getTime().count()*1.e-6
                                - contribution.first.getTravelTime());
            weights.push_back(contribution.first.getWeight());
        }
        // Shift
        double tShift = *std::min_element(originTimes.begin(),
                                          originTimes.end());
        std::transform(originTimes.begin(), originTimes.end(),
                       originTimes.begin(),
                       [=](const double t)
                       {
                           return t - tShift;
                       });
        try
        {
            pImpl->mClusterer->setData(originTimes.size(), 1, originTimes);
            pImpl->mClusterer->cluster();
        }
        catch (const std::exception &e)
        {
            pImpl->mLogger->debug("Clustering failed with: "
                                + std::string {e.what()});
            break;
        }
        auto nCandidateClusters = pImpl->mClusterer->getNumberOfClusters();
        auto labels = pImpl->mClusterer->getLabels();
//for (const auto &l : labels){std::cout << l << std::endl;}
//getchar();
        // Enforce causality in the clusters.  Tie-breaking between, say,
        // two P picks on CTU is done by comparing each contribution's to the
        // migration maximum.
        auto [nNewClusters, newLabels]
            = ::makeCausalClusters(nCandidateClusters,
                                   getMinimumNumberOfArrivalsToNucleate(),
                                   labels,
                                   contributions);
        if (nNewClusters == 0)
        {
            pImpl->mLogger->debug("No new clusters created - breaking...");
            break;
        }
//for (int i =0; i < newLabels.size(); ++i)
//{
// std::cout << labels[i] << " " << newLabels[i] << std::endl;
//}
#ifndef NDEBUG
        assert(newLabels.size() == labels.size());
#endif
        // For each cluster...
        // TODO: There's an optimization to be done.  Basically, what happens
        //       is if the new labels do not deviate from the original labels
        //       and there is only one cluster then we do not need to migrate
        //       again 
        for (int iCluster = 0; iCluster < nNewClusters; ++iCluster)
        {
            // Sufficient number of phase arrivals present to build the cluster?
            if (std::count(newLabels.begin(), newLabels.end(), iCluster) <
                getMinimumNumberOfArrivalsToNucleate())
            {
                pImpl->mLogger->debug("Will not refine cluster "
                                    + std::to_string(iCluster)
                                    + "; too few picks");
                continue;
            }
            auto nArrivalsInCluster = std::count(newLabels.begin(), newLabels.end(), iCluster);
            int nPArrivalsInCluster = ::countPArrivals(contributions);
            if (nArrivalsInCluster >= getMinimumNumberOfArrivalsToNucleate() &&
                nPArrivalsInCluster >= getMinimumNumberOfPArrivalsToNucleate())
            {
                // Extract the arrivals in this (causal) cluster 
                std::vector<MAssociate::Arrival> arrivalsInCluster;
                for (int i = 0; i < static_cast<int> (newLabels.size()); ++i)
                {
                    if (newLabels[i] == iCluster)
                    {
                        arrivalsInCluster.push_back(contributions.at(i).first);
                    }
                }
/*
            if (nPArrivalsInCluster < getMinimumNumberOfPArrivalsToNucleate())
            {
                pImpl->mLogger->debug("Will not refine cluster "
                                    + std::to_string(iCluster)
                                    + "; too few P picks");
                continue;
            }
*/
                try
                {
                    pImpl->mLogger->debug("Refining cluster " 
                                      + std::to_string(iCluster)
                                      + " with " 
                                      + std::to_string(arrivalsInCluster.size())
                                      + " arrivals");
                    pImpl->mOptimizer->setArrivals(arrivalsInCluster);
                    pImpl->mOptimizer->optimize();
                    if (!pImpl->mOptimizer->haveOptimum())
                    {
                        throw std::runtime_error("No optimum found");
                    }
                }
                catch (const std::exception &e)
                {
                    pImpl->mLogger->warn("Optimization for cluster "
                                       + std::to_string(iCluster)
                                       + " failed with: "
                                       + std::string {e.what()});
                    continue;
                }
                // One way or the other we are purging these arrivals
                auto clusterContributions
                     = pImpl->mOptimizer->getContributingArrivals();
                int nPArrivalsInEvent = countPArrivals(clusterContributions);
                for (const auto &arrivalsToRemove : clusterContributions)
                {
                    auto identifier = arrivalsToRemove.first.getIdentifier();
                    for (auto &pick : unassociatedPicks)
                    {
                        if (pick.getIdentifier() == identifier)
                        {
                            if (arrivalsToRemove.first.getPhase() == "P")
                            {
                                nPArrivalsInEvent = nPArrivalsInEvent + 1;
                            }
                            unassociatedPicks.erase(pick);
                            break;
                        }
                    }
                }
                // Does this event have too few arrivals to call it a proper event?
                auto eventType = Event::Type::Event;
                if (static_cast<int> (clusterContributions.size()) < 
                    getMinimumNumberOfArrivalsToNucleate() ||
                    nPArrivalsInEvent < getMinimumNumberOfPArrivalsToNucleate())
                {
                    eventType = Event::Type::Trigger;
                }
                // Now try to figure out the other event information
                std::vector<double> clusterArrivalTimes, clusterWeights;
                std::vector<double> clusterTravelTimes;
                clusterArrivalTimes.reserve(clusterContributions.size());
                clusterWeights.reserve(clusterContributions.size());
                clusterTravelTimes.reserve(clusterContributions.size());
                for (const auto &contribution : clusterContributions)
                {
                    clusterArrivalTimes.push_back(
                        contribution.first.getTime().count()*1.e-6);
                    clusterWeights.push_back(contribution.first.getWeight());
                    clusterTravelTimes.push_back(
                        contribution.first.getTravelTime());
                }
                double clusterOriginTime{0};
                if (pImpl->mOptimizer->getMigratorHandle()->getSignalType() ==
                    IMigrator::SignalType::DoubleDifference)
                {
                    try
                    {
                        originTimeCalculator.setArrivalTimes(clusterArrivalTimes,
                                                             clusterWeights);
                        originTimeCalculator.setTravelTimes(clusterTravelTimes);
                        originTimeCalculator.optimize();
                        clusterOriginTime = originTimeCalculator.getTime();
                    }
                    catch (const std::exception &e)
                    {
                        pImpl->mLogger->warn(
                            "Origin time calculation failed for cluster "
                          + std::to_string(iCluster));
                        continue;
                    }
                }
                else
                {
                    clusterOriginTime
                        = pImpl->mOptimizer->getOptimalOriginTime();
//std::cout << "lifted: " << clusterOriginTime  << std::endl;
                }
                if (clusterOriginTime < minimumOriginTime.count() ||
                    clusterOriginTime > maximumOriginTime.count())
                {
                    continue;
                }
                // Build the event
                auto [bestLatitude, bestLongitude, bestDepth]
                     = pImpl->mOptimizer->getOptimalHypocenter(); 
                Event event;
                try
                {
                    event.setLatitude(bestLatitude);
                    event.setLongitude(bestLongitude);
                    event.setDepth(bestDepth);
                    event.setOriginTime(clusterOriginTime);
                    event.setType(eventType);
                }
                catch (const std::exception &e)
                {
                    pImpl->mLogger->warn(
                         "Failed to set hypocentral info for cluster "
                        + std::to_string(iCluster) + ".  Failed with: "
                        + e.what());
                    continue;
                }
                // Sometimes, we can double bind a station/phase. Clean that up.
                std::vector<bool> stationPhaseMask(clusterContributions.size(), false);
                for (size_t iArrival = 0;
                     iArrival < clusterContributions.size(); ++iArrival)
                {
                    // Keep the best fitting
                    auto waveformIdentifier_i
                        = clusterContributions[iArrival].first.getWaveformIdentifier();
                    auto network_i = waveformIdentifier_i.getNetwork();
                    auto station_i = waveformIdentifier_i.getStation();
                    auto phase_i = clusterContributions[iArrival].first.getPhase();
                    double minimumAbsoluteResidual
                        = std::abs(clusterContributions[iArrival].first.getTime().count()*1.e-6
                                - (clusterOriginTime
                                 + clusterContributions[iArrival].first.getTravelTime()));
                    for (size_t jArrival = 0;
                         jArrival < clusterContributions.size(); ++jArrival)
                    {
                        if (iArrival == jArrival){continue;}
                        auto waveformIdentifier_j
                            = clusterContributions[jArrival].first.getWaveformIdentifier();
                        auto network_j = waveformIdentifier_j.getNetwork();
                        auto station_j = waveformIdentifier_j.getStation();
                        auto phase_j = clusterContributions[jArrival].first.getPhase();
                        // Match
                        if (network_i == network_j &&
                            station_i == station_j &&
                            phase_i == phase_j)
                        {
                            double residualJ
                                = std::abs(clusterContributions[jArrival].first.getTime().count()*1.e-6
                                        - (clusterOriginTime
                                         + clusterContributions[jArrival].first.getTravelTime()));
                            // Retain the best - mask the other
                            if (residualJ < minimumAbsoluteResidual)
                            {
                                minimumAbsoluteResidual = residualJ;
                                stationPhaseMask[iArrival] = true; // Mask i
                            }
                            else
                            {
                                stationPhaseMask[jArrival] = true; // Mask j
                            }
                        }
                    }
                }
                // Add the arrivals
                for (size_t iArrival = 0; iArrival < clusterContributions.size(); ++iArrival)
                {
                    if (!stationPhaseMask[iArrival])
                    {
                        try
                        {
                            event.addArrival(std::move(clusterContributions[iArrival].first));
                        }
                        catch (const std::exception &e)
                        {
                            pImpl->mLogger->warn(
                                "Could not add arrival to event.  Failed with "
                              + std::string {e.what()});
                        }
                    }
                    else
                    {
                        auto waveformIdentifier = clusterContributions[iArrival].first.getWaveformIdentifier(); 
                        auto phase = clusterContributions[iArrival].first.getPhase();
                        pImpl->mLogger->warn("Masked dupcliate arrival: "
                                           + waveformIdentifier.getNetwork() + "."
                                           + waveformIdentifier.getStation() + "."
                                           + phase);
                    }
                }
                if (event.getNumberOfArrivals() > 0)
                {
                    pImpl->mEvents.push_back(std::move(event));
                }
            }
        }
        unassociatedPicksPreviousIteration = unassociatedPicks;
        if (nNewClusters == 0)
        {
            pImpl->mLogger->debug("No new clusters created - breaking...");
            break;
        }
    }
    // Next, let's attempt to merge split events
    if (pImpl->mEvents.size() > 1)
    {
        auto nInitialEvents = static_cast<int> (pImpl->mEvents.size());
        std::vector<bool> deleteEvent(nInitialEvents, false);
        std::vector<Event> newEvents;
        bool truncateEvents{false};
        for (int i = 0; i < nInitialEvents; ++i)
        {
            if (deleteEvent[i]){continue;}
            auto originTime_i = pImpl->mEvents[i].getOriginTime().count();
            auto latitude_i   = pImpl->mEvents[i].getLatitude();
            auto longitude_i  = pImpl->mEvents[i].getLongitude(); 
            for (int j = i + 1; j < nInitialEvents; ++j)
            {
                if (deleteEvent[j]){continue;}
                auto originTime_j = pImpl->mEvents[j].getOriginTime().count();
                auto latitude_j   = pImpl->mEvents[j].getLatitude();
                auto longitude_j  = pImpl->mEvents[j].getLongitude();
                if (std::abs(originTime_i - originTime_j) < std::chrono::microseconds {1000000}.count() &&
                    std::abs(latitude_i - latitude_j) < 0.1 &&
                    std::abs(longitude_i - longitude_j) < 0.1)
                {
                    auto arrivals_i = pImpl->mEvents[i].getArrivals();
                    auto arrivals_j = pImpl->mEvents[j].getArrivals();
                    if (arrivals_i.size() >= arrivals_j.size())
                    { 
                        truncateEvents = true;
                        deleteEvent[j] = true; 
                        pImpl->mLogger->debug("Merging " + std::to_string(i) + " into " + std::to_string(j)); 
                        for (auto &arrival : arrivals_j)
                        {
                            if (pImpl->mEvents[i].canAddArrival(arrival, false) >= 0)
                            {
                                pImpl->mEvents[i].addArrival(arrival);
                            }
                        }
                    } 
                    else
                    {
                        truncateEvents = true;
                        deleteEvent[i] = true;
                        pImpl->mLogger->debug("Merging " + std::to_string(j) + " into " + std::to_string(i));
                        //auto event = pImpl->mEvents[j];
                        for (auto &arrival : arrivals_i)
                        {
                            if (pImpl->mEvents[j].canAddArrival(arrival, false) >= 0)
                            {
                                pImpl->mEvents[j].addArrival(arrival);
                            }
                        }
                    }
                }
            }
        }
        if (truncateEvents)
        {
            auto eventsCopy = pImpl->mEvents;
            pImpl->mEvents.clear();
            for (int i = 0; i < nInitialEvents; ++i)
            {
                if (!deleteEvent[i]){pImpl->mEvents.push_back(std::move(eventsCopy[i]));}
            }
            pImpl->mLogger->debug("Merged " + std::to_string(nInitialEvents -  pImpl->mEvents.size()) + " events");
        } 
    }
    // Finally build the definitive list of unassociated picks (i.e.,
    // picks that were not attached to events)
    for (const auto &event : pImpl->mEvents)
    {
        const auto &arrivalsInEvent = event.getArrivalsReference();
        for (const auto &arrival : arrivalsInEvent)
        {
            auto identifier = arrival.getIdentifier();
            pImpl->mUnassociatedPicks.erase(
                std::remove_if(pImpl->mUnassociatedPicks.begin(),
                               pImpl->mUnassociatedPicks.end(),
                               [=](const Pick &pick)
                               {
                                   return pick.getIdentifier() == identifier;
                               })); 
        }
    }
    for (const auto &pick : pImpl->mPicks)
    {
        if (pick.getTime().count()*1.e-6 < minimumOriginTime.count())
        {
            pImpl->mUnassociatedPicks.push_back(pick);
        }
    }
    auto nAssociated = pImpl->mPicks.size() - pImpl->mUnassociatedPicks.size();
    if (!pImpl->mEvents.empty())
    {
        std::sort(pImpl->mEvents.begin(), pImpl->mEvents.end(), 
                  [](const Event &lhs, const Event &rhs)
                  {
                     return lhs.getOriginTime() < rhs.getOriginTime();
                  });
        auto nTriggers = std::count_if(pImpl->mEvents.begin(),
                                       pImpl->mEvents.end(),
                                       [](const auto &event)
                                       {
                                          return event.getType() ==
                                                 Event::Type::Trigger;
                                       });
        auto nEarthquakes = static_cast<int> (pImpl->mEvents.size())
                          - nTriggers;
        pImpl->mLogger->info("Associated "
                           + std::to_string(nAssociated)
                           + " picks into "
                           + std::to_string(nEarthquakes)
                           + " event(s) and "
                           + std::to_string(nTriggers)
                           + " trigger(s)");
    }
}

/// The unassociated picks
std::vector<Pick> Associator::getUnassociatedPicks() const
{
    return pImpl->mUnassociatedPicks;
}

/// The events
std::vector<Event> Associator::getEvents() const
{
    return pImpl->mEvents;
}

/// The number of events
int Associator::getNumberOfEvents() const noexcept
{
    return static_cast<int> (pImpl->mEvents.size());
}
