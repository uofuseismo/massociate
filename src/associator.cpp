#include <string>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include "massociate/associator.hpp"
#include "massociate/associatorParameters.hpp"
#include "massociate/migrationParameters.hpp"
#include "massociate/event.hpp"
#include "massociate/pick.hpp"
#include "massociate/arrival.hpp"
#include "massociate/waveformIdentifier.hpp"
#include "massociate/migrate.hpp"
#include "massociate/mesh/spherical/points3d.hpp"
#include "massociate/mesh/cartesian/points3d.hpp"
#include "private/weightedStatistics.hpp"
#include "private/dbscan.hpp"

using namespace MAssociate;

namespace
{

struct Points3D
{
    double x = 0;
    double y = 0;
    double z = 0;
};

std::vector<MAssociate::Pick>
    getPicksInWindow(const std::vector<MAssociate::Pick> &picks,
                     const double T0, const double T1,
                     const bool removeT0 = true)
{
    std::vector<MAssociate::Pick> result;
    result.reserve(256);
    for (const auto &pick : picks)
    {
        auto t = pick.getTime();
        if (t >= T0 && t < T1)
        {
            result.push_back(pick);
        }
    }
    if (removeT0)
    {
        for (int i=0; i<static_cast<int> (result.size()); ++i)
        {
            result[i].setTime(result[i].getTime() - T0);
        }
    }
    return result;
}

/// @brief DBSCAN clusters on origin time however there is no requirement that
///        a station with two picks of the same phase cannot be clustered or
///        or a P arrival precede an S arrival for a given station.
///        This routine fixes that by doing the following:
///         (1) Ensures each cluster has picks that are causal.
///         (2) Prevents the same station from contributing the same phase
///             arrival more than once to a cluster.
///         (3) Reassigns cluster labels so that the cluster whose cumulative
///             contribution is labeled 0, then second highest is 1, etc.
///             so that it the `0' cluster is first. 
///         (4) Ensures each cluster has a minimum number of picks for
///             nucleation.
/// @param[in] nClusters       The number of clusters found by DBSCAN.
/// @param[in] minClusterSize  The minimum cluster size (size to nucleate).
/// @param[in] labels          The cluster labels created by DBSCAN.  Note that
///                            -1 indicates the pick is not assigned to a
///                            cluster.
/// @param[in] contributions   The contribution of each pick to the migration
///                            image.
/// @param[in] picks           The picks to causally cluster.
/// @param[out] newLabels      The new labels for each pick.  Here the 0'th
///                            group has the largest sum of contributions
///                            to the migration image max.  The 1'st group
///                            has the second largest contribution, etc.
///                            As before, -1 indicates the pick is not assigned
///                            to a group.
/// @param[out] nNewClusters   The number of new clusters.
void makeDBSCANClustersCausal(const int nClusters,
                              const int  minClusterSize,
                              const std::vector<int> &labels,
                              const std::vector<double> &contributions,
                              const std::vector<MAssociate::Pick> &picks,
                              std::vector<int> *newLabels,
                              int *nNewClusters)
{
    auto nPicks = static_cast<int> (picks.size());
    *nNewClusters = 0;
    newLabels->resize(nPicks, -1);
    std::vector<int> tempLabels(labels.size(), -1);
    if (nClusters == 0){return;} // Nothing to do
    std::vector<bool> lchecked(nPicks, false);
    std::vector<std::pair<int, double>> ranks(nClusters);
    // Fix each cluster:
    //   (1) Ensure the cluster is causal
    //   (2) Note the cluster's cumulative contribution sum 
    for (int ic=0; ic<nClusters; ++ic)
    {
        std::fill(lchecked.begin(), lchecked.end(), false);     
        for (int i=0; i<nPicks; ++i)
        {
            // If the pick is in the cluster and the label hasnt been created
            if (labels[i] == ic && !lchecked[i])
            {
                auto waveid1 = picks[i].getWaveformIdentifier();
                auto phase1 = picks[i].getPhaseName();
                int bestIndex = i;
                auto maxContribution = contributions[i];
                for (int j=i+1; j<nPicks; ++j)
                {
                    // Check this label is in the cluster
                    if (labels[j] == ic && !lchecked[j])
                    {
                        // If the SNLC/phase match then run the tie-breaker
                        auto waveid2 = picks[j].getWaveformIdentifier();
                        auto phase2 = picks[j].getPhaseName();
                        if (phase1 == phase2 && waveid1 == waveid2)
                        {
                            if (contributions[j] > maxContribution)
                            {
                                bestIndex = j;
                                maxContribution = contributions[j];
                            }
                            else
                            {
                                lchecked[j] = true;
                            }
                        } // Check check on phase and waveid match
                    } // Ehd check on another pick in this cluster
                } // Loop on `upper triangle' of pick matrix
                tempLabels[bestIndex] = ic; // Add pick that contributed most
                lchecked[i] = true;
            }
        } // Loop on picks 
        // Ensure that P arrivals come before S arrivals
        for (int i=0; i<nPicks; ++i)
        {
            if (tempLabels[i] == ic)
            {
                auto waveid1 = picks[i].getWaveformIdentifier();
                auto phase1 = picks[i].getPhaseName();
                auto time1 = picks[i].getTime();
                for (int j=i+1; j<nPicks; ++j)
                {
                    if (tempLabels[j] == ic)
                    {
                        auto waveid2 = picks[j].getWaveformIdentifier();
                        auto phase2 = picks[j].getPhaseName();
                        auto time2 = picks[j].getTime();
                        if (waveid1 == waveid2)
                        {
                            // Make sure P precedes S
                            if (phase1 == "P" && phase2 == "S")
                            {
                                if (time1 >= time2) // P arrives after S
                                {
std::cout << "p after s - removing tempLabels" << std::endl;
                                    if (contributions[i] > contributions[j])
                                    {
                                        tempLabels[j] =-1;
                                    }
                                    else
                                    {
                                        tempLabels[i] =-1;
                                    }
                                } 
                            }
                            else if (phase1 == "S" && phase2 == "P")
                            {
std::cout << "p after s - removing tempLabels 2" << std::endl;
                                if (time1 <= time2) // S arrives before P
                                {
                                    if (contributions[i] > contributions[j])
                                    {
                                        tempLabels[j] =-1;
                                    }
                                    else
                                    {
                                        tempLabels[i] =-1;
                                    }
                                }
                            }
                        } // End check on waveid match
                    } // End check on tempLabel j in this cluster
                } // Loop on j picks
            } // End check on templabel i in this cluster
        } // Loop on picks
        // Ensure this cluster has enough picks to nucleate
        auto nSubCluster = std::count(tempLabels.begin(), tempLabels.end(), ic);
        double contributionSum = 0;
        if (nSubCluster < minClusterSize)
        {
            nSubCluster = 0;
            contributionSum =-1;
        }
        else
        {
            contributionSum = 0;
            for (int i=0; i<nPicks; ++i)
            {
                if (tempLabels[i] == ic)
                {
                    contributionSum = contributionSum + contributions[i];
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
    //  (2) If a causal cluster is too small to nucleate then release all piocks
    *nNewClusters = 0;
    for (int ic=0; ic<nClusters; ++ic)
    {
        int newLabel = ranks[ic].first;
        if (ranks[ic].second < 0)
        {
            newLabel =-1;
        }
        else
        {
            *nNewClusters = *nNewClusters + 1;
        }
        for (int i=0; i<nPicks; ++i)
        {
            if (tempLabels[i] == ic){newLabels->at(i) = newLabel;}
        }
    }
}

/*
/// @brief Attempts to steal picks from other clusters.
///  
template<class T>
void stealPicks(std::vector<MAssociate::Event> *mEvents,
                const MAssociate::Migrate<T> &migrate)
{
    if (mEvents->empty()){return;}
    auto nEvents = static_cast<int> (mEvents->size());
    for (int iev=0; iev<nEvents; ++iev)
    {
        for (int jev=0; jev<nEvents; ++jev)
        {
        }
    }
}
*/

/*
int createCausalClusterFromContribution(
    const int cluster,
    const std::vector<int> &labels,
    const std::vector<double> &contributions,
    const std::vector<MAssociate::Pick> &picks,
    const int minClusterSize,
    std::vector<int> *newLabels,
    double *contributionSum)
{
    bool lProcess = true;
    *contributionSum = 0;
    // Initialize
    auto nPicks = static_cast<int> (labels.size());
    newLabels->resize(nPicks, -1);
    auto nSubCluster = std::count(labels.begin(), labels.end(), cluster);
    // Not enough picks to nucleate so don't bother
    if (nSubCluster < minClusterSize)
    {
        nSubCluster = 0;
        return nSubCluster;
    }
    // Investigate each pick 
    std::vector<bool> lchecked(nPicks, false);
    for (int i=0; i<nPicks; ++i)
    {
        // If the pick is in this cluster and the new label has not been creaetd
        if (labels[i] == cluster && !lchecked[i])
        {
            auto waveid1 = picks[i].getWaveformIdentifier();
            auto phase1 = picks[i].getPhaseName();
            int bestIndex = i;
            double maxContribution = contributions[i];
            for (int j=i+1; j<nPicks; ++j)
            {
                // Check this label is in the cluster
                if (labels[j] == cluster && !lchecked[j])
                {
                    // If the SNLC/phase match then run the tie-breaker
                    auto waveid2 = picks[j].getWaveformIdentifier();
                    auto phase2 = picks[j].getPhaseName(); 
                    if (phase1 == phase2 && waveid1 == waveid2)
                    {
                        if (contributions[j] > maxContribution) 
                        {
                            bestIndex = j;
                            maxContribution = contributions[j];
                        }
                        else
                        {
                            lchecked[j] = true;
                        }
                    }
                }
            } // End other phases
            newLabels->at(bestIndex) = 0; // Add pick that contributed most
            lchecked[i] = true;
        } // End check on if this pick should be scrutinized
    } // Loop on picks
    // If there are too few picks to nucleate then leave all picks unassociated
    nSubCluster = std::count(newLabels->begin(), newLabels->end(), 0);
    if (nSubCluster < minClusterSize)
    {
        lProcess = false;
        nSubCluster = 0;
        *contributionSum = 0;
        std::fill(newLabels->begin(), newLabels->end(), -1);
    }
    else
    {
        // Otherwise get the sum of the contributions
        double sum = 0; 
        for (int i=0; i<nPicks; ++i)
        {
            if (newLabels->at(i) == 0){sum = sum + contributions[i];}
        }
        *contributionSum = sum;
    }
    return nSubCluster;
} 
*/

}
///--------------------------------------------------------------------------///
///                               Implementation                             ///
///--------------------------------------------------------------------------///
template<class T>
class Associator<T>::AssociatorImpl
{
public:
    Points3D indexToPoints3D(const int index)
    {
        Points3D result;
        if (mGeometry == MAssociate::Geometry::SPHERICAL_POINTS_3D)
        {
            result.x = mSphericalPoints3D->getLongitude(index);
            result.y = mSphericalPoints3D->getLatitude(index);
            result.z = mSphericalPoints3D->getDepth(index);
        }
        else if (mGeometry == MAssociate::Geometry::CARTESIAN_POINTS_3D)
        {
            result.x = mCartesianPoints3D->getXPosition(index);
            result.y = mCartesianPoints3D->getYPosition(index);
            result.z = mCartesianPoints3D->getZPosition(index);
        }
        return result;
    }
    /// The nucleated events
    std::vector<MAssociate::Event> mEvents;
    /// The optimal indices in the grid for each event.
    std::vector<int> mOptimalIndices; 
    /// The assocation parameters
    MAssociate::AssociatorParameters mParameters;
    /// The class responsible for migrating differential pick times
    Migrate<T> mMigrate;
    /// The picks to associate
    std::vector<MAssociate::Pick> mPicks;
    /// Setting a lot of picks can be very time consuming because I check for
    /// repeat pick IDs.  This speeds up that activity.
    std::vector<uint64_t> mPickIDs;
    //std::vector<int> mAttemptedAssociations; // TODO is this used?
    /// The geometry
    std::unique_ptr<MAssociate::Mesh::Spherical::Points3D<T>>
        mSphericalPoints3D = nullptr;
    std::unique_ptr<MAssociate::Mesh::Cartesian::Points3D<T>>
        mCartesianPoints3D = nullptr;
    MAssociate::Geometry mGeometry = MAssociate::Geometry::UNKNOWN;
    /// Event ID counter
    uint64_t mEventID = 0;
    /// Max differential travel time in seconds that can be migrated    
    T mMaxDifferentialTime = 0;
    /// Flag indicating the class is initialized
    bool mInitialized = false;
};

/// C'tor
template<class T>
Associator<T>::Associator() :
    pImpl(std::make_unique<AssociatorImpl> ())
{
}

/// Destructor
template<class T>
Associator<T>::~Associator() = default;

/// Clears the class
template<class T>
void Associator<T>::clear() noexcept
{
    pImpl->mEvents.clear();
    pImpl->mOptimalIndices.clear();
    pImpl->mMigrate.clear();
    pImpl->mPicks.clear();
    pImpl->mPickIDs.clear();
    //pImpl->mAttemptedAssociations.clear();
    if (pImpl->mSphericalPoints3D)
    {
        pImpl->mSphericalPoints3D->clear();
    }
    if (pImpl->mCartesianPoints3D)
    {
        pImpl->mCartesianPoints3D = nullptr;
    }
    pImpl->mSphericalPoints3D = nullptr;
    pImpl->mCartesianPoints3D = nullptr;
    pImpl->mGeometry = MAssociate::Geometry::UNKNOWN;
    pImpl->mMaxDifferentialTime =-1;
    pImpl->mEventID = 0;
    pImpl->mInitialized = false;
} 

/// Initialize the class
template<class T>
void Associator<T>::initialize(const AssociatorParameters &parameters,
                               const Mesh::IMesh<T> &mesh)
{
    clear();
    // Check some of the parameters
    MAssociate::AssociatorParameters parmsWork(parameters);
    if (!parmsWork.haveNumberOfTravelTimeTables())
    {
        throw std::invalid_argument(
           "Number of travel time tables must be specified");
    }
    auto nTables = parmsWork.getNumberOfTravelTimeTables();
    if (nTables < 2)
    {
        throw std::invalid_argument(
           "There must be at least 2 travel time tables");
    }
    // Check the geometry
    int nPoints = 0;
    pImpl->mGeometry = mesh.getGeometry();
    if (pImpl->mGeometry == MAssociate::Geometry::SPHERICAL_POINTS_3D)
    {
        pImpl->mSphericalPoints3D = mesh.cloneSphericalPoints3D();
        if (!pImpl->mSphericalPoints3D->haveLatitudes() ||
            !pImpl->mSphericalPoints3D->haveLongitudes() ||
            !pImpl->mSphericalPoints3D->haveDepths())
        {
            throw std::invalid_argument(
               "spherical points mesh must have lats, lons, and depths");
        }
        nPoints = pImpl->mSphericalPoints3D->getNumberOfPoints();
    }
    else if (pImpl->mGeometry == MAssociate::Geometry::CARTESIAN_POINTS_3D)
    {
        pImpl->mCartesianPoints3D = mesh.cloneCartesianPoints3D();
        if (!pImpl->mCartesianPoints3D->haveXPositions() ||
            !pImpl->mCartesianPoints3D->haveYPositions() ||
            !pImpl->mCartesianPoints3D->haveZPositions())
        {
            throw std::invalid_argument(
               "cartesian points mesh must have x, y, and z");
        }
        nPoints = pImpl->mCartesianPoints3D->getNumberOfPoints();
    }
    else
    {
        throw std::runtime_error("Unhandled geometry");
    }
    // Set the migration parameters and initialize that engine
    if (nPoints < 1)
    {
        throw std::invalid_argument("No points in geometry");
    }
    // Create the migration parameters
    MAssociate::MigrationParameters migrationParms;
    migrationParms.setNumberOfPointsInTravelTimeTable(nPoints);
    migrationParms.setNumberOfTravelTimeTables(nTables);
    auto tileSize = nPoints;
    tileSize = std::min(tileSize, parmsWork.getTileSize());
    migrationParms.setTileSize(tileSize);
    migrationParms.setAnalyticCorrelationFunction(
        parmsWork.getAnalyticCorrelationFunction());
    pImpl->mMigrate.initialize(migrationParms);
    // Finish this off
    pImpl->mParameters = parmsWork;
    pImpl->mInitialized = true;
}

template<class T>
bool Associator<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Travel time tables
template<class T>
int Associator<T>::getNumberOfPointsInTravelTimeTable() const
{
    return pImpl->mMigrate.getNumberOfPointsInTravelTimeTable();
}

template<class T>
template<typename U>
void Associator<T>::setTravelTimeTable(const std::string &network,
                                       const std::string &station,
                                       const std::string &phase,
                                       int nPoints, const U times[])
{
    pImpl->mMaxDifferentialTime =-1;
    if (!isInitialized()){throw std::runtime_error("Class not inititialized");}
    pImpl->mMigrate.setTravelTimeTable(network, station, phase, nPoints, times);
    if (haveAllTravelTimeTables())
    {
        pImpl->mMaxDifferentialTime
           = getMaximumDifferentialTravelTime();
    }
}

template<class T>
bool Associator<T>::haveTravelTimeTable(
    const std::string &network, const std::string &station,
    const std::string &phase) const noexcept
{
    return pImpl->mMigrate.haveTravelTimeTable(network, station, phase);
}

template<class T>
bool Associator<T>::haveAllTravelTimeTables() const noexcept
{
    return pImpl->mMigrate.haveAllTravelTimeTables();
}
 
template<class T>
T Associator<T>::getMaximumDifferentialTravelTime() const noexcept
{
    if (pImpl->mMaxDifferentialTime < 0)
    {
        pImpl->mMaxDifferentialTime
            = pImpl->mMigrate.getMaximumDifferentialTravelTime();
    }
    return pImpl->mMaxDifferentialTime;
}

/*
/// Picks to bind
template<class T>
void Associator<T>::bindPickToEvent(const Pick &pick)
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (!pick.haveTime())
    {
        throw std::invalid_argument("time not set");
    }
    if (!pick.haveIdentifier())
    {
        throw std::invalid_argument("pick id not set");
    }
    // Nothing to do
    auto nEvents = getNumberOfEvents();
    if (nEvents < 1){return;}
    // If the migration engine doesn't have the table then quit early
    auto sncl = pick.getWaveformIdentifier();
    auto network = sncl.getNetwork();
    auto station = sncl.getStation();
    auto phase = pick.getPhaseName();
    auto pickID = pick.getIdentifier();
    if (!pImpl->mMigrate.haveTravelTimeTable(network, station, phase))
    {
        throw std::runtime_error("No travel time table for " + network
                               + "." + station + "." + phase);
    }
    // Find best fitting event
    MAssociate::Arrival arrival(pick);
    auto resMin = std::numeric_limits<double>::max();
    int eventIndex =-1;
    for (int iev=0; iev<nEvents; ++iev)
    {
        // Arrival already exists - don't add
        if (pImpl->mEvents[iev].haveArrival(pickID)){continue;}
        const bool overWriteIfExists = false;
        if (pImpl->mEvents[iev].canAddArrival(arrival, overWriteIfExists))
        {
            continue;
        }
        auto ot = pImpl->mEvents[iev].getOriginTime();
        auto idx = pImpl->mOptimalIndices[iev];
        auto tEst = ot
                  + pImpl->mMigrate.getTravelTime(network, station, phase,
                                                  idx);
         
    } 
}
*/

/// Picks
template<class T>
void Associator<T>::addPick(const Pick &newPick)
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (!pImpl->mMigrate.isInitialized())
    {
        throw std::runtime_error("migration engine not initialized");
    }
    if (!newPick.haveTime())
    {
        throw std::invalid_argument("time not set");
    }
    if (!newPick.haveIdentifier())
    {
        throw std::invalid_argument("pick id not set");
    }
    if (pImpl->mPicks.capacity() == 0)
    {
        pImpl->mPicks.reserve(8192);
        pImpl->mPickIDs.reserve(8192);
    }
    // If the migration doesn't have the table then quit early
    auto newSNCL = newPick.getWaveformIdentifier();
    auto networkNew = newSNCL.getNetwork();
    auto stationNew = newSNCL.getStation();
    auto phaseNew = newPick.getPhaseName();
    if (!pImpl->mMigrate.haveTravelTimeTable(networkNew,
                                             stationNew,
                                             phaseNew))
    {
        throw std::runtime_error("No travel time table for " + networkNew
                               + "." + stationNew + "." + phaseNew);
    }
    // Look for this pick
    auto pickID = newPick.getIdentifier();
    auto it = std::find(pImpl->mPickIDs.begin(), pImpl->mPickIDs.end(), pickID);
    if (it == pImpl->mPickIDs.end())
    {
        pImpl->mPicks.push_back(newPick);
        pImpl->mPickIDs.push_back(pickID);
    }
    else
    {
        std::cerr << "Pick already exists - overwriting" << std::endl;
        auto ia = std::distance(pImpl->mPickIDs.begin(), it);
        pImpl->mPicks[ia] = newPick;
    }
/*
    bool isNew = true;
    int ia = 0;
    for (const auto &pick : pImpl->mPicks)
    {
        if (pick.getIdentifier() == pickID)
        {
            std::cerr << "Pick already exists - overwriting" << std::endl;
            pImpl->mPicks[ia] = newPick;
            isNew = false;
            break;
        }
        ia = ia + 1;
    }
    if (isNew)
    {
        pImpl->mPicks.push_back(newPick);
        //pImpl->mPickIDs.push_back(pickID);
    }
*/
}

/// Clear the picks
template<class T>
void Associator<T>::clearPicks() noexcept
{
    pImpl->mPicks.clear();
    pImpl->mPickIDs.clear();
    
}

/// Gets the number of picks
template<class T>
int Associator<T>::getNumberOfPicks() const noexcept
{
    return static_cast<int> (pImpl->mPicks.size());
}

/// Associate
template<class T>
void Associator<T>::associate()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    auto nPicks = getNumberOfPicks();
    if (getNumberOfPicks() < 1)
    {
        std::cerr << "No picks" << std::endl;
        return;
    }
    // Sort the picks temporally
    std::sort(pImpl->mPicks.begin(), pImpl->mPicks.end(),
              [](const MAssociate::Pick &a, const MAssociate::Pick &b)
              {
                 return a.getTime() < b.getTime();
              });
    std::vector<bool> isAssociated(nPicks, false);
    // Break the problem up into temporal chunks.  Effectively, when picks 
    // are too far spaced in time they can't meaningfully interact through our
    // modeled differential travel times.
    auto picks = pImpl->mPicks;
    auto maxDT = getMaximumDifferentialTravelTime();
    auto minArrivalsToNucleate
        = pImpl->mParameters.getMinimumNumberOfArrivalsToNucleate();
    double rootT0 = picks[0].getTime();
    auto lastPickTime = picks.back().getTime();
    // Extract picks from window:
    //    [T0 - 3*maxDT : T0 + maxDT : T0 + 2.2*maxDT]
    auto T0 = rootT0;
    //int nClusters0 =-1; // Number of clusters in previous iteration
    while (true) //for (int kwin=0; kwin<nPicks; ++kwin)
    {
        //auto T1 = T0 + maxDT*2; // TODO make 2 a factor
        // There's a bit of a tradeoff here.  If we grab too many irrelavant
        // picks then our migration image will be corrupted and it will make
        // the greedy algorithm's job more difficult.  However, if the window
        // is too tight then we potentially miss relevant picks.
        auto maxOriginTime = T0 + maxDT; // Don't want events at end of window
        auto minPickTime = T0 - 2*maxDT;
        auto maxPickTime = T0 + 2*maxDT + 0.5; // 0.5 is like max noise on pick
        // Quitting time?
        if (maxOriginTime > lastPickTime + 2*maxDT){break;}
        // Get the picks in this chunk of time.
        auto localPicks = getPicksInWindow(picks, minPickTime, maxPickTime, true); 
        // Insufficient number of picks in window -> advance window
        if (localPicks.size() < minArrivalsToNucleate)
        {
            T0 = maxOriginTime;
            continue;
        }
        // Load the unassociated picks into the migration engine. 
        pImpl->mMigrate.clearPicks();
        for (const auto &localPick : localPicks)
        {
            pImpl->mMigrate.addPick(localPick);
        }
        // Migrate
        pImpl->mMigrate.migrate();
        if (!pImpl->mMigrate.haveImage())
        {
            T0 = maxOriginTime;
            continue;
        }
        // Attempt to cluster (origin times) given the maximum of the
        // migration image
        auto travelTimesToMaximum = pImpl->mMigrate.getTravelTimesToMaximum();
        auto contributions = pImpl->mMigrate.getContributionToMaximum();
        std::vector<double> originTimes(localPicks.size());
        for (int ip=0; ip<localPicks.size(); ++ip)
        {
            originTimes[ip] = localPicks[ip].getTime()
                            - travelTimesToMaximum[ip];
            //std::cout << ip << " " << originTimes[ip] << std::endl;
        }
        // Cluster origin times.  There is no causality here.  It is just an 
        // approximation of how the picks should be distributed.
        DBSCAN dbscan;
        dbscan.initialize(pImpl->mParameters.getDBSCANEpsilon(),
                          pImpl->mParameters.getDBSCANMinimumClusterSize()); 
        dbscan.setData(originTimes.size(), 1, originTimes.data());
        dbscan.cluster();
        auto nClusters = dbscan.getNumberOfClusters();
        auto labels = dbscan.getLabels();
        // Enforce causality in the clusters.  Tie-breaking between, say,
        // two P picks on CTU is done by comparing each contribution's to the
        // migration maximum.
        int nNewClusters;
        std::vector<int> newLabels;
        makeDBSCANClustersCausal(nClusters, minArrivalsToNucleate,
                                 labels, contributions, localPicks,
                                 &newLabels, &nNewClusters);
/*
for (int ic=0; ic<nNewClusters; ++ic)
{
for (int i=0; i<newLabels.size(); ++i)
{
 if (newLabels[i] == ic)
 {
     std::cout << "newlabels: " << ic << " " << picks[i].getIdentifier() << " " << picks[i].getWaveformIdentifier() << " " << picks[i].getPhaseName() << std::endl;
 } 
}
}
*/
        /// Next let's attempt to associate the largest cluster
        auto earliestOriginTime = std::numeric_limits<double>::max();
        for (int ic=0; ic<nNewClusters; ++ic)
        {
            std::vector<MAssociate::Pick> picksInCluster;
            pImpl->mMigrate.clearPicks();
            for (int ip=0; ip<static_cast<int> (labels.size()); ++ip)
            {
                if (newLabels[ip] == ic)
                {
                    pImpl->mMigrate.addPick(localPicks[ip]); 
                    picksInCluster.push_back(localPicks[ip]);
                }
            }
            // Locate this batch of picks and solve for an origin time
            if (pImpl->mMigrate.getNumberOfPicks() > minArrivalsToNucleate)
            {
                // Re-migrate
                pImpl->mMigrate.migrate();
                // Get the travel times and location.  Note, this has a
                // static correction.
                travelTimesToMaximum
                    =  pImpl->mMigrate.getTravelTimesToMaximum();
                // Compute origin time.  For least-squares this is the weighted
                // average.  For L1 this is the weighted median.
                std::vector<double> weights
                    = pImpl->mMigrate.getContributionToMaximum(true);
                std::vector<double> residuals(weights.size(), 0);
                for (int ip=0; ip<static_cast<int>(picksInCluster.size()); ++ip)
                {
                    // I define a residual as the observed - predicted time.
                    // Note that the static corrections have already been added
                    // to travelTimesToMaximum. 
                    residuals[ip] = picksInCluster[ip].getTime()
                                  - travelTimesToMaximum[ip]; 
                }
                double originTime = 0;
                if (pImpl->mParameters.getOriginTimeObjectiveFunction() ==
                    MAssociate::OriginTimeObjectiveFunction::L1)
                {
                    originTime = weightedMedian(residuals.size(),
                                                residuals.data(),
                                                weights.data());
                }
                else
                {
                    originTime = weightedMean(residuals.size(),
                                              residuals.data(),
                                              weights.data());
                }
                originTime = originTime + minPickTime; // Add in origin time
                // If the origin time is too close to the end of the window
                // then leave the picks as unassociated and continue
                earliestOriginTime = std::min(earliestOriginTime, originTime);

std::cout << std::setprecision(12);
                if (originTime > maxOriginTime)
                {
                    //std::cout << "skipping: " << maxOriginTime << std::endl;
                    continue;
                }
/*
            double originTime = 0;
            for (int ip=0; ip<static_cast<int> (picksInCluster.size()); ++ip)
            {
                    std::cout << ip << ","
                      << picksInCluster[ip].getWaveformIdentifier()
                      << "," << picksInCluster[ip].getIdentifier()
                      << "," << picksInCluster[ip].getPhaseName()
                      << "," << picksInCluster[ip].getTime() << std::endl;
                // The residual is the observed time - predicted time.
                // Note that the static corrections have already been added
                // to travelTimesToMaximum.
                double dt = picksInCluster[ip].getTime()
                          - travelTimesToMaximum[ip];
std::cout << "residual " << dt << std::endl;
                originTime = originTime + weights[ip]*dt;
            }
*/
//            originTime = originTime + T0; // Add in pick shift
                // Now create the event
                auto locationPair = pImpl->mMigrate.getImageMaximum();
                auto location = pImpl->indexToPoints3D(locationPair.first);
                MAssociate::Event event;
                event.setOriginTime(originTime);
                event.setXPosition(location.x);
                event.setYPosition(location.y);
                event.setZPosition(location.z);
#ifndef NDEBUG
std::cout << "Origin time: " << originTime << "(x,y,z)" << location.x << "," << location.y << "," << location.z << std::endl;
#endif
                for (int ip=0; ip<static_cast<int> (picksInCluster.size()); ++ip)
                {
#ifndef NDEBUG
                    std::cout << ip
                              << "," << picksInCluster[ip].getWaveformIdentifier()
                              << "," << picksInCluster[ip].getIdentifier()
                              << "," << picksInCluster[ip].getPhaseName()
                              << "," << picksInCluster[ip].getTime() << std::endl;
#endif
                    Arrival arrival(picksInCluster[ip]);
                    arrival.setTime(arrival.getTime() + minPickTime);
                    arrival.setTravelTime(travelTimesToMaximum[ip]
                                        - arrival.getStaticCorrection());
                    event.addArrival(arrival);
                }
                pImpl->mEventID = pImpl->mEventID + 1;
                event.setIdentifier(pImpl->mEventID);
                pImpl->mEvents.push_back(event);
                pImpl->mOptimalIndices.push_back(locationPair.first);
            }
            // Remove the picks in this cluster
            for (const auto &pickToRemove : picksInCluster)
            {
                int ip = 0;
                for (auto &p : picks)
                {
                    if (p.getIdentifier() == pickToRemove.getIdentifier())
                    {
                        picks.erase(picks.begin() + ip);
                        break;
                    }
                    ip = ip + 1;
                }
            }
        } // Loop
        // Update the time window
        if (nNewClusters == 0 || maxOriginTime <= earliestOriginTime)
        {
            T0 = maxOriginTime;
        }
/*
        // This loop will attempt to refine the above clusters in a 
        // causal fashion.  Effectively, this loop looks to assign
        // picks to earthquakes in the same location but with different
        // origin times.
        for (int kCluster=0; kCluster<nClusters; ++kCluster)
        {
            std::vector<std::pair<int, double>> ranks(nClusters);
            std::vector<std::vector<int>> newLabels(nClusters);
            // Enforce causality by tie-breaking with the size of the
            // contribution.  Then greedily strip out subclusters.
            for (int ic=0; ic<nClusters; ++ic)
            {
                double contributionSum = 0;
                auto nSubCluster = createCausalClusterFromContribution(
                                     ic, labels, contributions, localPicks,
                                     minArrivalsToNucleate,
                                     &newLabels[ic],
                                     &contributionSum);
                if (nSubCluster == 0){contributionSum =-1;}
                ranks[ic] = std::pair(ic, contributionSum);
            }
            // Sort in descending order based on cumulative contributions.
            std::sort(ranks.begin(), ranks.end(), 
                      [](const std::pair<int, double> &a,
                         const std::pair<int, double> &b)
                      {
                         return a.second > b.second;
                      });
            // Relocate all picks in largest cluster.
            int icMax = ranks[0].first;
            pImpl->mMigrate.clearPicks();
            std::vector<MAssociate::Pick> picksInCluster;
            picksInCluster.reserve(newLabels[icMax].size());
            for (int ip=0; ip<static_cast<int> (newLabels[icMax].size()); ++ip)
            {
                if (newLabels[icMax][ip] == 0)
                {
                    pImpl->mMigrate.addPick(localPicks[ip]);
                    picksInCluster.push_back(localPicks[ip]);
                }
            }
            // Whether or not this results in an event we have to remove
            // these picks.
            for (const auto p : picksInCluster)
            {
                for (int ip=0; ip<static_cast<int> (localPicks.size()); ++ip)
                {
                    if (p.getIdentifier() == localPicks[ip].getIdentifier())
                    {
                        localPicks.erase(localPicks.begin() + ip);
                        labels.erase(labels.begin() + ip);
                        break;
                    }
                }
            }

            // It's possible that the largest cluster is too small to migrate.
         pImpl->mMigrate.getNumberOfPicks(); 
            pImpl->mMigrate.migrate();
            // Get the travel times and location.  Note, this has a correction.
            travelTimesToMaximum = pImpl->mMigrate.getTravelTimesToMaximum();
            // Compute origin time.  For least-squares this is the weighted 
            // average.  Note, the weights are normalized such that they
            // sum to unity hence we don't divide in a subsequent step
            auto weights = pImpl->mMigrate.getContributionToMaximum(true);
            double originTime = 0;
            for (int ip=0; ip<static_cast<int> (picksInCluster.size()); ++ip)
            {
                    std::cout << ip << ","
                      << picksInCluster[ip].getWaveformIdentifier()
                      << "," << picksInCluster[ip].getIdentifier()
                      << "," << picksInCluster[ip].getPhaseName()
                      << "," << picksInCluster[ip].getTime() << std::endl;
                // The residual is the observed time - predicted time.
                // Note that the static corrections have already been added
                // to travelTimesToMaximum. 
                auto dt = picksInCluster[ip].getTime()
                        - travelTimesToMaximum[ip];
std::cout << "residual " << dt << std::endl;
                originTime = originTime + weights[ip]*dt;
            }
            originTime = originTime + T0; // Add in pick shift
            // Now create the event
            auto locationPair = pImpl->mMigrate.getImageMaximum();
            auto location = pImpl->indexToPoints3D(locationPair.first);
            MAssociate::Event event;
            event.setOriginTime(originTime); 
            event.setXPosition(location.x);
            event.setYPosition(location.y);
            event.setZPosition(location.z);
std::cout << std::setprecision(12);
std::cout << "Origin time: " << originTime << "(x,y,z)" << location.x << "," << location.y << "," << location.z << std::endl;
            for (int ip=0; ip<static_cast<int> (picksInCluster.size()); ++ip)
            {
                Arrival arrival(picksInCluster[ip]);
                arrival.setTravelTime(travelTimesToMaximum[ip]
                                    - arrival.getStaticCorrection());
                event.addArrival(arrival);
            }
            pImpl->mEventID = pImpl->mEventID + 1;
            event.setIdentifier(pImpl->mEventID);
            pImpl->mEvents.push_back(event);
            // And remove the associated picks
            for (const auto pickInCluster : picksInCluster)
            {
                for (int ip=0; ip<static_cast<int> (picks.size()); ++ip)
                {
                    if (picks[ip].getIdentifier() ==
                        pickInCluster.getIdentifier())
                    {
                        picks.erase(picks.begin() + ip);
                    }
                }
            }
        } // Loop on clusters
        // If there were no clusters in the window then iterate.  Otherwise,
        // attempt to process the window event.
        if (nClusters == 0)
        {
            T0 = T1;
        }
*/
    }
}

/// Gets the events
template<class T>
std::vector<MAssociate::Event> Associator<T>::getEvents() const
{
    return pImpl->mEvents;
}

/// Clears out the events
template<class T>
void Associator<T>::clearEvents(const bool resetEventCounter) noexcept
{
    pImpl->mEvents.clear();
    if (resetEventCounter){pImpl->mEventID = 0;}
}

/// Gets the number of events
template<class T>
int Associator<T>::getNumberOfEvents() const noexcept
{
    return static_cast<int> (pImpl->mEvents.size());
}

///--------------------------------------------------------------------------///
///                          Template Instantiation                          ///
///--------------------------------------------------------------------------///
template class MAssociate::Associator<double>;
template class MAssociate::Associator<float>;

template void MAssociate::Associator<double>::setTravelTimeTable(
    const std::string &network, const std::string &station,
    const std::string &phase, int nPoints, const double times[]);
template void MAssociate::Associator<double>::setTravelTimeTable(
    const std::string &network, const std::string &station,
    const std::string &phase, int nPoints, const float times[]);

template void MAssociate::Associator<float>::setTravelTimeTable(
    const std::string &network, const std::string &station,
    const std::string &phase, int nPoints, const double times[]);
template void MAssociate::Associator<float>::setTravelTimeTable(
    const std::string &network, const std::string &station,
    const std::string &phase, int nPoints, const float times[]);
