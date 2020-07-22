#include <string>
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

std::vector<int> createCausalClusterFromContribution(
    const int cluster,
    const std::vector<int> &labels,
    const std::vector<double> &contributions,
    const std::vector<MAssociate::Pick> &picks,
    const int minClusterSize = 4)
{
    // Do I have enough?
    auto nPicks = static_cast<int> (labels.size());
    std::vector<int> newLabels(nPicks, -1);
    auto nInCluster = std::count(labels.begin(), labels.end(), cluster);
    if (nInCluster < minClusterSize)
    {
        return newLabels;
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
            newLabels[bestIndex] = 0; // Add pick that contributed most
            lchecked[i] = true;
        } // End check on if this pick should be scrutinized
    } // Loop on picks
    // If there are too few picks to nucleate then leave all picks unassociated
    nInCluster = std::count(newLabels.begin(), newLabels.end(), 0);
    if (nInCluster < minClusterSize)
    {
        std::fill(newLabels.begin(), newLabels.end(), -1);
    }
    return newLabels;
} 

}

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
    std::vector<MAssociate::Event> mEvents;
    MAssociate::AssociatorParameters mParameters;
    Migrate<T> mMigrate;
    std::vector<MAssociate::Pick> mPicks;
    std::vector<int> mAttemptedAssociations;
    //std::unique_ptr<MAssociate::Mesh::IMesh<T>> mMesh = nullptr;
    std::unique_ptr<MAssociate::Mesh::Spherical::Points3D<T>>
        mSphericalPoints3D = nullptr;
    std::unique_ptr<MAssociate::Mesh::Cartesian::Points3D<T>>
        mCartesianPoints3D = nullptr;
    MAssociate::Geometry mGeometry = MAssociate::Geometry::UNKNOWN;
    uint64_t mEventID = 0;
    T mMaxDifferentialTime = 0;
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
    pImpl->mMigrate.clear();
    pImpl->mPicks.clear();
    pImpl->mAttemptedAssociations.clear();
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
    // If the table doesn't have a table then quit early
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
    bool isNew = true;
    auto pickID = newPick.getIdentifier();
    int ia = 0;
    for (const auto &pick : pImpl->mPicks)
    {
        if (pick.getIdentifier() == pickID)
        {
            std::cerr << "Pick already exists - overwriting" << std::endl;
            pImpl->mPicks[ia] = newPick;
            isNew = false;
        }
        ia = ia + 1;
    }
    if (isNew)
    {
        pImpl->mPicks.push_back(newPick);
        pImpl->mAttemptedAssociations.push_back(0);
    }
}

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
    auto rootT0 = picks[0].getTime();
    auto T0 = rootT0;
    int nClusters0 =-1; // Number of clusters in previous iteration
    for (int kwin=0; kwin<nPicks; ++kwin)
    {
        auto T1 = T0 + maxDT*2; // TODO make 2 a factor
        // Get the picks in this chunk of time.
        auto localPicks = getPicksInWindow(picks, T0, T1); 
        if (localPicks.size() < 4) // TODO fix this
        { 
            T0 = T1;
            continue;
        }
        // Load the unassociated picks into the migration engine. 
        pImpl->mMigrate.clearPicks();
        for (const auto &localPick : localPicks)
        {
            pImpl->mMigrate.addPick(localPick);
        }
int minPicksToNucleate = 6;
        pImpl->mMigrate.migrate();
        auto travelTimesToMaximum = pImpl->mMigrate.getTravelTimesToMaximum();
        auto contributions = pImpl->mMigrate.getContributionToMaximum();
        std::vector<double> originTimes(localPicks.size());
        for (int ip=0; ip<localPicks.size(); ++ip)
        {
            originTimes[ip] = localPicks[ip].getTime()
                            - travelTimesToMaximum[ip];  
        } 
        // Cluster origin times.  There is no causality here.  It is just an 
        // approximation of how the picks should be distributed.
        DBSCAN dbscan;
        dbscan.initialize(0.2, 5); 
        dbscan.setData(originTimes.size(), 1, originTimes.data());
        dbscan.cluster();
        auto nClusters = dbscan.getNumberOfClusters();
        auto labels = dbscan.getLabels();
        // This loop will attempt to refine the above clusters in a 
        // causal fashion.  Effectively, this loop looks to assign
        // picks to earthquakes in the same location but have different
        // origin times.
        for (int kCluster=0; kCluster<nClusters; ++kCluster)
        {
            std::vector<std::pair<int, double>> ranks(nClusters);
            std::vector<std::vector<int>> newLabels(nClusters);
            // Enforce causality by tie-breaking with the size of the
            // contribution.  Then greedily strip out subclusters.
            for (int ic=0; ic<nClusters; ++ic)
            {
                newLabels[ic] = createCausalClusterFromContribution(
                                     ic, labels, contributions, localPicks,
                                     minPicksToNucleate);
                int nSubCluster = 0;
                double sum = 0;
                for (int i=0; i<nPicks; ++i)
                {
                    if (newLabels[ic][i] == 0)
                    {
                        sum = sum + contributions[i];
                        nSubCluster = nSubCluster + 1;
                    }
                }
                if (nSubCluster == 0){sum =-1;}
                ranks[ic] = std::pair(ic, sum);
            }
            // Sort in descending order
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
                originTime = originTime + weights[ip]*dt;
            }
            originTime = originTime + T0; // Add in pick shift
std::cout << "Origin time: " << originTime << std::endl;
            // Now create the event
            auto locationPair = pImpl->mMigrate.getImageMaximum();
            auto location = pImpl->indexToPoints3D(locationPair.first);
            MAssociate::Event event;
            event.setOriginTime(originTime); 
            event.setXPosition(location.x);
            event.setYPosition(location.y);
            event.setZPosition(location.z);
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
        }
        // If there were no clusters in the window then iterate.  Otherwise,
        // attempt to process the window event.
        if (nClusters == 0)
        {
            T0 = T1;
        }
    }
/*
    // Sort the picks temporally

    // Perform initial migration
    pImpl->mMigrate.migrate();
    // Now cluster
*/
}

/// Gets the events
template<class T>
std::vector<MAssociate::Event> Associator<T>::getEvents() const
{
    return pImpl->mEvents;
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
