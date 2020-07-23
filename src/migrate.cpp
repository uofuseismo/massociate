#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <boost/align/aligned_allocator.hpp>
#ifdef USE_TBB
#include <tbb/tbb.h>
#endif
#include "massociate/migrate.hpp"
#include "massociate/migrationParameters.hpp"
#include "massociate/pick.hpp"
#include "massociate/waveformIdentifier.hpp"
#include "private/analyticSignals.hpp"
#include "private/dbscan.hpp"

#define ALIGNMENT 64
using namespace MAssociate;

namespace
{
/// @brief Computes the matrix padding length so that we get bit aligned memory.
/// @param[in] n              The array length.  This must be positive.
/// @param[in] precisionSize  The size of the precision, ex: sizeof(float) = 4.
/// @param[in] alignment      The byte alignment.  This should be a power of 2.
/// @result The padded array length.  This will be greater than or equal to n.
int padLength(const int n,
              const size_t precisionSize = sizeof(double),
              const int alignment = ALIGNMENT)
{
    auto size = static_cast<int> (precisionSize);
    int padLength = 0;
    auto xmod = (n*size)%alignment;
    if (xmod != 0){padLength = (alignment - xmod)/size;}
    auto nPointsPadded = n + padLength;
    return nPointsPadded;
}

/*!
 * @brief Defines the table name, e.g., NN.SS.PP
 * @param network   The network name.
 * @param station   The station name.
 * @param phase     The phase name.
 * @return The travel time table name.
 */
std::string makeTableName(const std::string &network,
                          const std::string &station,
                          const std::string &phase)
{
    auto result = network + "." + station + "." + phase;
    return result;
}

struct PickTablePair
{
    /// The travel IDs of the first and second table.
    std::pair<int, int> tableIDs;
    /// The index of the first and second pick in mPicks.
    std::pair<int, int> pickIndices;
};
/*!
 * @brief Makes the travel time pairs.
 * @param ids   Travel time table IDs.
 * @return A list of picks to migrate.
 */
template<class T>
std::vector<struct PickTablePair> makePairs(const std::vector<T> &ids)
{
    auto n = static_cast<int> (ids.size());
    auto nPairs = (n*(n - 1))/2;
    std::vector<PickTablePair> pairs;
    PickTablePair pair;
    pairs.reserve(nPairs);
    for (int i=0; i<n; ++i)
    {
        auto id = ids[i];
        for (int j=i+1; j<n; ++j)
        {
            auto jd = ids[j];
            // For example, a station may get two subsequent P wave triggers.
            // We don't want to migrate that pair since those picks cannot
            // be associated.
            if (id != jd)
            {
                pair.tableIDs = std::make_pair(id, jd);
                pair.pickIndices = std::make_pair(i, j);
                pairs.push_back(pair);
            }
        }
    }
    return pairs;
}

std::vector<struct PickTablePair>
    makePairs(const std::vector<MAssociate::Pick> &picks,
              const std::vector<int> &ids)
{
    auto n = static_cast<int> (ids.size());
    auto nPairs = (n*(n - 1))/2;
    std::vector<PickTablePair> pairs;
    PickTablePair pair;
    pairs.reserve(nPairs);
    for (int i=0; i<n; ++i)
    {
        auto id = ids[i];
        auto iOnsetTime = picks[i].getTime();
        auto iPhase = picks[i].getPhaseName();
        auto iSNCL = picks[i].getWaveformIdentifier();
        for (int j=i+1; j<n; ++j)
        {
            auto jd = ids[j];
            auto jOnsetTime = picks[j].getTime(); 
            auto jPhase = picks[j].getPhaseName();
            auto jSNCL = picks[j].getWaveformIdentifier();
            // Causality checks:
            //  (1) id == jd means that the station/phase match.  In this case, 
            //            the picks cannot be associated.
            if (id == jd){continue;}
            // An S on a station wave cannot show up prior to a P wave
            if (iSNCL.getNetwork() == jSNCL.getNetwork() &&
                iSNCL.getStation() == jSNCL.getStation())
            {
                if (iPhase == jPhase){continue;} // Should be handled by id==jd
                // P arrives after S
                if (iPhase == "P" && jPhase == "S")
                {
                    if (iOnsetTime >= jOnsetTime){continue;}    
                }
                // S arrives before P 
                if (iPhase == "S" && jPhase == "P")
                {
                    if (iOnsetTime <= jOnsetTime){continue;}
                }
            }
            pair.tableIDs = std::make_pair(id, jd);
            pair.pickIndices = std::make_pair(i, j);
            pairs.push_back(pair); 
        }
    }
    return pairs;
}

/*!
 * @brief Argument sort in descending order.
 * @param[in] v  The values to argument sort.
 * @result A permutation array such that v[perm[:]] is in descending order.
 */
/*
template <typename T>
std::vector<int> argSortDescending(const std::vector<T> &v)
{
    // initialize original index locations
    std::vector<int> perm(v.size());
    std::iota(perm.begin(), perm.end(), 0);
    // Sort indexes based on comparing values in v using std::stable_sort
    // instead of std::sort to avoid unnecessary index re-orderings
    // when v contains elements of equal values 
    std::stable_sort(perm.begin(), perm.end(),
                     [&v](const int i1, const int i2)
                     {
                          return v[i1] > v[i2];
                     });
    return perm;
}
*/

}

template<class T>
class Migrate<T>::MigrateImpl
{
public:
    [[nodiscard]] size_t getTableIndex(const size_t tile,
                                       const size_t tableID) const
    {
        auto result
           = tile*static_cast<size_t> (mNumberOfTables)
                 *static_cast<size_t> (mPaddedTileSize)
           + tableID*static_cast<size_t> (mPaddedTileSize);
        return result;
    }
    [[nodiscard]] int getTableIdentifier(const std::string &network,
                                         const std::string &station,
                                         const std::string &phase) const
    {
        auto name = makeTableName(network, station, phase);
        auto it = std::find(mTableNames.begin(),mTableNames.end(), name);
        if (it == mTableNames.end()){return -1;}
        auto id = static_cast<int> (std::distance(mTableNames.begin(), it));
        return id;
    }
    /// Holds the travel time tables.  This has dimension:
    /// [mNumberOfTiles x mNumberOfTables x mPaddedTileSize]
    [[maybe_unused]]
    std::vector<T, boost::alignment::aligned_allocator<T, 64>> mTables;
    [[maybe_unused]]
    std::vector<T, boost::alignment::aligned_allocator<T, 64>> mImage;
    [[maybe_unused]] std::vector<std::string> mTableNames;
    /// Contains the pick information and table identifier
    [[maybe_unused]] std::vector<MAssociate::Pick> mPicks;
    // For the ia'th pick this returns the corresponding table identifier.
    [[maybe_unused]] std::vector<int> mPickTableIdentifier;
    // The ia'th's picks contribution to the maximum objective function.
    [[maybe_unused]] std::vector<double> mContribution;
    // The travel time from the location of the objective function's
    // maximum to the ia'th pick.
    [[maybe_unused]] std::vector<double> mTravelTimesToMaximum;
    // Defines the connectivity between picks 
    std::vector<double> mSimilarityMatrix;
    std::vector<std::pair<int, int>> mAdjacency;
    std::vector<double> mEdgeWeights;
    MigrationParameters mParameters;
    T mMaxDifferentialTravelTime =-1;
    [[maybe_unused]] int mNumberOfTiles = 0;
    [[maybe_unused]] int mNumberOfTables = 0;
    [[maybe_unused]] int mTileSize = 0;
    [[maybe_unused]] int mPaddedTileSize = 0;
    [[maybe_unused]] int mMaxIndex =-1;
    [[maybe_unused]] bool mHaveImage = false;
    [[maybe_unused]] bool mInitialized = false;
};

/// Constructor
template<class T>
Migrate<T>::Migrate() :
    pImpl(std::make_unique<MigrateImpl> ())
{
}

/// Copy c'tor
template<class T>
[[maybe_unused]]
Migrate<T>::Migrate(const Migrate &migrate)
{
    *this = migrate;
}

/// Move c'tor
template<class T>
[[maybe_unused]]
Migrate<T>::Migrate(Migrate &&migrate) noexcept
{
    *this = std::move(migrate);
}

/// Copy assignment
template<class T>
Migrate<T>& Migrate<T>::operator=(const Migrate<T> &migrate)
{
    if (&migrate == this){return *this;}
    pImpl = std::make_unique<MigrateImpl> (*migrate.pImpl);
    return *this;
}

/// Move assignment
template<class T>
Migrate<T>& Migrate<T>::operator=(Migrate<T> &&migrate) noexcept
{
    if (&migrate == this){return *this;}
    pImpl = std::move(migrate.pImpl);
    return *this;
}

template<class T>
void Migrate<T>::clear() noexcept
{
    pImpl->mParameters.clear();
    pImpl->mTables.clear();
    pImpl->mImage.clear();
    pImpl->mTableNames.clear();
    pImpl->mPicks.clear();
    pImpl->mPickTableIdentifier.clear();
    pImpl->mSimilarityMatrix.clear();
    pImpl->mAdjacency.clear();
    pImpl->mEdgeWeights.clear();
    pImpl->mTravelTimesToMaximum.clear();
    pImpl->mContribution.clear();
    pImpl->mMaxDifferentialTravelTime =-1;
    pImpl->mNumberOfTiles = 0;
    pImpl->mNumberOfTables = 0;
    pImpl->mTileSize = 0;
    pImpl->mPaddedTileSize = 0;
    pImpl->mMaxIndex =-1;
    pImpl->mHaveImage = false;
    pImpl->mInitialized = false;
}

/// Destructor
template<class T>
Migrate<T>::~Migrate() = default;

template<class T>
void Migrate<T>::initialize(const MigrationParameters &parameters)
{
    // Reset the class
    clear();
    // Check the inputs
    if (!parameters.haveNumberOfTravelTimeTables())
    {
        throw std::invalid_argument("Number of tables not set");
    }
    if (!parameters.haveNumberOfPointsInTravelTimeTable())
    {
        throw std::invalid_argument("Number of points in table not set");
    }
    // Determine how much space to allocate
    auto nPoints = parameters.getNumberOfPointsInTravelTimeTable(); // Throws
    auto nTables = parameters.getNumberOfTravelTimeTables(); // Throws
    auto tileSize = parameters.getTileSize();
    pImpl->mNumberOfTables = nTables;
    pImpl->mPaddedTileSize = padLength(tileSize, sizeof(T), ALIGNMENT);
    pImpl->mTileSize = tileSize;
    pImpl->mNumberOfTiles
        = static_cast<int> (std::ceil(static_cast<double> (nPoints)
                                     /static_cast<double> (tileSize)));
    if (pImpl->mTileSize*pImpl->mNumberOfTiles < nPoints)
    {
        throw std::runtime_error("Algorithm failure");
    }
    auto workSpace = static_cast<size_t> (pImpl->mNumberOfTiles)
                    *static_cast<size_t> (pImpl->mNumberOfTables)
                    *static_cast<size_t> (pImpl->mPaddedTileSize);
    pImpl->mTables.resize(workSpace, 0);
    pImpl->mImage.resize(tileSize*pImpl->mNumberOfTiles, 0);
    pImpl->mTableNames.reserve(pImpl->mNumberOfTables);
    // Checks out
    pImpl->mParameters = parameters;
    pImpl->mInitialized = true;
}

/// Sets a travel time table
template<class T>
template<typename U>
[[maybe_unused]]
void Migrate<T>::setTravelTimeTable(const std::string &network,
                                    const std::string &station,
                                    const std::string &phase,
                                    const int nPoints, const U *times)
{
    auto nPts = getNumberOfPointsInTravelTimeTable(); //Throws if uninitialized
    pImpl->mMaxDifferentialTravelTime =-1;
    pImpl->mHaveImage = false;
    if (nPts != nPoints)
    {
        throw std::invalid_argument("nPoints = " + std::to_string(nPoints)
                                 + " must equal " + std::to_string(nPts));
    }
    if (times == nullptr){throw std::invalid_argument("Times is NULL");}
    if (network.empty()){throw std::invalid_argument("Network is empty");}
    if (station.empty()){throw std::invalid_argument("Station is empty");}
    if (phase.empty()){throw std::invalid_argument("Phase is empty");}
    auto minTime = std::min_element(times, times + nPoints);
    if (*minTime < 0)
    {
        throw std::invalid_argument("Negative travel time detected");
    }
    // Figure out the tableID
    int tableID = pImpl->getTableIdentifier(network, station, phase);
    auto tableName = makeTableName(network, station, phase);
    //
    if (tableID >= 0)
    {
        std::cerr << "Overwriting table: " << tableName << std::endl;
    }
    else
    {
        if (haveAllTravelTimeTables())
        {
            auto mTables
                = static_cast<size_t>(
                     pImpl->mParameters.getNumberOfTravelTimeTables());
            throw std::invalid_argument("All tables set - max is: "
                                       + std::to_string(mTables) );
        }
        pImpl->mTableNames.push_back(tableName);
        tableID = pImpl->getTableIdentifier(network, station, phase);
    }
    // Copy the table tile by tile
    auto numberOfTiles = pImpl->mNumberOfTiles;
    auto tileSize = pImpl->mTileSize;
    for (int tile=0; tile<numberOfTiles; ++tile)
    {
        // Source
        auto it1 = tile*tileSize;
        auto it2 = std::min(nPts, it1 + tileSize);
        // Destination
        auto jt1 = pImpl->getTableIndex(tile, tableID);
        T *tPtr __attribute__((aligned(64))) = pImpl->mTables.data() + jt1;
        std::copy(times + it1, times + it2, tPtr);
    }
}

/// Gets a travel time 
template<class T>
T Migrate<T>::getTravelTime(const std::string &network,
                            const std::string &station,
                            const std::string &phase, const int index) const
{
    auto nPts = getNumberOfPointsInTravelTimeTable(); //Throws if uninitialized
    if (!haveTravelTimeTable(network, station, phase))
    {
        throw std::invalid_argument("Table does not exist for: "
                                  + network + "," + station + "," + phase);
    }
    if (index < 0 || index >= nPts)
    {
        throw std::invalid_argument("index = " + std::to_string(index)
                                  + " must be in the range [0,"
                                  + std::to_string(nPts - 1) + "]");
    }
    auto tableID = pImpl->getTableIdentifier(network, station, phase);
    int tile = index/pImpl->mTileSize;
    int gridIndexInTile = index - tile*pImpl->mTileSize;
    assert(tile*pImpl->mTileSize + gridIndexInTile == index);
    auto it = pImpl->getTableIndex(tile, tableID);
    T tEst = pImpl->mTables.at(it + gridIndexInTile);
    return tEst;
}

/// Gets a specific travel time table
template<class T>
std::vector<T> Migrate<T>::getTravelTimeTable(const std::string &network,
                                              const std::string &station,
                                              const std::string &phase) const
{
    auto nPts = getNumberOfPointsInTravelTimeTable(); //Throws if uninitialized
    if (!haveTravelTimeTable(network, station, phase))
    {
        throw std::invalid_argument("Table does not exist for: "
                                  + network + "," + station + "," + phase);
    }
    // Copy table tile by tile
    auto tableID = pImpl->getTableIdentifier(network, station, phase);
    auto numberOfTiles = pImpl->mNumberOfTiles;
    auto tileSize = pImpl->mTileSize;
    std::vector<T> times(nPts);
    for (int tile=0; tile<numberOfTiles; ++tile)
    {
        // Source
        auto it1 = pImpl->getTableIndex(tile, tableID);
        // Destination
        auto jt1 = tile*tileSize;
        auto jt2 = std::min(nPts, jt1 + tileSize);
        auto nCopy = jt2 - jt1;
        T *tPtr __attribute__((aligned(64))) = pImpl->mTables.data() + it1;
        std::copy(tPtr, tPtr + nCopy,  times.data() + jt1);
    }
    return times;
}

/// Have the travel time table?
template<class T>
bool Migrate<T>::haveTravelTimeTable(const std::string &network,
                                     const std::string &station,
                                     const std::string &phase) const noexcept
{
    if (!isInitialized()){return false;}
    auto id = pImpl->getTableIdentifier(network, station, phase);
    return (id >= 0);
}

/// Get the current max differential travel time
template<class T>
T Migrate<T>::getMaximumDifferentialTravelTime() const noexcept
{
    T dtMax = 0;
    if (!isInitialized()){return dtMax;}
    auto nPts = getNumberOfPointsInTravelTimeTable();
    if (pImpl->mMaxDifferentialTravelTime >= 0)
    {
        return pImpl->mMaxDifferentialTravelTime;
    }
    auto nTables = static_cast<int> (pImpl->mTableNames.size());
    if (nTables < 2){return dtMax;} // Can't compute anything
    auto numberOfTiles = pImpl->mNumberOfTiles;
    auto tileSize = pImpl->mTileSize;
    auto mTables = pImpl->mTables.data();
    for (int tile=0; tile<numberOfTiles; ++tile)
    {
        auto jt1 = tile*tileSize;
        auto jt2 = std::min(nPts, jt1 + tileSize);
        auto nLocal = jt2 - jt1;
        for (int it=0; it<nTables; ++it)
        {
            for (int jt=it+1; jt<nTables; ++jt)
            {
                // Source
                auto it1 = pImpl->getTableIndex(tile, it);
                auto it2 = pImpl->getTableIndex(tile, jt);
                // Destination
                const T *tPtr1 __attribute__((aligned(64))) = mTables + it1;
                const T *tPtr2 __attribute__((aligned(64))) = mTables + it2;
                #pragma omp simd reduction(max:dtMax) aligned(tPtr1, tPtr2: 64)
                for (int i=0; i<nLocal; ++i)
                {
                    dtMax = std::max(dtMax, std::abs(tPtr1[i] - tPtr2[i]));
                }
            }
        }
    }
    pImpl->mMaxDifferentialTravelTime = dtMax;
    return pImpl->mMaxDifferentialTravelTime;
}

/// Have all the tables?
template<class T>
bool Migrate<T>::haveAllTravelTimeTables() const noexcept
{
    if (!isInitialized()){return false;}
    return (static_cast<int> (pImpl->mTableNames.size()) ==
            pImpl->mNumberOfTables);
}

/// Picks
template<class T>
void Migrate<T>::addPick(const Pick &newPick)
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    pImpl->mHaveImage = false;
    // Check pick has a corresponding travel time table and pick time
    auto newSNCL = newPick.getWaveformIdentifier();
    auto networkNew = newSNCL.getNetwork();
    auto stationNew = newSNCL.getStation();
    auto phaseNew = newPick.getPhaseName();
    auto id = pImpl->getTableIdentifier(networkNew, stationNew, phaseNew);
    if (id < 0)
    {
        auto errmsg = "Travel time table for: "
                    + makeTableName(networkNew, stationNew, phaseNew)
                    + " does not exist";
        throw std::invalid_argument(errmsg);
    }
    if (!newPick.haveTime())
    {
        throw std::invalid_argument("Pick does not have time");
    }
    pImpl->mPicks.push_back(newPick);
    pImpl->mPickTableIdentifier.push_back(id);
}

template<class T>
void Migrate<T>::clearPicks() noexcept
{
    pImpl->mPicks.clear();
    pImpl->mPickTableIdentifier.clear();
}

template<class T>
int Migrate<T>::getNumberOfPicks() const noexcept
{
    return static_cast<int> (pImpl->mPicks.size());
}

/// Do the hard work and migrate
template<class T>
void Migrate<T>::migrate()
{
    const T sqrt12 = 3.464101615137754587054892683011744733;
    if (!isInitialized())
    { throw std::runtime_error("Class not initialized"); }
    // Null out result
    pImpl->mHaveImage = false;
    pImpl->mMaxIndex =-1;
    std::fill(pImpl->mImage.begin(), pImpl->mImage.end(), 0);
    auto nPicks = getNumberOfPicks();
    if (nPicks < 2)
    {
        throw std::runtime_error("At least two picks required");
    }
    std::vector<int> tableIDs(nPicks);
    for (int ia=0; ia<nPicks; ++ia)
    {
        tableIDs[ia] = pImpl->mPickTableIdentifier[ia];
    }
    // Figure out the travel time table index pairs
    auto pickTablePairs = makePairs(pImpl->mPicks,
                                       pImpl->mPickTableIdentifier); //tableIDs);
    auto nPairs = static_cast<int> (pickTablePairs.size());
    if (nPairs < 1)
    {
        std::cerr << "No pick pairs to migrate" << std::endl;
        return;
    }
    // Get some other constants needed for the loop
    auto method = pImpl->mParameters.getAnalyticCorrelationFunction();
    auto nPts = getNumberOfPointsInTravelTimeTable();
    auto tileSize = pImpl->mParameters.getTileSize();
    auto numberOfTiles = pImpl->mNumberOfTiles;
    const T *mTables = pImpl->mTables.data();
    // Loop on the tiles.  This prevents us from having each thread allocate
    // space for an image which would have to be reduced.
    #pragma omp parallel for \
     shared(pickTablePairs, method, mTables, numberOfTiles) \
     shared(nPicks, nPairs, nPts, tileSize) \
     default(none)
    for (int tile = 0; tile < numberOfTiles; ++tile)
    {
        // Destination
        auto jt1 = tile*tileSize;
        auto jt2 = std::min(nPts, jt1 + tileSize);
        auto nLocal = jt2 - jt1;
        auto imageLocal = pImpl->mImage.data() + jt1;
        // Loop on travel times pairs to migrate
        for (int iPair=0; iPair<nPairs; ++iPair)
        {
            auto t1 = pickTablePairs[iPair].tableIDs.first;
            auto t2 = pickTablePairs[iPair].tableIDs.second;
            auto ip1 = pickTablePairs[iPair].pickIndices.first;
            auto ip2 = pickTablePairs[iPair].pickIndices.second;
            T p1 = pImpl->mPicks[ip1].getTime();
            T std1 = pImpl->mPicks[ip1].getStandardDeviation();
            T tCorr1 = pImpl->mPicks[ip1].getStaticCorrection();
            T p2 = pImpl->mPicks[ip2].getTime();
            T std2 = pImpl->mPicks[ip2].getStandardDeviation();
            T tCorr2 = pImpl->mPicks[ip2].getStaticCorrection();
            // Source
            auto it1 = pImpl->getTableIndex(tile, t1);
            auto it2 = pImpl->getTableIndex(tile, t2);
            const T *tPtr1 __attribute__((aligned(64))) = mTables + it1;
            const T *tPtr2 __attribute__((aligned(64))) = mTables + it2;
            // Perform Gaussian migration
            if (method == AnalyticCorrelationFunction::GAUSSIAN)
            {
                T normalizeAmp
                    = analyticGaussianCorrelationAmplitudeNormalization(std1, std2);
                T normalizeExp
                    = analyticGaussianCorrelationExponentNormalization(std1, std2);
                #pragma omp simd aligned(tPtr1, tPtr2 : 64)
                for (int i = 0; i < nLocal; ++i)
                {
                    auto deltaT = (tPtr2[i] + tCorr2)
                                - (tPtr1[i] + tCorr1);
                    imageLocal[i] = imageLocal[i]
                        + analyticGaussianCorrelation(deltaT, p1, p2,
                                                      normalizeAmp,
                                                      normalizeExp);
                     //if (jt1 + i == 13873){std::cout << deltaT << "," << p2 - p1 << std::endl;}
                } // Loop on points in tile
            }
            else
            {
                T w1 = sqrt12*std1;
                T w2 = sqrt12*std2;
                //std::cout << w1 << "," << w2 << std::endl;
                #pragma omp simd aligned(tPtr1, tPtr2 : 64)
                for (int i = 0; i < nLocal; ++i)
                {
                    auto deltaT = (tPtr2[i] + tCorr2)
                                - (tPtr1[i] + tCorr1);
                    imageLocal[i] = imageLocal[i]
                        + analyticBoxcarCorrelation(deltaT, p1, p2, w1, w2);
                     //if (jt1 + i == 13873){std::cout << deltaT << "," << p2 - p1 << std::endl;}
//if (jt1 +i == 13873){std::cout << analyticBoxcarCorrelation(deltaT, p1, p2, w1, w2) << "," << w1 << ", " << w2 << "," << deltaT << "," << p2 - p1 << std::endl;}
                } // Loop on points in tile
            }
        } // Loop on pairs
    } // Loop on tiles
    // Get the maximum point
    auto xMax = std::max_element(pImpl->mImage.begin(), pImpl->mImage.end());
    int maxIndex = std::distance(pImpl->mImage.begin(), xMax);
    //std::cout << maxIndex << "," << *xMax << std::endl;
    int tile = maxIndex/tileSize;
    int optimumGridInTile = maxIndex - tile*tileSize;
    assert(tile*tileSize + optimumGridInTile == maxIndex);
    // Figure out how much each pick contributed to the max 
    std::vector<T> contribution(nPicks, 0);
    std::vector<double> travelTime(nPicks, 0);
    std::vector<bool> haveTravelTime(nPicks, false);
    pImpl->mSimilarityMatrix.resize(nPicks*nPicks, 0);
    pImpl->mAdjacency.resize(nPairs);
    pImpl->mEdgeWeights.resize(nPairs, 0);
    T optSum = 0;
    T half = 0.5;
    for (int iPair=0; iPair<nPairs; ++iPair)
    {
        auto t1 = pickTablePairs[iPair].tableIDs.first;
        auto t2 = pickTablePairs[iPair].tableIDs.second;
        auto ip1 = pickTablePairs[iPair].pickIndices.first;
        auto ip2 = pickTablePairs[iPair].pickIndices.second;
        T p1 = pImpl->mPicks[ip1].getTime();
        T std1 = pImpl->mPicks[ip1].getStandardDeviation();
        T tCorr1 = pImpl->mPicks[ip1].getStaticCorrection();
        T p2 = pImpl->mPicks[ip2].getTime();
        T std2 = pImpl->mPicks[ip2].getStandardDeviation();
        T tCorr2 = pImpl->mPicks[ip2].getStaticCorrection();
        // Source
        auto it1 = pImpl->getTableIndex(tile, t1);
        auto it2 = pImpl->getTableIndex(tile, t2);
        T tEst1 = mTables[it1 + optimumGridInTile] + tCorr1;
        T tEst2 = mTables[it2 + optimumGridInTile] + tCorr2;
        auto deltaT = tEst2 - tEst1;
        travelTime[ip1] = tEst1;
        travelTime[ip2] = tEst2;
        haveTravelTime[ip1] = true;
        haveTravelTime[ip2] = true;
        T pdf = 0;
        if (method == AnalyticCorrelationFunction::GAUSSIAN)
        {
            T normalizeAmp
               = analyticGaussianCorrelationAmplitudeNormalization(std1, std2);
            T normalizeExp
               = analyticGaussianCorrelationExponentNormalization(std1, std2);
            pdf = analyticGaussianCorrelation(deltaT, p1, p2,
                                              normalizeAmp,
                                              normalizeExp);
            contribution[ip1] = contribution[ip1] + pdf;
            contribution[ip2] = contribution[ip2] + pdf;
        }
        else
        {
            T w1 = sqrt12*std1;
            T w2 = sqrt12*std2;
            pdf = analyticBoxcarCorrelation(deltaT, p1, p2, w1, w2);
        }
        pImpl->mAdjacency[iPair] = std::make_pair(ip1, ip2);
        pImpl->mEdgeWeights[iPair] = pdf;
        // Both pick 1 and pick 2 contribute equally.  
        // This ensures that sum(contributions) = value in mImage[maxIndex] 
        auto indx = ip1*nPicks + ip2;
        auto jndx = ip2*nPicks + ip1;
        pImpl->mSimilarityMatrix[indx] = pdf;
        pImpl->mSimilarityMatrix[jndx] = pdf;
        contribution[ip1] = contribution[ip1] + half*pdf;
        contribution[ip2] = contribution[ip2] + half*pdf;
        optSum = optSum + pdf;
    }
    pImpl->mContribution.resize(contribution.size());
    std::copy(contribution.begin(), contribution.end(),
              pImpl->mContribution.begin());
    pImpl->mTravelTimesToMaximum = travelTime;
    // Create the main diagonal
    for (int ip=0; ip<nPicks; ++ip)
    {
        T std1 = pImpl->mPicks[ip].getStandardDeviation();
        T pdf = 0;
        T zero = 0;
        if (method == AnalyticCorrelationFunction::GAUSSIAN)
        {
            T normalizeAmp
               = analyticGaussianCorrelationAmplitudeNormalization(std1, std1);
            T normalizeExp
               = analyticGaussianCorrelationExponentNormalization(std1, std1);
            pdf = analyticGaussianCorrelation(zero, zero, zero,
                                              normalizeAmp, normalizeExp);
        }
        else
        {
            T w1 = sqrt12*std1;
            pdf = analyticBoxcarCorrelation(zero, zero, zero, w1, w1);
        }
        auto indx = ip*nPicks + ip;
        //std::cout << pdf << ","  << indx << std::endl;
        pImpl->mSimilarityMatrix[indx] = pdf;
    }

    // Normalize and sort into descending order
    if (optSum == 0){optSum = 1;}
    for (auto &c : contribution){c = c/optSum;}
/*
    auto perm = argSortDescending(contribution);
    // Perform some basic clustering
    std::vector<double> originTimes;
    std::vector<double> weights;
    originTimes.reserve(nPicks);
    weights.reserve(nPicks);
    for (int ia=0; ia<nPicks; ++ia)
    {
        if (haveTravelTime[ia])
        {
            originTimes.push_back(pImpl->mPicks[ia].getTime()
                                - travelTime[ia]);
            weights.push_back(std::max(std::numeric_limits<T>::epsilon(),
                                       contribution[ia]));
        }
    }
    // Cluster
    DBSCAN dbscan;
    dbscan.initialize(pImpl->mParameters.getDBSCANEpsilon(),
                      pImpl->mParameters.getDBSCANMinimumClusterSize());
    dbscan.setWeightedData(originTimes.size(), 1,
                           originTimes.data(), weights.data());
    dbscan.cluster();
    auto nClusters = dbscan.getNumberOfClusters();
    if (nClusters > 0)
    {
        auto labels = dbscan.getLabels();
        if (nClusters == 1)
        {
            for (int ia=0; ia<nPicks; ++ia)
            {
                if (labels[ia] >-1)
                {
                }
            } 
        }
        else
        {
            auto uniqueLabels = labels;
            std::sort(uniqueLabels.begin(), uniqueLabels.end());
            auto it = std::unique(uniqueLabels.begin(), uniqueLabels.end());
            uniqueLabels.resize(std::distance(uniqueLabels.begin(), it));
       }
    }
for (auto &ot : originTimes){std::cout << ot << std::endl;}
std::cout << contribution[perm[0]] << std::endl;
std::cout << nClusters << "," << optSum << std::endl;
*/
    pImpl->mMaxIndex = maxIndex;
    pImpl->mHaveImage = true;
}

/// Gets the max index in the image
template<class T>
std::pair<int, T> Migrate<T>::getImageMaximum() const
{
    if (!haveImage())
    {
        throw std::invalid_argument("Migration image not yet computed");
    }
    return std::pair(pImpl->mMaxIndex, pImpl->mImage[pImpl->mMaxIndex]);
}

/// Gets the travel times to all picks from the maximum index
template<class T>
std::vector<double> Migrate<T>::getTravelTimesToMaximum() const
{
    if (!haveImage()){throw std::runtime_error("Image not yet computed");}
    return pImpl->mTravelTimesToMaximum;
}

/// Gets the picks contribution to maximum index
template<class T>
std::vector<double>
    Migrate<T>::getContributionToMaximum(const bool normalize) const
{
    if (!haveImage()){throw std::runtime_error("Image not yet computed");}
    if (normalize)
    {
        const double zero = 0;
        std::vector<double> result(pImpl->mContribution.size(), zero);
        auto xnorm = std::accumulate(pImpl->mContribution.begin(),
                                     pImpl->mContribution.end(), zero);
        if (xnorm > 0)
        {
            xnorm = 1./xnorm;
            std::transform(pImpl->mContribution.begin(),
                           pImpl->mContribution.end(),
                           result.begin(),
                           [xnorm](const double x) -> double
                           {
                               return x*xnorm;
                           });
        }
        return result;
    }  
    return pImpl->mContribution;
}


/// Migration image computed?
template<class T>
bool Migrate<T>::haveImage() const noexcept
{
    return pImpl->mHaveImage;
}

/// Initialized?
template<class T>
bool Migrate<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

template<class T>
int Migrate<T>::getNumberOfPointsInTravelTimeTable() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mParameters.getNumberOfPointsInTravelTimeTable();
}

///--------------------------------------------------------------------------///
///                         Instantiate the Templates                        ///
///--------------------------------------------------------------------------///
template class MAssociate::Migrate<double>;
template class MAssociate::Migrate<float>;

template void MAssociate::Migrate<double>::setTravelTimeTable(
    const std::string &network, const std::string &station,
    const std::string &phase, int nPoints, const double times[]);
template void MAssociate::Migrate<double>::setTravelTimeTable(
    const std::string &network, const std::string &station,
    const std::string &phase, int nPoints, const float times[]);

template void MAssociate::Migrate<float>::setTravelTimeTable(
    const std::string &network, const std::string &station,
    const std::string &phase, int nPoints, const double times[]);
template void MAssociate::Migrate<float>::setTravelTimeTable(
    const std::string &network, const std::string &station,
    const std::string &phase, int nPoints, const float times[]);
