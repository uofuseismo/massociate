#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include "massociate/associator.hpp"
#include "massociate/associatorParameters.hpp"
#include "massociate/migrationParameters.hpp"
#include "massociate/pick.hpp"
#include "massociate/waveformIdentifier.hpp"
#include "massociate/migrate.hpp"
#include "massociate/mesh/spherical/points3d.hpp"
#include "massociate/mesh/cartesian/points3d.hpp"

using namespace MAssociate;

namespace
{

}

template<class T>
class Associator<T>::AssociatorImpl
{
public:
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
    if (getNumberOfPicks() < 1)
    {
        std::cerr << "No picks" << std::endl;
    }
    // Sort the picks temporally

    // Perform initial migration
    pImpl->mMigrate.migrate();
    // Now cluster
}

///--------------------------------------------------------------------------///
///                          Template Instantiation                          ///
///--------------------------------------------------------------------------///
template class MAssociate::Associator<double>;
template class MAssociate::Associator<float>;
