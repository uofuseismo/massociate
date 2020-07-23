#include <iostream>
#include <algorithm>
#include "massociate/migrationParameters.hpp"

using namespace MAssociate;

class MigrationParameters::MigrationParametersImpl
{
public:
    int mTables =-1;
    int mPoints =-1;
    int mTileSize = 1024;
    AnalyticCorrelationFunction mAnalyticFunction
       = AnalyticCorrelationFunction::GAUSSIAN;
};

/// Constructor
MigrationParameters::MigrationParameters() :
    pImpl(std::make_unique<MigrationParametersImpl> ())
{
}

/// Copy c'tor
MigrationParameters::MigrationParameters(const MigrationParameters &parameters)
{
    *this = parameters;
}

/// Move c'tor
[[maybe_unused]]
MigrationParameters::MigrationParameters(
    MigrationParameters &&parameters) noexcept
{
    *this = std::move(parameters);
}

/// Destructor
MigrationParameters::~MigrationParameters() = default;

/// Resets the class
[[maybe_unused]] [[maybe_unused]]
void MigrationParameters::clear() noexcept
{
    pImpl->mPoints =-1;
    pImpl->mTables =-1;
    pImpl->mTileSize = 1024;
    pImpl->mAnalyticFunction = AnalyticCorrelationFunction::GAUSSIAN;
}

/// Copy assignment
MigrationParameters&
MigrationParameters::operator=(const MigrationParameters &parameters)
{
    if (&parameters == this){return *this;}
    pImpl = std::make_unique<MigrationParametersImpl> (*parameters.pImpl);
    return *this;
}

/// Move assignment
MigrationParameters&
MigrationParameters::operator=(MigrationParameters &&parameters) noexcept
{
    if (&parameters == this){return *this;}
    pImpl = std::move(parameters.pImpl);
    return *this;
}

/// Sets/gets number of points in the travel time table
void MigrationParameters::setNumberOfPointsInTravelTimeTable(const int nPoints)
{
    if (nPoints < 1)
    {
        throw std::invalid_argument("nPoints must be positive");
    }
    pImpl->mPoints = nPoints;
    pImpl->mTileSize = std::min(pImpl->mPoints, pImpl->mTileSize);
}

int MigrationParameters::getNumberOfPointsInTravelTimeTable() const
{
    if (!haveNumberOfPointsInTravelTimeTable())
    {
        throw std::runtime_error("Points in travel time table not yet set");
    }
    return pImpl->mPoints;
}

bool MigrationParameters::haveNumberOfPointsInTravelTimeTable() const noexcept
{
    return (pImpl->mPoints > 0);
}

/// Sets/gets number of tables
void MigrationParameters::setNumberOfTravelTimeTables(const int nTables)
{
    if (nTables < 2)
    {
        throw std::invalid_argument("At least 2 tables required");
    }
    pImpl->mTables = nTables;
}

int MigrationParameters::getNumberOfTravelTimeTables() const
{
    if (!haveNumberOfTravelTimeTables())
    {
        throw std::runtime_error("Number of travel time tables not yet set");
    }
    return pImpl->mTables;
}

bool MigrationParameters::haveNumberOfTravelTimeTables() const noexcept
{
    return (pImpl->mTables > 0);
}

/// Sets/gets the tile size
void MigrationParameters::setTileSize(const int tileSize)
{
    if (tileSize < 1)
    {
        throw std::invalid_argument("Tile size must be positive");
    }
    int nPoints = pImpl->mTileSize;
    try
    {
        nPoints = getNumberOfPointsInTravelTimeTable();
    }
    catch (const std::exception &e)
    {
        std::cerr << "WARNING: Number of points not yet set" << std::endl;
    }
    pImpl->mTileSize = std::min(nPoints, tileSize);
}

int MigrationParameters::getTileSize() const noexcept
{
    return pImpl->mTileSize;
}

/// Sets/gets correlation function
void MigrationParameters::setAnalyticCorrelationFunction(
    MAssociate::AnalyticCorrelationFunction function) noexcept
{
    pImpl->mAnalyticFunction = function;
}

AnalyticCorrelationFunction
MigrationParameters::getAnalyticCorrelationFunction() const noexcept
{
    return pImpl->mAnalyticFunction;
}

