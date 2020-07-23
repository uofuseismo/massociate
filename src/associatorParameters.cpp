#include "massociate/associatorParameters.hpp"
using namespace MAssociate;

class AssociatorParameters::AssociatorParametersImpl
{
public:
    double mDBSCANEpsilon = 0.2;
    double mPageRankDamping = 0.85;
    int mArrivalsToNucleate = 6;
    int mDBSCANMinClusterSize = 6;
    int mPageRankIterations = 20;
    int mTables =-1;
    int mTileSize = 1024;
    AnalyticCorrelationFunction mAnalyticFunction
        = AnalyticCorrelationFunction::BOXCAR;
};

/// C'tor
AssociatorParameters::AssociatorParameters() :
    pImpl(std::make_unique<AssociatorParametersImpl>())
{
}

/// Copy c'tor
AssociatorParameters::AssociatorParameters(
    const AssociatorParameters &parameters)
{
    *this = parameters;
}

/// Move c'tor
[[maybe_unused]]
AssociatorParameters::AssociatorParameters(
    AssociatorParameters &&parameters) noexcept
{
    *this = std::move(parameters);
}

/// Copy assignment
AssociatorParameters&
AssociatorParameters::operator=(const AssociatorParameters &parameters)
{
    if (&parameters == this){return *this;}
    pImpl = std::make_unique<AssociatorParametersImpl> (*parameters.pImpl);
    return *this;
}

/// Move assignment
AssociatorParameters&
AssociatorParameters::operator=(AssociatorParameters &&parameters) noexcept
{
    if (&parameters == this){return *this;}
    pImpl = std::move(parameters.pImpl);
    return *this;
}

/// Destructor
AssociatorParameters::~AssociatorParameters() = default;

/// Resets the class
void AssociatorParameters::clear() noexcept
{
    pImpl->mDBSCANEpsilon = 0.2;
    pImpl->mPageRankDamping = 0.85;
    pImpl->mArrivalsToNucleate = 6;
    pImpl->mDBSCANMinClusterSize = 6;
    pImpl->mPageRankIterations = 20;
    pImpl->mTables =-1;
    pImpl->mTileSize = 1024;
    pImpl->mAnalyticFunction = AnalyticCorrelationFunction::BOXCAR;
}

/// Sets/gets number of tables
void AssociatorParameters::setNumberOfTravelTimeTables(const int nTables)
{
    if (nTables < 2)
    {
        throw std::invalid_argument("At least 2 tables required");
    }
    pImpl->mTables = nTables;
}

int AssociatorParameters::getNumberOfTravelTimeTables() const
{
    if (!haveNumberOfTravelTimeTables())
    {
        throw std::runtime_error("Number of travel time tables not yet set");
    }
    return pImpl->mTables;
}

bool AssociatorParameters::haveNumberOfTravelTimeTables() const noexcept
{
    return (pImpl->mTables > 0);
}


/// Sets/gets minimum number of arrivals to nucleate
void AssociatorParameters::setMinimumNumberOfArrivalsToNucleate(
    const int minArrivals)
{
    if (minArrivals < 1)
    {
        throw std::invalid_argument("minArrivals must be positive");
    }
    pImpl->mArrivalsToNucleate = minArrivals;
}

int AssociatorParameters::getMinimumNumberOfArrivalsToNucleate() const noexcept
{
    return pImpl->mArrivalsToNucleate;
}

/// Sets/gets DBSCAN epsilon
void AssociatorParameters::setDBSCANEpsilon(const double epsilon)
{
    if (epsilon < 0){throw std::invalid_argument("epsilon cannot be negative");}
    pImpl->mDBSCANEpsilon = epsilon;
}

double AssociatorParameters::getDBSCANEpsilon() const noexcept
{
    return pImpl->mDBSCANEpsilon;
}

/// Sets gets DBSCAN min cluster size
void AssociatorParameters::setDBSCANMinimumClusterSize(
    const int minimumClusterSize)
{
    if (minimumClusterSize < 2)
    {
        throw std::invalid_argument("minimum cluster size must be at least 2");
    }
    pImpl->mDBSCANMinClusterSize = minimumClusterSize;
}

int AssociatorParameters::getDBSCANMinimumClusterSize() const noexcept
{
    return pImpl->mDBSCANMinClusterSize;
}

/// Sets/gets PageRank damping
void AssociatorParameters::setPageRankDampingFactor(const double damping)
{
    if (damping < 0 || damping > 1)
    {
        throw std::invalid_argument("Damping must be in range [0,1]");
    }
    pImpl->mPageRankDamping = damping;
}

double AssociatorParameters::getPageRankDampingFactor() const noexcept
{
    return pImpl->mPageRankDamping;
}

/// Sets/gets PageRank iterations
void AssociatorParameters::setPageRankNumberOfIterations(
    const int nIterations)
{
    if (nIterations < 0)
    {
        throw std::invalid_argument(
            "Number of power iterations must be positive");
    }
    pImpl->mPageRankIterations = nIterations;
}

int AssociatorParameters::getPageRankNumberOfIterations() const noexcept
{
    return pImpl->mPageRankIterations;
}

/// Sets/gets the tile size
void AssociatorParameters::setTileSize(const int tileSize)
{
    if (tileSize < 1)
    {
        throw std::invalid_argument("tileSize must be positive");
    }
    pImpl->mTileSize = tileSize;
}

int AssociatorParameters::getTileSize() const noexcept
{
    return pImpl->mTileSize;
}

/// Sets/gets correlation function
void AssociatorParameters::setAnalyticCorrelationFunction(
    MAssociate::AnalyticCorrelationFunction function) noexcept
{
    pImpl->mAnalyticFunction = function;
}

AnalyticCorrelationFunction
AssociatorParameters::getAnalyticCorrelationFunction() const noexcept
{
    return pImpl->mAnalyticFunction;
}
