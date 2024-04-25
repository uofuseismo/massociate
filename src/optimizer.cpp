#include "massociate/optimizer.hpp"
#include "massociate/migrator.hpp"

using namespace MAssociate;

class IOptimizer::IOptimizerImpl
{
public:
    std::unique_ptr<IMigrator> mMigrator{nullptr};
};

/// Constructor
IOptimizer::IOptimizer() :
    pImpl(std::make_unique<IOptimizerImpl> ())
{
}

/// Destructor
IOptimizer::~IOptimizer() = default;

/// Migrator
void IOptimizer::setMigrator(std::unique_ptr<IMigrator> &&migrator)
{
    if (!migrator->haveGeographicRegion())
    {
        throw std::invalid_argument("Geographic region not set on migrator");
    }
    if (!migrator->haveTravelTimeCalculatorMap())
    {
        throw std::invalid_argument(
           "Travel time calculator map not set on migrator");
    }
    pImpl->mMigrator = std::move(migrator);
}

std::unique_ptr<IMigrator> IOptimizer::releaseMigrator()
{
    if (!haveMigrator()){throw std::runtime_error("Migrator not set");}
    auto result = std::move(pImpl->mMigrator);
    pImpl->mMigrator = nullptr;
    return result;
}

IMigrator *IOptimizer::getMigratorHandle()
{
    if (!haveMigrator()){throw std::runtime_error("Migrator not set");}
    return pImpl->mMigrator.get();
}

bool IOptimizer::haveMigrator() const noexcept
{
    return pImpl->mMigrator != nullptr;
}

/// Sets the arrivals
void IOptimizer::setArrivals(const std::vector<Arrival> &arrivals)
{
    if (!haveMigrator()){throw std::runtime_error("Arrivals not set");}
    pImpl->mMigrator->setArrivals(arrivals);
}

bool IOptimizer::haveArrivals() const noexcept
{
    if (!haveMigrator()){return false;}
    return pImpl->mMigrator->haveArrivals();
}
