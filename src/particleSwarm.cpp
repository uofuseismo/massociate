#include <cmath>
#include <vector>
#include <limits>
#include <functional>
#include <pagmo/pagmo.hpp>
#include <uLocator/position/geographicRegion.hpp>
#include <uLocator/optimizers/originTime.hpp>
#include <umps/logging/standardOut.hpp>
#include "massociate/particleSwarm.hpp"
#include "massociate/migrator.hpp"
#include "massociate/arrival.hpp"
#include "massociate/event.hpp"

using namespace MAssociate;

namespace
{
struct FixedDepthDoubleDifferenceMigration
{
    pagmo::vector_double fitness(const pagmo::vector_double &dv) const
    {
        double fitness = std::numeric_limits<double>::max();
#ifndef NDEBUG
        assert(dv.size() == static_cast<int> (mParameters));
#endif
        try
        {
            // Make this negative b/c the migration operator looks to maximize
            // but Pagmo wants to minimize - i.e., the `best' image point is
            // the most negative image point
            fitness =-mMigrationFunction(dv.at(0), dv.at(1));
        }
        catch (const std::exception &e) 
        {
            auto errorMessage = "Problem for source at f.s. (x,y) = ("
                              + std::to_string(dv.at(0)) + "," 
                              + std::to_string(dv.at(1))
                              + ").  Failed with: "  + std::string {e.what()};
            std::cerr << errorMessage << std::endl;
        }
        return pagmo::vector_double {fitness};
    }   
    /// @result The number of equality constraints
    pagmo::vector_double::size_type get_eic() const
    {   
        return 0;
    }   
    /// @result The number of inequality constraints.
    pagmo::vector_double::size_type get_nic() const
    {   
        return 0;
    }   
    /// @result The hard bounds on the search region.
    std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const
    {
        return std::pair {mLowerBounds, mUpperBounds};
    }   
    /// @brief Sets the search boundaries.
    void setSearchBoundaries(
        const std::vector<double> &lowerBoundaries,
        const std::vector<double> &upperBoundaries)
    {   
        if (static_cast<int> (lowerBoundaries.size()) != mParameters)
        {
            throw std::invalid_argument("lowerBoundaries size != 2");
        }
        if (lowerBoundaries.size() != upperBoundaries.size())
        {
            throw std::invalid_argument("upperBoundaries size != 2");
        }
        for (int i = 0; i < static_cast<int> (lowerBoundaries.size()); ++i)
        {
            if (lowerBoundaries[i] >= upperBoundaries[i])
            {
                throw std::invalid_argument("l[i] >= u[i]");
            }
        }
        mLowerBounds = lowerBoundaries;
        mUpperBounds = upperBoundaries;
    } 

    MAssociate::IMigrator *mMigrator{nullptr}; 
    std::function< double(const double, const double) >
        mMigrationFunction;
    pagmo::vector_double mLowerBounds;
    pagmo::vector_double mUpperBounds;
    int mParameters{2}; 
};
}

class ParticleSwarm::ParticleSwarmImpl
{
public:
    ParticleSwarmImpl(std::shared_ptr<UMPS::Logging::ILog> logger = nullptr) :
        mLogger(logger)
    {
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> ();
        }
    }
    [[nodiscard]] double migrateFixedDepth(const double x, const double y) const
    {
        return mMigrator->evaluate(x, y, mDepth);
    }
    std::function<double (const double, const double)>
    mFixedDepthMigrationFunction
    {
        std::bind(&ParticleSwarmImpl::migrateFixedDepth,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2)
    };
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    std::unique_ptr<IMigrator> mMigrator{nullptr};
    Event mEvent;
    std::pair<double, double> mExtentInX;
    std::pair<double, double> mExtentInY;
    double mOptimalMigrationValue{0};
    double mDepth{6000};
    int mGenerations{50};
    int mParticles{20};
    int mMinimumNumberOfArrivalsToBuildEvent{4};
    bool mHaveExtentInX{false};
    bool mHaveExtentInY{false};
    bool mHaveEvent{false};
};

/// Constructor
ParticleSwarm::ParticleSwarm() :
    pImpl(std::make_unique<ParticleSwarmImpl> ())
{
}

/// Number of particles
void ParticleSwarm::setNumberOfParticles(const int nParticles)
{
    if (nParticles < 1)
    {
        throw std::runtime_error("Number of particles must be positive");
    }
    pImpl->mParticles = nParticles;
}

int ParticleSwarm::getNumberOfParticles() const noexcept
{
    return pImpl->mParticles;
}

/// Number of generations
void ParticleSwarm::setNumberOfGenerations(const int nGenerations)
{
    if (nGenerations < 1)
    {
        throw std::runtime_error("Number of generations must be postiive");
    }
    pImpl->mGenerations = nGenerations;
}

int ParticleSwarm::getNumberOfGenerations() const noexcept
{
    return pImpl->mGenerations;
}

/// Migrator
void ParticleSwarm::setMigrator(std::unique_ptr<IMigrator> &&migrator)
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

std::unique_ptr<IMigrator> ParticleSwarm::releaseMigrator()
{
    if (!haveMigrator()){throw std::runtime_error("Migrator not set");}
    auto result = std::move(pImpl->mMigrator);
    pImpl->mMigrator = nullptr;
    return result;
}

bool ParticleSwarm::haveMigrator() const noexcept
{
    return pImpl->mMigrator != nullptr;
}

/// The depth at which to perform the migration
void ParticleSwarm::setDepth(const double depth)
{
    if (depth < -8500 || depth > 800000)
    {
        throw std::invalid_argument("Depth must be in range [-8500, 800,000]");
    }
    pImpl->mDepth = depth;
}

/// The depth
double ParticleSwarm::getDepth() const noexcept
{
    return pImpl->mDepth;
}

/// Sets the extent
void ParticleSwarm::setExtentInX(
    const std::pair<double, double> &extent)
{
    if (extent.second <= extent.first)
    {
        throw std::invalid_argument("extent.second <= extent.first in x");
    }
    pImpl->mExtentInX = extent;
    pImpl->mHaveExtentInX = true;
}


void ParticleSwarm::setExtentInY(
    const std::pair<double, double> &extent)
{
    if (extent.second <= extent.first)
    {
        throw std::invalid_argument("extent.second <= extent.first in y");
    }
    pImpl->mExtentInY = extent;
    pImpl->mHaveExtentInY = true;
}

std::pair<double, double> ParticleSwarm::getExtentInX() const
{
    if (!pImpl->mHaveExtentInX)
    {
        throw std::runtime_error("Extent in x not set");
    }
    return pImpl->mExtentInX;
}

std::pair<double, double> ParticleSwarm::getExtentInY() const
{
    if (!pImpl->mHaveExtentInY)
    {
        throw std::runtime_error("Extent in y not set");
    }
    return pImpl->mExtentInY;
}

/// Sets the arrivals
void ParticleSwarm::setArrivals(const std::vector<Arrival> &arrivals)
{
    if (!haveMigrator()){throw std::runtime_error("Arrivals not set");}
    pImpl->mMigrator->setArrivals(arrivals);
}

bool ParticleSwarm::haveArrivals() const noexcept
{
    if (!haveMigrator()){return false;}
    return pImpl->mMigrator->haveArrivals();
}

/// Performs the optimization
void ParticleSwarm::optimize()
{
    pImpl->mHaveEvent = false;
    if (!haveMigrator())
    {
        throw std::runtime_error("Migration engine not set");
    }
    if (!haveArrivals())
    {
        throw std::runtime_error("Arrivals not set");
    }
    if (pImpl->mMigrator->getNumberOfArrivals() <
        getMinimumNumberOfArrivalsToBuildEvent() )
    {
        throw std::runtime_error(
           "Insufficient number of arrivals set to build event");
    }
    // Ensure extent is set
    auto region = pImpl->mMigrator->getGeographicRegion();
    if (!pImpl->mHaveExtentInX)
    {
        setExtentInX(region->getExtentInX());
    }   
    if (!pImpl->mHaveExtentInY)
    {   
        setExtentInY(region->getExtentInY());
    }   
    auto [x0, x1] = getExtentInX();
    auto [y0, y1] = getExtentInY();
    auto evaluationDepth = getDepth();
    // Make sure the migrator will not save the contribution history
    pImpl->mMigrator->disableSaveArrivalContributionList();
    // Instantiate the fitness class 
    ::FixedDepthDoubleDifferenceMigration fitness;
    fitness.mMigrationFunction = pImpl->mFixedDepthMigrationFunction;
    // Set bounds 
    fitness.setSearchBoundaries(std::vector<double> {x0, y0},
                                std::vector<double> {x1, y1} );
    // Instantiate the problem
    pagmo::problem problem{std::move(fitness)};
    // Instantiate the particle swarm algorithm
    pagmo::algorithm algorithm{pagmo::pso(getNumberOfGenerations())};
    // Instantiate the population
    auto populationSize = static_cast<size_t> (getNumberOfParticles());
    pagmo::population population{problem, populationSize};
    // Evolve the population
    pImpl->mLogger->debug("Beginning PSO for 2D");
    auto newPopulation = algorithm.evolve(population);
    pImpl->mLogger->debug("PSO finished!");
/*
    const auto &xHistory = newPopulation.get_x();
    const auto &fHistory = newPopulation.get_f();
    for (size_t i = 0; i < newPopulation.size(); ++i)
    {
        auto [latitude, longitude]
            = region->localToGeographicCoordinates(xHistory[i][0],
                                                   xHistory[i][1]);
        ParticleSwarmImpl::Point2D point2d{static_cast<float> (longitude),
                                           static_cast<float> (latitude),
                                           static_cast<float> (-fHistory[i][0])};
std::cout << point2d.longitude << " " << point2d.latitude << " " << point2d.value << std::endl;
    }
*/
    // Pick a winner and extract the hypocenter and origin time
    auto optimumLocation = newPopulation.champion_x();
    auto xOptimum = optimumLocation.at(0);
    auto yOptimum = optimumLocation.at(1);
    auto [latitude, longitude]
         = region->localToGeographicCoordinates(xOptimum, yOptimum);
    Event event;
    event.setLatitude(latitude);
    event.setLongitude(longitude);
    event.setDepth(evaluationDepth);
    // Compute the corresponding origin time
    try
    {
        pImpl->mMigrator->enableSaveArrivalContributionList();
        pImpl->mOptimalMigrationValue
            = pImpl->mMigrator->evaluate(xOptimum, yOptimum, evaluationDepth);
        auto contributingArrivals = pImpl->mMigrator->getContributingArrivals();
        if (!contributingArrivals.empty())
        {
            // Try to build the event
            for (const auto &arrival : contributingArrivals)
            {
                constexpr bool overWriteIfExists{false};
                if (event.canAddArrival(arrival, overWriteIfExists) >= 0)
                {
                    event.addArrival(arrival);
                }
            }
            const auto &eventArrivals = event.getArrivalsReference();
            // Get the origin time
            auto nEventArrivals = static_cast<int> (eventArrivals.size());
            if (nEventArrivals >= getMinimumNumberOfArrivalsToBuildEvent())
            {
                pImpl->mMigrator->disableSaveArrivalContributionList();
                std::vector<double> arrivalTimes(nEventArrivals);
                std::vector<double> travelTimes(nEventArrivals);
                std::vector<double> weights(nEventArrivals);
                for (int i = 0; i < nEventArrivals; ++i)
                {
                    arrivalTimes[i] = eventArrivals[i].getTime().count()*1.e-6;
                    travelTimes[i] = eventArrivals[i].getTravelTime();
                    weights[i] = eventArrivals[i].getWeight();
                }
                ULocator::Optimizers::OriginTime originTime;
                originTime.setNorm(ULocator::Optimizers::IOptimizer::Norm::L1,
                                   1);
                originTime.enableTimeReduction();
                originTime.setArrivalTimes(arrivalTimes, weights);
                originTime.setTravelTimes(travelTimes);
                originTime.optimize();
                event.setOriginTime(originTime.getTime());

                pImpl->mEvent = std::move(event);
                pImpl->mHaveEvent = true;
            }
        }
    }
    catch (const std::exception &e)
    {
        pImpl->mLogger->warn("PSO could not build event: "
                           + std::string {e.what()});
    }
}

/// Get the event
Event ParticleSwarm::getEvent() const
{
    if (!haveEvent()){throw std::runtime_error("Event not yet built");}
    return pImpl->mEvent;
}

/// Have event?
bool ParticleSwarm::haveEvent() const noexcept
{
    return pImpl->mHaveEvent;
}

/// Destructor
ParticleSwarm::~ParticleSwarm() = default; 

/// Number of arrivals to build an event
void ParticleSwarm::setMinimumNumberOfArrivalsToBuildEvent(int n)
{
    if (n < 1)
    {
        throw std::invalid_argument(
           "At least one arrival required to build an event");
    }
    pImpl->mMinimumNumberOfArrivalsToBuildEvent = n;
}

int ParticleSwarm::getMinimumNumberOfArrivalsToBuildEvent() const noexcept
{
    return pImpl->mMinimumNumberOfArrivalsToBuildEvent;
}

