#ifndef MASSOCIATE_PARTICLE_SWARM_HPP
#define MASSOCIATE_PARTICLE_SWARM_HPP
#include <memory>
namespace MAssociate
{
 class IMigrator;
 class Arrival;
 class Event;
}
namespace MAssociate
{
class GoldenSectionSearchOptions
{
public:
 
private:

};
class ParticleSwarm
{
public:
    /// @brief Constructor.
    ParticleSwarm();

    /// @brief Sets the migration engine.
    void setMigrator(std::unique_ptr<IMigrator> &&migrator);
    /// @result Releases the migration engine.
    [[nodiscard]] std::unique_ptr<IMigrator> releaseMigrator(); 
    /// @result True indicates the migration engine was set.
    [[nodiscard]] bool haveMigrator() const noexcept;

    /// @brief Sets the arrivals.
    /// @throws std::invalid_argument if \c haveMigrator() is false.
    void setArrivals(const std::vector<Arrival> &arrivals);
    /// @result True indicates the arrivals were set.
    [[nodiscard]] bool haveArrivals() const noexcept;

    /// @brief To build an event, this number of arrivals is required to
    ///        contribute to the image maximum.
    void setMinimumNumberOfArrivalsToBuildEvent(int nArrivals);
    /// @result After migrating at at a point, this number of arrivals is
    ///         required to contribute in order to build the event.
    [[nodiscard]] int getMinimumNumberOfArrivalsToBuildEvent() const noexcept;

    /// @brief Sets the number of particles.
    /// @throws std::invalid_argument if nParticles is not positive.
    void setNumberOfParticles(int particles);
    /// @result The number of particles.
    [[nodiscard]] int getNumberOfParticles() const noexcept;
    /// @brief Sets the number of generations to evolve the particle
    ///        population.
    /// @throws std::invalid_argument if nGenerations is not positive.
    void setNumberOfGenerations(int nGenerations);
    /// @result The number of generations.
    [[nodiscard]] int getNumberOfGenerations() const noexcept;

    /// @brief Sets the depth at which the migration optimization will
    ///        be performed.
    /// @param[in] depth   The depth in meters.
    void setDepth(double depth);
    /// @result The depth at which the migration will be performed.
    [[nodiscard]] double getDepth() const noexcept;

    /// @brief Sets the search extent in X.
    /// @param[in] extentInX   The lower and upper extent to search in x in
    ///                        meters.  This should be a subset of the region's
    ///                        extent.
    void setExtentInX(const std::pair<double, double> &extentInX);
    /// @result The extent to search in x.  By default this is the region's
    ///         extent.
    [[nodiscard]] std::pair<double, double> getExtentInX() const;
    /// @brief Sets the search extent in Y.
    /// @param[in] extentInY   The lower and upper extent to search in y in
    ///                        meters.  This should be a subset of the region's
    ///                        extent.
    void setExtentInY(const std::pair<double, double> &extentInY);
    /// @result The extent to search in y.  By default this is the region's
    ///         extent.
    [[nodiscard]] std::pair<double, double> getExtentInY() const;

    /// @brief Performs the particle swarm optimization.
    void optimize();

    /// @result True indicates we were able to migrate the arrivals and create
    ///         an event.
    [[nodiscard]] bool haveEvent() const noexcept;
    /// @result The event with associated arrivals.
    /// @throws std::runtime_error if \c haveEvent() is false.
    [[nodiscard]] Event getEvent() const;

    /// @brief Destructor.
    ~ParticleSwarm();

    ParticleSwarm(const ParticleSwarm &) = delete;
    ParticleSwarm& operator=(const ParticleSwarm &) = delete;
private:
    class ParticleSwarmImpl;
    std::unique_ptr<ParticleSwarmImpl> pImpl; 
};
}
#endif
