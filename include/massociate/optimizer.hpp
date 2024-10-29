#ifndef MASSOCIATE_OPTIMIZER_HPP
#define MASSOCIATE_OPTIMIZER_HPP
#include <vector>
#include <memory>
#include <tuple>
#include <umps/logging/log.hpp>
namespace MAssociate
{
 class Arrival;
 class Event;
 class IMigrator;
}
namespace MAssociate
{
/// @brief An optimizer is any utility that can optimize a migration image
///        and build a corresponding event from that optimum.
class IOptimizer
{
public:
    struct OptimalLocation
    {
        double latitude{0};
        double longitude{0};
        double depth{0};
    };
public:
    IOptimizer();

    virtual void setMigrator(std::unique_ptr<IMigrator> &&migrator);
    virtual std::unique_ptr<IMigrator> releaseMigrator();
    virtual IMigrator *getMigratorHandle();
    virtual bool haveMigrator() const noexcept;

    /// @brief Sets the arrivals.
    /// @throws std::invalid_argument if \c haveMigrator() is false.
    virtual void setArrivals(const std::vector<Arrival> &arrivals);
    /// @result True indicates the arrivals were set.
    [[nodiscard]] virtual bool haveArrivals() const noexcept;


    virtual void optimize() = 0;
    //[[nodiscard]] virtual Event getEvent() const = 0;
    //[[nodiscard]] virtual bool haveEvent() const noexcept = 0;

    [[nodiscard]] virtual double getOptimalValue() const = 0;
    [[nodiscard]] virtual bool haveOptimum() const noexcept = 0;
    /// @result The optimal (latitude, longitude, depth) where the
    ///         migration image is maximal.
    [[nodiscard]] virtual std::tuple<double, double, double> getOptimalHypocenter() const = 0;
    /// @result The optimal origin time where the standard migration is maximal.
    [[nodiscard]] virtual double getOptimalOriginTime() const = 0;
    /// @result The arrivals that contributed to the optimum.
    [[nodiscard]] virtual std::vector<std::pair<Arrival, double>> getContributingArrivals() const = 0;

    virtual ~IOptimizer();
private:
    class IOptimizerImpl;
    std::unique_ptr<IOptimizerImpl> pImpl;
};
}
#endif
