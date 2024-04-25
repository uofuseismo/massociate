#ifndef MASSOCIATE_MIGRATOR_HPP
#define MASSOCIATE_MIGRATOR_HPP
#include <memory>
#include <vector>
#include <uLocator/position/knownLocalLocation.hpp>
namespace UMPS::Logging
{
 class ILog;
}
namespace ULocator
{
 class TravelTimeCalculatorMap;
 namespace Position 
 {
  class IGeographicRegion;
 }
}
namespace MAssociate
{
 class Arrival;
}
namespace MAssociate
{
/// @class IMigrator "migrator.hpp" "massociate/migrator.hpp"
/// @brief This is a base class on which we can build utilities to optimize
///        the migration image - i.e., find a set of `brightest' points at
///        a given depth slice.
/// @copyright Ben Baker (UUSS) distributed under teh MIT license.
class IMigrator
{
public:
    /// @brief This defines the optimization signal type.  For example, in the 
    ///        double-difference approach we migrate cross-correlated pick
    ///        signals whereas in the absolute approach we migrate pick signals
    ///        but increase the problem dimensionality to include a temporal
    ///        component.
    enum class SignalType
    {
        DoubleDifference, /*!< This migrates cross-correlated signals which annihalates the origin time dependence. */
        Absolute /*!< This migrates pick signals but adds a temporal dimension to the optimization. */
    };
    /// @brief Defines the pick signal to migrate.  Note, the cross-correlation
    ///        of these signls is performed anlytically.
    enum class PickSignal
    {
        Boxcar, /*!< A boxcar probability density function. */
        Gaussian, /*!< A Gaussian probability density function.  Note, this does not have compact support. */
        TruncatedGaussian /*!< A brief a truncated Gaussian probability density function. */
    };
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    IMigrator();
    /// @brief Constructor with a logger.
    explicit IMigrator(std::shared_ptr<UMPS::Logging::ILog> &logger);
    /// @brief Move constructor
    IMigrator(IMigrator &&migrator) noexcept;
    /// @}

    /// @name Region
    /// @{
 
    /// @param[in] region  This defines the search region as a simple rectangle 
    ///                    and the mapping from this simple rectangle to and
    ///                    from geographic coordinates.
    virtual void setGeographicRegion(const ULocator::Position::IGeographicRegion &region);
    /// @result The geographic region.
    /// @throws std::runtime_error \c haveGeographicRegion() is false.
    [[nodiscard]] virtual std::unique_ptr<ULocator::Position::IGeographicRegion> getGeographicRegion() const;
    /// @result True indicates the region was set.
    [[nodiscard]] virtual bool haveGeographicRegion() const noexcept;
    /// @}

    /// @name Travel Time Calculator
    /// @{

    /// @brief Sets the travel time calculator map.
    /// @param[in,out] calculatorMap  The travel time calculator map.  On exit,
    ///                               calculatorMap's behavior is undefined.
    virtual void setTravelTimeCalculatorMap(std::unique_ptr<ULocator::TravelTimeCalculatorMap> &&calculatorMap);
    /// @result True indicates the travel time calculator map was set.
    virtual bool haveTravelTimeCalculatorMap() const noexcept; 
    /// @result A pointer to the travel time calculator map.
    [[nodiscard]] virtual const ULocator::TravelTimeCalculatorMap *getTravelTimeCalculatorMap() const;
    /// @result Releases the travel time calculator map.
    [[nodiscard]] std::unique_ptr<ULocator::TravelTimeCalculatorMap> releaseTravelTimeCalculatorMap();
    /// @param[in] network  The network code - e.g., UU.
    /// @param[in] station  The station name - e.g., KNB.
    /// @param[in] phase    The phase identifier - e.g., P.
    /// @result True indicates the travel time calculator for this
    ///         network/station/phase exists in the map.
    [[nodiscard]] virtual bool haveTravelTimeCalculator(const std::string &network,
                                                        const std::string &station,
                                                        const std::string &phase) const noexcept;
    /// @}

    void setDefaultSearchLocations(const std::vector<std::unique_ptr<ULocator::Position::IKnownLocalLocation>> &locations);
    [[nodiscard]] std::unique_ptr<ULocator::Position::IKnownLocalLocation> getKnownSearchLocation(size_t index) const;

    /// @brief Sets the type of pick signal (absolute or differential)
    ///        to migrate.
    /// @param[in] signalType  The pick signal to migrate.
    void setSignalType(SignalType signalType) noexcept;
    /// @result The pick signal to migrate - i.e., absolute or differential.
    [[nodiscard]] SignalType getSignalType() const noexcept;

    /// @brief This will enable the saving of the arrival contribution list.
    void enableSaveArrivalContributionList() noexcept;
    /// @brief This will disable the saving of the arrival contribution list. 
    void disableSaveArrivalContributionList() noexcept;
    /// @result True indicates the pick contribution list should be tabulated.
    ///         By default this is false.
    [[nodiscard]] bool saveArrivalContributionList() const noexcept;

    /// @result True indicates the search history will be saved.
    [[nodiscard]] bool saveHistory() const noexcept;

    /// @brief Defines how we will define the pick signal.
    /// @param[in] signal  The pick signal.
    /// @note While a pick is technically a Dirac-delta function, it would be
    ///       highly unlikely that that a migration point produces a travel
    ///       time that exactly predicts this observed time.
    void setPickSignalToMigrate(PickSignal signal) noexcept;
    /// @result The pick signal to migrate.
    [[nodiscard]] PickSignal getPickSignalToMigrate() const noexcept;

    /// @brief Defines the maximum epicentral distance for a station to be
    ///        able to contribute to the migration.
    /// @param[in] distance   The maximum distance in meters.
    /// @throws std::invalid_argument if the distance is not positive.
    void setMaximumEpicentralDistance(double distance);
    /// @result The maximum epicentral distance in meters for a station
    ///         to be able to to contribute to the migration.
    [[nodiscard]] double getMaximumEpicentralDistance() const noexcept;

    /// @brief Sets the arrivals to migrate.
    void setArrivals(const std::vector<Arrival> &arrivals);
    /// @result The arrivals to migrate.
    [[nodiscard]] std::vector<Arrival> getArrivals() const;
    /// @result The number of arrivals.
    [[nodiscard]] int getNumberOfArrivals() const noexcept;
    /// @result The a constant reference to the arrivals to migrate.
    [[nodiscard]] const std::vector<Arrival> &getArrivalsReference() const;
    /// @result True indicates the arrivals were not set.
    [[nodiscard]] bool haveArrivals() const noexcept;

    /// @brief Evaluates at known locations.
    [[nodiscard]] virtual std::vector<double> evaluateAtKnownLocations() const;
    /// @throws std::runtime_error if \c haveTravelTimeCalculatorMap() is false
    ///         or \c haveArrivals() is false is false.
    [[nodiscard]] virtual double evaluate(double x, double y, double z) const;

    /// @result The list of arrivals that contributed to the stack.
    /// @note If the signal to migrate has global support then this will be
    ///       true for all arrivals.
    /// @throws std::runtime_error if evaluate was not performed with
    ///         \c saveArrivalContributionList was false.
    [[nodiscard]] std::vector<std::pair<Arrival, double>> getContributingArrivals() const;

    /// @brief Performs the optimization over the region.
    //virtual void optimize() = 0;

    /// @name Destructors
    /// @{

    /// @brief Destructor.
    virtual ~IMigrator();
    /// @}

    IMigrator& operator=(IMigrator &&migrator) noexcept;

    IMigrator(const IMigrator &) = delete;
    IMigrator& operator=(const IMigrator &) = delete;
private:
    class IMigratorImpl;
    std::unique_ptr<IMigratorImpl> pImpl;
};
}
#endif
