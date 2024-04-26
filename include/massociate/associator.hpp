#ifndef MASSOCIATOR_ASSOCIATOR_HPP
#define MASSOCIATOR_ASSOCIATOR_HPP
#include <memory>
#include <chrono>
#include <vector>
#include <umps/logging/log.hpp>
namespace MAssociate
{
 class Arrival;
 class AssociatorParameters;
 class Event;
 class IClusterer;
 class IOptimizer;
 class Pick;
}
namespace MAssociate
{
/// @brief This class performs the association.
class Associator
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Associator();
    /// @brief Constructor with a given logger.
    explicit Associator(std::shared_ptr<UMPS::Logging::ILog> &logger);
    /// @}

    /// @name Step 1: Initiliazation
    /// @{
  
    /// @brief Initializes the associator class.
    /// @param parameters   The parameters from which to initialize this class.
    void initialize(int nArrivalsToNucleate,
                    int nPArrivalsToNucleate = 3); //const AssociatorParameters &parameters);
    /// @result The minimum number of arrivals required to associate.
    [[nodiscard]] int getMinimumNumberOfArrivalsToNucleate() const noexcept;    
void setMinimumNumberOfPArrivalsToNucleate(int nArrivals);
    /// @result The minimum number of P arrivals required to associate.
    [[nodiscard]] int getMinimumNumberOfPArrivalsToNucleate() const noexcept;
    /// @return True indicates that the class is initialized.
    [[nodiscard]] bool isInitialized() const noexcept;
    /// @}

    /// @name Step 2: Set the optimizer and clustering engines
    /// @{

    /// @brief Sets the optimizer.
    void setOptimizer(std::unique_ptr<IOptimizer> &&optimizer);
    /// @result The optimizer.
    [[nodiscard]] std::unique_ptr<IOptimizer> releaseOptimizer();
    /// @result True indicates the optimizer was set.
    [[nodiscard]] bool haveOptimizer() const noexcept;

    /// @brief Sets the clusterer.
    void setClusterer(std::unique_ptr<IClusterer> &&clusterer);
    /// @result The clusterer.
    [[nodiscard]] std::unique_ptr<IClusterer> releaseClusterer();
    /// @result True indicates the clusterer was set.
    [[nodiscard]] bool haveClusterer() const noexcept;
    /// @}

    /// @name Step 3: Set Picks
    /// @{

    void setPicks(const std::vector<Pick> &pick);
    void setPicks(std::vector<Pick> &&picks);

    /// @brief Adds an pick whose arrival travel time will be migrated.
    /// @param[in] pick   The pick whose time will be migrated.
    /// @throws std::invalid_argument if the pick's network/station/phase
    ///         do not correspond to an existing travel time table or the
    ///         pick time is not set.
    /// @throws std::runtime_error if the class is not initialized.
    void addPick(const Pick &pick);
    /// @return The number of picks.
    [[nodiscard]] int getNumberOfPicks() const noexcept;
    /// @brief Clears the picks.
    /// @note This should be used when seeking to migrate a new batch of
    ///       picks using the existing travel time tables.
    void clearPicks() noexcept;
    /// @}

    /// @name Step 4: Create Events
    /// @{

    /// @brief Creates events whose origin times are in the given interval.
    /// @param[in] minimumOriginTime
    ///  
    void associate(const std::chrono::seconds &minimumOriginTime = std::chrono::seconds{-2208988800},
                   const std::chrono::seconds &maximumOriginTime = std::chrono::seconds{4102444800});
    /// @}

    /*!
     * @brief This will scan through the available events created by calling
     *        \c associate() and attempt to bind a pick to an appropriate event.
     * @param[in] pick   The pick to attempt to bind to an event.
     * @throws std::invalid_argument if the pick does not have a time, phase,
     *         waveform identifier, pick ID, or corresponding travel time table.
     * @throws std::runtime_error if \c isInitialized() is false.
     * @note This is useful for when you have a very noisy station that isn't
     *       appropriate for association but still may offer useful information
     *       to the locator.
     */
    //void bindPickToEvent(const MAssociate::Pick &pick);
    /// @}

    /// @name Step 5: Get Events
    /// @{
 
    /// @result A vector of events.
    [[nodiscard]] std::vector<MAssociate::Event> getEvents() const;
    /*!
     * @brief Clears the events.  This is useful for when you want to release
     *        memory on old events.
     * @param[in] resetCounter  If true then the event ID counter will be reset.
     *                          This behavior is not typically recommended. 
     * @note Call \c getEvents() prior to calling this otherwise you will lose
     *       your events.
     */
    //void clearEvents(bool resetCounter = false) noexcept;
    /// @result The number of events in the associator.
    [[nodiscard]] int getNumberOfEvents() const noexcept;
    /// @result The unassociated picks.
    [[nodiscard]] std::vector<Pick> getUnassociatedPicks() const;
    /// @}

    /// @name Destructors
    /// @{
 
    /// @brief Releases all memory and resets the class.
    void clear() noexcept;
    /// @brief Destructor.
    ~Associator();
    /// @}

    Associator(const Associator &) = delete;
    Associator& operator=(const Associator &) = delete;
private:
    class AssociatorImpl;
    std::unique_ptr<AssociatorImpl> pImpl;
};
}
#endif
