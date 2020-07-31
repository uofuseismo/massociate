#ifndef MASSOCIATOR_ASSOCIATOR_HPP
#define MASSOCIATOR_ASSOCIATOR_HPP
#include <memory>
namespace MAssociate
{
class Pick;
class Event;
class Arrival;
class AssociatorParameters;
namespace Mesh
{
template<class T> class IMesh;
}
};
namespace MAssociate
{
/*!
 * @brief This class performs the association.
 */
template<class T>
class Associator
{
public:
    /*!
     * @brief Constructor.
     */
    Associator();

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~Associator();
    /*!
     * @brief Releases all memory and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Step 1: Initiliazation
     * @{
     */
    /*!
     * @brief Initializes the associator class.
     * @param parameters   The parameters from which to initialize this class.
     * @throws std::invalid_argument if required parameters are not set.
     */
    void initialize(const AssociatorParameters &parameters,
                    const Mesh::IMesh<T> &geometry);
    /*!
     * @return True indicates that the class is initialized.
     */
    [[nodiscard]] bool isInitialized() const noexcept;
    /*! @} */

    /*! @name Step 2: Set Travel Time Tables
     * @{
     */
    /*!
     * @return The expected number of points in a travel time table.
     * @throws std::runtime_error if the class is not initialized.
     */
    [[nodiscard]] int getNumberOfPointsInTravelTimeTable() const;
    /*!
     * @brief Sets the travel time table for the network/station/phase.
     * @param network   The network name.
     * @param station   The station name.
     * @param phase     The seismic phase, e.g., P or S.
     * @param nPoints   The number of points in the field.  This must match
     *                  \c getNumberOfPointsInTravelTimeTable().
     * @param times
     */
    template<typename U>
    void setTravelTimeTable(const std::string &network,
                            const std::string &station,
                            const std::string &phase,
                            int nPoints, const U times[]);
    /*!
     * @brief This is a debugging routine for getting the travel time table
     *        corresponding to the network, station, phase.
     * @return The travel times (seconds) from the given network/station/phase
     *         to all points.
     * @throws std::runtime_error if the travel time table doesn't exist.
     * @sa \c haveTravelTimeTable()
     */
    [[nodiscard]]
    std::vector<T> getTravelTimeTable(const std::string &network,
                                      const std::string &station,
                                      const std::string &phase) const;
    /*!
     * @param[in] network  The network name.
     * @param[in] station  The station name.
     * @param[in] phase    The seismic phase.
     * @return True indicates the table for this network/station/phase tuple
     *         exists.
     */
    [[nodiscard]]
    bool haveTravelTimeTable(const std::string &network,
                             const std::string &station,
                             const std::string &phase) const noexcept;
    /*!
     * @return True indicates that all travel time tables have been set.
     */
    [[nodiscard]] bool haveAllTravelTimeTables() const noexcept;
    /*!
     * @brief Gets the travel time for the network/station/phase given a
     *        source at the specified index.
     * @param[in] network   The network name.
     * @param[in] station   The station name.
     * @param[in] phase     The phase name.
     * @param[in] index     The index.  This must be in the range
     *                      [0, \c getNumberOfPointsInTravelTimeTable()].
     * @result The travel time of a phase from the `source' at the index
     *         to the given station in seconds.
     */
    [[nodiscard]] T getTravelTime(const std::string &network,
                                  const std::string &station,
                                  const std::string &phase,
                                  int index) const;
    /*!
     * @result The maximum differential travel time in seconds.
     */
    [[nodiscard]] T getMaximumDifferentialTravelTime() const noexcept;
    /*! @} */

    /*! @name Step 3: Set Picks
     * @{
     */
    /*!
     * @brief Adds an pick whose differential travel time will be migrated.
     * @param[in] pick   The pick whose time will be migrated.
     * @throws std::invalid_argument if the pick's network/station/phase
     *         do not correspond to an existing travel time table or the
     *         pick time is not set.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c haveTravelTimeTable().
     */
    void addPick(const Pick &pick);
    /*!
     * @return The number of picks.
     */
    [[nodiscard]] int getNumberOfPicks() const noexcept;
    /*!
     * @brief Clears the picks.
     * @note This should be used when seeking to migrate a new batch of
     *       picks using the existing travel time tables.
     */
    void clearPicks() noexcept;
    /*! @} */

    /*! @name Step 4: Create Events
     * @{
     */
    void associate();
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
    void bindPickToEvent(const MAssociate::Pick &pick);
    /*! @} */

    /*! @name Step 5: Get Events
     * @{
     */
    /*!
     * @result A vector of events.
     */
    [[nodiscard]] std::vector<MAssociate::Event> getEvents() const;
    /*!
     * @brief Clears the events.  This is useful for when you want to release
     *        memory on old events.
     * @param[in] resetCounter  If true then the event ID counter will be reset.
     *                          This behavior is not typically recommended. 
     * @note Call \c getEvents() prior to calling this otherwise you will lose
     *       your events.
     */
    void clearEvents(bool resetCounter = false) noexcept;
    /*!
     * @result The number of events in the associator.
     */
    [[nodiscard]] int getNumberOfEvents() const noexcept;
private:
    class AssociatorImpl;
    std::unique_ptr<AssociatorImpl> pImpl;
};
}
#endif
