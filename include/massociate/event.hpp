#ifndef MASSOCIATE_EVENT_HPP
#define MASSOCIATE_EVENT_HPP
#include <memory>
namespace MAssociate
{
class Arrival;
/*!
 * @brief Defines the preliminary event information.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class Event
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor
     */
    Event();
    /*!
     * @brief Copy constructor.
     * @param[in] event   The event class from which to initialize this class.
     */
    Event(const Event &event);
    /*!
     * @brief Move constructor.
     * @param[in,out] event  The event class from which to initialize this
     *                       class.  On exit, event's behavior is undefined.
     */
    Event(Event &&event) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment.
     * @param[in] event   The event class to copy to this.
     * @result A deep copy of the event class.
     */
    Event& operator=(const Event &event);
    /*!
     * @brief Move assignment.
     * @param[in,out] event  The event class whose memory will be moved to this.
     *                       On exit, event's behavior is undefined.
     * @result The memory from event moved to this.
     */
    Event& operator=(Event &&event) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor
     */
    ~Event();
    /*!
     * @brief Resets the class and releases memory.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Event Identifier
     * @{
     */
    /*!
     * @brief Sets the event identifier.
     * @param[in] evid   The event identifier.
     */
    void setIdentifier(uint64_t evid) noexcept;
    /*!
     * @result The event identifier.
     * @throws std::runtime_error if the \c haveIdentifier() is false.
     */
    [[nodiscard]] uint64_t getIdentifier() const;
    /*!
     * @result True indicates that the event identifier is set.
     */
    [[nodiscard]] bool haveIdentifier() const noexcept;
    /*! @} */

    /*! @name Hypocenter and Origin Time
     * @{
     */
    /*!
     * @brief Sets the event's x position.  This could also be the longitude.
     * @param[in] x   The source's x position in meters (or degrees).
     */
    void setXPosition(double x) noexcept;
    /*!
     * @result The event's x position.
     */ 
    [[nodiscard]] double getXPosition() const;
    /*!
     * @result True indicates that the x position has been set.
     */
    [[nodiscard]] bool haveXPosition() const noexcept;

    /*!
     * @brief Sets the event's y position.  This could also be the latitude.
     * @param[in] x   The source's y position in meters (or degrees).
     */
    void setYPosition(double y) noexcept;
    /*!
     * @result The event's y position.
     */
    [[nodiscard]] double getYPosition() const;
    /*!
     * @result True indicates that the y position has been set.
     */
    [[nodiscard]] bool haveYPosition() const noexcept;

    /*!
     * @brief Sets the event's z position.
     * @param[in] z   The source's z position in meters.
     */
    void setZPosition(double x) noexcept;
    /*!
     * @result The event's z position.
     */
    [[nodiscard]] double getZPosition() const;
    /*!
     * @result True indicates that the z position has been set.
     */
    [[nodiscard]] bool haveZPosition() const noexcept;

    /*!
     * @brief Sets the event's latitude.
     * @param[in] latitude   The event's latitude in degrees  where latitude
     *                       increases positive north.
     * @throws std::invalid_argument if the latitude is not in the
     *         range [-90,90].
     */
    void setLatitude(double latitude);
    /*!
     * @result The event's latitude in degrees.
     * @throws std::runtime_error if \c haveLatitude() is false.
     */
    [[nodiscard]] double getLatitude() const;
    /*!
     * @result True indicates that the event's latitude is set.
     */
    [[nodiscard]] bool haveLatitude() const noexcept;

    /*!
     * @brief Sets the event's longitude.
     * @param[in] longitude  The event's longitude in degrees where longitude
     *                       increases positive east.
     * @throws std::invalid_argument if the longitude is not in the range
     *         [-540,540).
     */
    void setLongitude(double longitude);
    /*!
     * @result The event's longitude.
     * @note This will be in the range [0,360).
     * @throws std::runtime_error if \c haveLongitude() is false.
     */
    [[nodiscard]] double getLongitude() const;
    /*!
     * @brief True indicates that the longitude is set.
     */
    [[nodiscard]] bool haveLongitude() const noexcept;

    /*!
     * @brief The event depth in meters.
     * @param[in] depth   The event depth.  This increases positive down from
     *                    some reference such as sea-level..
     */
    void setDepth(double depth) noexcept;
    /*!
     * @result The event's depth in meters.
     * @throws std::runtime_error of \c haveDepth() is false.
     */
    [[nodiscard]] double getDepth() const;
    /*!
     * @result True indicates that the event depth is set.
     */
    [[nodiscard]] bool haveDepth() const noexcept;
   
    /*!
     * @brief Sets the event's origin time.
     * @param[in] originTime  The UTC origin time in seconds since the epoch.
     */
    void setOriginTime(double originTime) noexcept; 
    /*!
     * @result The origin time.
     * @throws std::runtime_error of \c haveOriginTime() is false.
     */
    [[nodiscard]] double getOriginTime() const;
    /*!
     * @result True indicates that the origin time is set.
     */
    [[nodiscard]] bool haveOriginTime() const noexcept;
    /*! @} */

    /*! @name Arrivals
     * @{
     */
    /*!
     * @brief Adds an arrival to the event.
     * @param[in] arrival  The arrival to add to the event.  This must have
     *                     a waveform identifier, a phase name, an arrival time,
     *                     and an arrival identifier.
     * @note If the station/phase pair already exist then the old arrival
     *       will be overwritten. 
     */
    void addArrival(const Arrival &arrival);
    /*!
     * @result A vector of all arrivals set for this event.
     * @note To sort in increasing order of arrival time try:
     *       std::sort(arrivals.begin(), arrivals.end(),
     *                 [](const Arrival &a, const Arrival &b)
     *                 {
     *                    return a.getTime() < b.getTime();
     *                 });
     */
    [[nodiscard]] std::vector<Arrival> getArrivals() const;
    /*!
     * @brief Clears all the arrivals from the class.
     */
    void clearArrivals() noexcept;
    /*!
     * @result The number of arrivals corresponding to this event.
     */
    [[nodiscard]] int getNumberOfArrivals() const noexcept;
    /*! @} */
private:
    class EventImpl;
    std::unique_ptr<EventImpl> pImpl;
};
}
#endif
