#ifndef MASSOCIATE_ARRIVAL_HPP
#define MASSOCIATE_ARRIVAL_HPP
#include <memory>
#include "massociate/enums.hpp"
namespace MAssociate
{
class Pick;
class WaveformIdentifier;
/*!
 * @class Arrival "arrival.hpp" "massociate/arrival.hpp"
 * @brief Defines an arrival.  This is an unassociated pick.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class Arrival
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor
     */
    Arrival();
    /*!
     * @brief Copy constructor.
     * @param[in] arrival  The class from which to initialize this class.
     */
    Arrival(const Arrival &arrival);
    /*!
     * @brief Creats an arrival from a pick.
     * @param[in] arrival  The pick from which to create this arrival.
     */
    explicit Arrival(const Pick &pick);
    /*!
     * @brief Move constructor.
     * @param[in,out] arrival  The class from which to initialize thi class.
     *                         On exit, arrival's behavior is undefined.
     */
    Arrival(Arrival &&arrival) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] arrival  The arrival class to copy to this.
     * @result A deep copy of the arrival.
     */
    Arrival& operator=(const Arrival &arrival);
    /*!
     * @brief Creates an arrival from a pick.
     * @param[in] pick   The pick to convert to an arrival.
     * @result The resulting arrival.
     * @note You may still have to add additional fields to the arrival.
     */
    Arrival& operator=(const Pick &pick);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] arrival  The arrival class whose memory will be moved
     *                         to this.  On exit, arrival's behavior is
     *                         undefined.
     * @result The memory from arrival moved to this.
     */
    Arrival& operator=(Arrival &&arrival) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    ~Arrival();
    /*!
     * @brief Resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Station, Network, Channel, Location
     * @{
     */
    /*!
     * @brief Sets the station, network, channel, and location code on which
     *        the pick was observed.
     * @param[in] waveformIdentifier   The waveform identifie.  This
     *                                 should at least contain the station 
     *                                 and network but cannot be empty.
     * @throws std::invalid_argument if \c waveformIdentifier.isEmpty() is true.
     */
    void setWaveformIdentifier(
        const MAssociate::WaveformIdentifier &waveformIdentifier);
    /*!
     * @result Gets the waveform identifier on which the arrival was observed.
     * @throws std::runtime_error if the waveform identifier was not set.
     * @sa \c haveWaveformIdentifier()
     */
    [[nodiscard]] MAssociate::WaveformIdentifier getWaveformIdentifier() const;
    /*!
     * @result True indicates that the waveform identifier was set.
     */
    [[nodiscard]] bool haveWaveformIdentifier() const noexcept;
    /*! @} */
 
    /*! @name Phase Name
     * @{
     */
    /*!
     * @brief Sets the phase name.
     * @param[in] phaseName  The name of the seismic phase.  
     * @throws std::invalid_argument if the phaseName is empty.
     */
    void setPhaseName(const std::string &phaseName);
    /*!
     * @brief Gets the name of the seismic phase for this arrival.
     * @result The name of the seismic phase.
     * @throws std::invalid_argument if the phase name is not set.
     * @sa \c havePhaseName()
     */
    [[nodiscard]] std::string getPhaseName() const;
    /*!
     * @result True indicates that the phase name is set.
     */
    [[nodiscard]] bool havePhaseName() const noexcept;
    /*! @} */

    /*! @name Arrival Time
     * @{
     */
    /*!
     * @brief Sets the arrival time.
     * @param[in] arrivalTime  The arrival time in UTC seconds since the epoch.
     */
    void setTime(double arrivalTime) noexcept;
    /*!
     * @brief Gets the arrival time.
     * @result The arrival time in UTC seconds since the epoch.
     * @throws std::runtime_error if the arrival time was not set.
     * @sa \c haveTime()
     */
    [[nodiscard]] double getTime() const;
    /*!
     * @result True indicates that the arrival time was set.
     */
    [[nodiscard]] bool haveTime() const noexcept;
    /*! @} */

    /*! @name Identifier
     * @{
     */
    /*!
     * @brief Sets a unique arrival identifier number.  This likely will
     *        be assigned by the real-time system but an algorithm as simple
     *        as counting picks can be used.
     * @param[in] identifier   The unique arrival identifier.
     */
    void setIdentifier(uint64_t identifier) noexcept;
    /*!
     * @return The arrival's unique identification number.
     * @throws std::runtime_error if \c haveIdentifier() is false.
     */
    [[nodiscard]] uint64_t getIdentifier() const;
    /*!
     * @return True indicates that the identifier has been set.
     */
    [[nodiscard]] bool haveIdentifier() const noexcept;
    /*! @} */

    /*! @name Weight
     * @{
     */
    /*!
     * @brief Defines the arrival's standard deviation (error) in seconds.
     * @param[in] std   The standard deviation in seconds.
     * @throws std::invalid_argument if the standard deviation is not positive.
     * @note When migrating with a Gaussian distribution this is the canonical
     *       standard deviation.  When migrating with a boxcar function the
     *       standard deviation is \f$ \sigma = \sqrt{ \frac{(B-A)^2}{12} } \f$
     *       where $W = B - A$ is the boxcar's width.  I find it convenient to
     *       think of parameterizing boxcars based on width, so, in this case I
     *       use \f$ \sigma = \frac{W}{2 \sqrt{3}} \f$. 
     */
    void setStandardDeviation(double std);
    /*!
     * @result The arrival's standard deviation in seconds.
     */
    [[nodiscard]] double getStandardDeviation() const noexcept;
    /*!
     * @result The arrival's weight in 1/seconds.  This is computed from
     *         the standard deviation i.e., \f$ \frac{1}{\sigma} \f$.
     */
    [[nodiscard]] double getWeight() const noexcept;
    /*! @} */

    /*! @name Polarity
     * @{
     */
    /*!
     * @brief Sets the polarity.
     * @param[in] polarity  The arrival's polarity.
     */
    void setPolarity(MAssociate::Polarity polarity) noexcept;
    /*!
     * @brief Gets the arrival's polarity.
     */
    [[nodiscard]] MAssociate::Polarity getPolarity() const noexcept;
    /*!
     * @brief Sets the polarities weight.  This is useful when the polarity
     *        is produced by an ML algorithm and has an accompanying posterior
     *        probability.
     * @param[in] weight  The polarity's weight.
     * @throws std::invalid_argument if the weight is not positive.
     */
    void setPolarityWeight(double weight);
    /*!
     * @result The polarity's weight.
     */
    [[nodiscard]] double getPolarityWeight() const noexcept;
    /*! @} */

    /*! @name Travel Time
     * @{
     */
    /*!
     * @param[in] travelTime  The source to receiver travel time in seconds.
     * @throws std::invalid_argument if this is negative
     */
    void setTravelTime(double travelTime);
    /*!
     * @result The travel time in seconds from the source to the receiver.
     * @throws std::runtime_error if \c haveTravelTime() is false.
     */
    [[nodiscard]] double getTravelTime() const;
    /*!
     * @result True indicates that the travel time was set.
     */
    [[nodiscard]] bool haveTravelTime() const noexcept;

    /*! @name Static Correction
     * @{
     */
    /*!
     * @brief Sets the static currection such that the modeled pick time is
     *        \f$ T_{modeled} = T_0 + T + T_s \f$
     *        where \f$ T_0 \f$ is the origin time, \f$ T \f$ is the travel time
     *        and \f$ T_s \f$ is the static correction.
     * @param[in] correction   The static correction in seconds.
     */
    void setStaticCorrection(double correction) noexcept;
    /*!
     * @result The static correction in seconds.
     */
    [[nodiscard]] double getStaticCorrection() const noexcept;
    /*! @} */
private:
    class ArrivalImpl;
    std::unique_ptr<ArrivalImpl> pImpl;
};
}
#endif
