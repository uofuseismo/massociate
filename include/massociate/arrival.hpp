#ifndef MASSOCIATE_ARRIVAL_HPP
#define MASSOCIATE_ARRIVAL_HPP
#include <memory>
#include <chrono>
//#include "massociate/enums.hpp"
namespace MAssociate
{
 class Pick;
 class WaveformIdentifier;
}
namespace MAssociate
{
/// @class Arrival "arrival.hpp" "massociate/arrival.hpp"
/// @brief Defines an arrival.  This is an unassociated pick.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Arrival
{
public:
    /// @brief Some typical seismic phases.
    enum class Phase
    {
        P, /*!< A first-arrival compressional arrival. */
        S  /*!< A shear-wave arrival. */
    };
    /// @brief Defines the first motion of an arrival.
    /// @copyright Ben Baker (University of Utah) distributed under the MIT license.
    enum class FirstMotion
    {
        Up = 1,      /*!< Upwards (compressional) first motion. */
        Unknown = 0, /*!< Unknown first motion. */
        Down  =-1    /*!< Downwards (dilitational) first motion. */
    };
public:
    /// @name Constructors
    /// @{
    /// @brief Constructor
    Arrival();
    /// @brief Copy constructor.
    /// @param[in] arrival  The class from which to initialize this class.
    Arrival(const Arrival &arrival);
    /// @brief Creats an arrival from a pick.
    /// @param[in] arrival  The pick from which to create this arrival.
    explicit Arrival(const Pick &pick);
    /// @brief Move constructor.
    /// @param[in,out] arrival  The class from which to initialize this class.
    ///                         On exit, arrival's behavior is undefined.
    Arrival(Arrival &&arrival) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy assignment operator.
    /// @param[in] arrival  The arrival class to copy to this.
    /// @result A deep copy of the arrival.
    Arrival& operator=(const Arrival &arrival);
    /// @brief Creates an arrival from a pick.
    /// @param[in] pick   The pick to convert to an arrival.
    /// @result The resulting arrival.
    /// @note You may still have to add additional fields to the arrival.
    Arrival& operator=(const Pick &pick);
    /// @brief Move assignment operator.
    /// @param[in,out] arrival  The arrival class whose memory will be moved
    ///                         to this.  On exit, arrival's behavior is
    ///                         undefined.
    /// @result The memory from arrival moved to this.
    Arrival& operator=(Arrival &&arrival) noexcept;
    /// @} 

    /// @name Destructors
    /// @{
    ~Arrival();
    /// @brief Resets the class.
    void clear() noexcept;
    /// @}

    /// @name Waveform Identifier
    /// @{
    /// @brief Sets the station, network, channel, and location code on which
    ///        the pick was observed.
    /// @param[in] waveformIdentifier  The waveform identifier.  This
    ///                                should at least contain the station 
    ///                                and network but cannot be empty.
    /// @throws std::invalid_argument if \c waveformIdentifier.isEmpty()
    ///         is true.
    void setWaveformIdentifier(
        const MAssociate::WaveformIdentifier &waveformIdentifier);
    /// @result Gets the waveform identifier on which the arrival was observed.
    /// @throws std::runtime_error if the waveform identifier was not set.
    /// @sa \c haveWaveformIdentifier()
    [[nodiscard]] MAssociate::WaveformIdentifier getWaveformIdentifier() const;
    /// @result True indicates that the waveform identifier was set.
    [[nodiscard]] bool haveWaveformIdentifier() const noexcept;
    /// @}
 
    /// @name Phase Name
    /// @{

    /// @brief Sets the phase name.
    /// @param[in] phase  The name of the seismic phase.
    void setPhase(const Arrival::Phase phase) noexcept;
    /// @brief Sets the phase name.
    /// @param[in] phase  The name of the seismic phase.  
    /// @throws std::invalid_argument if the phaseName is empty.
    void setPhase(const std::string &phaseName);
    /// @brief Gets the name of the seismic phase for this arrival.
    /// @result The name of the seismic phase.
    /// @throws std::invalid_argument if the phase name is not set.
    /// @sa \c havePhase()
    [[nodiscard]] std::string getPhase() const;
    /// @result True indicates that the phase name is set.
    [[nodiscard]] bool havePhase() const noexcept;
    /// @}

    /// @name Arrival Time
    /// @{

    /// @brief Sets the arrival time.
    /// @param[in] time  The arrival time (UTC) in micro-seconds since
    ///                  the epoch.
    void setTime(const std::chrono::microseconds &arrivalTime) noexcept;
    /// @brief Sets the arrival time.
    /// @param[in] time  The arrival time (UTC) in seconds since the epoch.
    void setTime(double time) noexcept;
    /// @brief Gets the arrival time.
    /// @result The arrival time (UTC) in microseconds since the epoch.
    /// @throws std::runtime_error if the arrival time was not set.
    /// @sa \c haveTime()
    [[nodiscard]] std::chrono::microseconds getTime() const;
    /// @result True indicates that the arrival time was set.
    [[nodiscard]] bool haveTime() const noexcept;
    /// @}

    /// @name Identifier
    /// @{

    /// @brief Sets a unique arrival identifier number.  This likely will
    ///        be assigned by the real-time system but an algorithm as simple
    ///        as counting picks can be used.
    /// @param[in] identifier   The unique arrival identifier.
    void setIdentifier(uint64_t identifier) noexcept;
    /// @result The arrival's unique identification number.
    /// @throws std::runtime_error if \c haveIdentifier() is false.
    [[nodiscard]] uint64_t getIdentifier() const;
    /// @result True indicates that the identifier has been set.
    [[nodiscard]] bool haveIdentifier() const noexcept;
    /// @}

    /// @name Weight
    /// @{

    /// @brief Defines the arrival's standard error (deviation) in seconds.
    /// @param[in] standardError   The standard deviation in seconds.
    /// @throws std::invalid_argument if the standard deviation is not positive.
    /// @note When migrating with a Gaussian distribution this is the canonical
    ///       standard deviation.  When migrating with a boxcar function the
    ///       standard deviation is \f$ \sigma = \sqrt{ \frac{(B-A)^2}{12} } \f$
    ///       where $W = B - A$ is the boxcar's width.  I find it convenient to
    ///       think of parameterizing boxcars based on width, so, in this case I
    ///       use \f$ \sigma = \frac{W}{2 \sqrt{3}} \f$. 
    void setStandardError(double standardError);
    /// @result The arrival's standard error in seconds.
    [[nodiscard]] double getStandardError() const noexcept;
    /// @result The arrival's weight in 1/seconds.  This is computed from
    ///         the standard error i.e., \f$ \frac{1}{\sigma} \f$.
    [[nodiscard]] double getWeight() const noexcept;
    /// @}

    /// @name First Motion
    /// @{

    /// @brief Sets the first motion.
    /// @param[in] firstMotion  The arrival's first motion.
    void setFirstMotion(Arrival::FirstMotion firstMotion) noexcept;
    /// @brief Gets the arrival's first motion.
    [[nodiscard]] MAssociate::Arrival::FirstMotion getFirstMotion() const noexcept;
    /// @brief Sets the first motion's weight.  This is useful when the first
    ///        motion is produced by an ML algorithm and has an accompanying
    ///        posterior probability.
    /// @param[in] weight  The first motion's weight.
    /// @throws std::invalid_argument if the weight is not positive.
    void setFirstMotionWeight(double weight);
    /// @result The first motoin's weight.
    [[nodiscard]] double getFirstMotionWeight() const noexcept;
    /// @}

    /// @name Travel Time
    /// @{

    /// @brief Sets the travel time.
    /// @param[in] travelTime  The source to receiver travel time in seconds.
    /// @throws std::invalid_argument if this is negative
    void setTravelTime(double travelTime);
    /// @result The travel time in seconds from the source to the receiver.
    /// @throws std::runtime_error if \c haveTravelTime() is false.
    [[nodiscard]] double getTravelTime() const;
    /// @result True indicates that the travel time was set.
    [[nodiscard]] bool haveTravelTime() const noexcept;

private:
    class ArrivalImpl;
    std::unique_ptr<ArrivalImpl> pImpl;
};
}
#endif
