#ifndef MASSOCIATE_PICK_HPP
#define MASSOCIATE_PICK_HPP
#include <memory>
#include "massociate/enums.hpp"
namespace MAssociate
{
class Arrival;
class WaveformIdentifier;
/// @class Pick "pick.hpp" "massociate/pick.hpp"
/// @brief Defines a pick.  This need not be associated with a seismic event.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Pick
{
public:
    /// @name Constructors
    /// @{
    /// @brief Constructor
    Pick();
    /// @brief Copy constructor.
    /// @param[in] pick  The class from which to initialize this class.
    Pick(const Pick &pick);
    /// @brief Move constructor.
    /// @param[in,out] pick  The class from which to initialize this class.
    ///                      On exit, pick's behavior is undefined.
    Pick(Pick &&pick) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy assignment operator.
    /// @param[in] pick  The pick class to copy to this.
    /// @result A deep copy of the pick.
    Pick& operator=(const Pick &pick);
    /// @brief Move assignment operator.
    /// @param[in,out] pick  The pick class whose memory will be moved
    ///                      to this.  On exit, pick's behavior is
    ///                      undefined.
    /// @result The memory from pick moved to this.
    Pick& operator=(Pick &&pick) noexcept;
    /// @}

    /// @name Destructors
    /// @{
    ~Pick();
    /// @brief Resets the class.
    void clear() noexcept;
    /// @}

    /// @name Waveform Identifier
    /// @{
    /// @brief Sets the station, network, channel, and location code on which
    ///        the pick was observed.
    /// @param[in] waveformIdentifier  The station, network, channel, and
    ///                                location.  This cannot be empty.
    /// @throws std::invalid_argument if \c waveformIdentifier.isEmpty()
    ///         is true.
    void setWaveformIdentifier(
        const MAssociate::WaveformIdentifier &waveformIdentifier);
    /// @result Gets the WaveformIdentifier on which the pick was observed.
    /// @throws std::runtime_error if the WaveformIdentifier was not set.
    /// @sa \c haveWaveformIdentifier()
    [[nodiscard]] MAssociate::WaveformIdentifier getWaveformIdentifier() const;
    /// @result True indicates that the WaveformIdentifier was set.
    [[nodiscard]] bool haveWaveformIdentifier() const noexcept;
    /// @}
 
    /// @name Phase Name
    /// @{
    /// @brief Sets the phase name.
    /// @param[in] phaseName  The name of the seismic phase.  
    /// @throws std::invalid_argument if the phaseName is empty.
    void setPhaseName(const std::string &phaseName);
    /// @brief Gets the name of the seismic phase for this pick.
    /// @result The name of the seismic phase.
    /// @throws std::invalid_argument if the phase name is not set.
    /// @sa \c havePhaseName()
    [[nodiscard]] std::string getPhaseName() const;
    /// @result True indicates that the phase name is set.
    [[nodiscard]] bool havePhaseName() const noexcept;
    /// @}

    /// @name Pick Time
    /// @{
    /// @brief Sets the pick time.
    /// @param[in] pickTime  The pick time in UTC seconds since the epoch.
    void setTime(double pickTime) noexcept;
    /// @brief Gets the pick time.
    /// @result The pick time in UTC seconds since the epoch.
    /// @throws std::runtime_error if the pick time was not set.
    /// @sa \c haveTime()
    [[nodiscard]] double getTime() const;
    /// @result True indicates that the pick time was set.
    [[nodiscard]] bool haveTime() const noexcept;
    /// @}

    /// @name Identifier
    /// @{
    /// @brief Sets a unique pick identifier number.  This likely will
    ///        be assigned by the real-time system but an algorithm as simple
    ///        as counting picks can be used.
    /// @param[in] identifier   The unique pick identifier.
    void setIdentifier(uint64_t identifier) noexcept;
    /// @result The pick's unique identification number.
    /// @throws std::runtime_error if \c haveIdentifier() is false.
    [[nodiscard]] uint64_t getIdentifier() const;
    /// @result True indicates that the identifier has been set.
    [[nodiscard]] bool haveIdentifier() const noexcept;
    /// @}

    /// @name Weight
    /// @{
    /// @brief Defines the pick's standard deviation (error) in seconds.
    /// @param[in] std   The standard deviation in seconds.
    /// @throws std::invalid_argument if the standard deviation is not positive.
    /// @note When migrating with a Gaussian distribution this is the canonical
    ///       standard deviation.  When migrating with a boxcar function the
    ///       standard deviation is \f$ \sigma = \sqrt{ \frac{(B-A)^2}{12} } \f$
    ///       where $W = B - A$ is the boxcar's width.  I find it convenient to
    ///       think of parameterizing boxcars based on width, so, in this case I
    ///       use \f$ \sigma = \frac{W}{2 \sqrt{3}} \f$. 
    void setStandardDeviation(double std);
    /// @result The pick's standard deviation in seconds.
    [[nodiscard]] double getStandardDeviation() const noexcept;
    /// @result The pick's weight in 1/seconds.  This is computed from
    ///         the standard deviation i.e., \f$ \frac{1}{\sigma} \f$.
    [[nodiscard]] double getWeight() const noexcept;
    /// @}

    /// @name Polarity
    /// @{
    /// @brief Sets the polarity.
    /// @param[in] polarity  The pick's polarity.
    void setPolarity(MAssociate::Polarity polarity) noexcept;
    /// @brief Gets the pick's polarity.
    [[nodiscard]] MAssociate::Polarity getPolarity() const noexcept;
    /// @brief Sets the polarities weight.  This is useful when the polarity
    ///        is produced by an ML algorithm and has an accompanying posterior
    ///        probability.
    /// @param[in] weight  The polarity's weight.
    /// @throws std::invalid_argument if the weight is not positive.
    void setPolarityWeight(double weight);
    /// @result The polarity's weight.
    [[nodiscard]] double getPolarityWeight() const noexcept;
    /// @} 

    /// @name Static Correction
    /// @{
    /// @brief Sets the static currection such that the modeled pick time is
    ///        \f$ T_{modeled} = T_0 + T + T_s \f$ 
    ///        where \f$ T_0 \f$ is the origin time, \f$ T \f$ is the travel
    ///        time and \f$ T_s \f$ is the static correction.
    /// @param[in] correction   The static correction in seconds.
    void setStaticCorrection(double correction) noexcept;
    /// @result The static correction in seconds.
    [[nodiscard]] double getStaticCorrection() const noexcept;
    /// @}
private:
    class PickImpl;
    std::unique_ptr<PickImpl> pImpl;
};
}
#endif
