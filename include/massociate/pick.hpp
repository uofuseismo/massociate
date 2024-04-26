#ifndef MASSOCIATE_PICK_HPP
#define MASSOCIATE_PICK_HPP
#include <memory>
#include <chrono>
#include <massociate/arrival.hpp>
namespace MAssociate
{
 class WaveformIdentifier;
}
namespace MAssociate
{
/// @class Pick "pick.hpp" "massociate/pick.hpp"
/// @brief Defines a pick.  This need not be associated with a seismic event.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Pick
{
public:
    /// @brief This defines a hint as to the phase of the pick.
    ///        This is usually informed by the machine learning detector.
    enum class PhaseHint
    {
         Unknown, /*!< There is no phase hint. */
         P,       /*!< The phase hint is a primary arrival. */
         S,       /*!< The phase hint is a shear-wave arrival. */
         Noise    /*!< The phase hint that this is a noise pick. */
    };
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

    /// @brief Destructor.
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

     /// @brief Sets the phase hint as informed by the detector.
    /// @param[in] phaseHint   The inferred seismic phase type.
    void setPhaseHint(PhaseHint phaseHint) noexcept;
    /// @brief Sets the phase hint as informed by the detector.
    /// @param[in] phaseHint  The inferred seismic phase type.
    /// @throws std::invalid_argument if the phaseName is empty.
    void setPhaseHint(const std::string &phaseHint);
    /// @brief Gets the name of the seismic phase hint for this pick.
    /// @result The name of the seismic phase hint.
    /// @throws std::invalid_argument if the phase name is not set.
    /// @sa \c havePhaseName()
    [[nodiscard]] std::string getPhaseHint() const;
    /// @result True indicates that the phase name is set.
    [[nodiscard]] bool havePhaseHint() const noexcept;
    /// @}

    /// @name Pick Time
    /// @{

    /// @brief Sets the pick time.
    /// @param[in] time  The pick time (UTC) in microseconds since the epoch.
    void setTime(const std::chrono::microseconds &time) noexcept;
    /// @brief Sets the pick time.
    /// @param[in] time  The pick time (UTC) in seconds since the epoch.
    void setTime(double time) noexcept;
    /// @brief Gets the pick time.
    /// @result The pick time (UTC) in microseconds since the epoch.
    /// @throws std::runtime_error if the pick time was not set.
    /// @sa \c haveTime()
    [[nodiscard]] std::chrono::microseconds getTime() const;
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

    /// @brief Defines the pick's standard error (deviation) in seconds.
    /// @param[in] standardError  The standard error in seconds.
    /// @throws std::invalid_argument if the standard deviation is not positive.
    /// @note When migrating with a Gaussian distribution this is the canonical
    ///       standard deviation.  When migrating with a boxcar function the
    ///       standard deviation is \f$ \sigma = \sqrt{ \frac{(B-A)^2}{12} } \f$
    ///       where $W = B - A$ is the boxcar's width.  I find it convenient to
    ///       think of parameterizing boxcars based on width, so, in this case I
    ///       use \f$ \sigma = \frac{W}{2 \sqrt{3}} \f$. 
    void setStandardError(double standardError);
    /// @result The pick's standard error in seconds.
    [[nodiscard]] double getStandardError() const noexcept;
    /// @result The pick's weight in 1/seconds.  This is computed from
    ///         the standard error i.e., \f$ \frac{1}{\sigma} \f$.
    [[nodiscard]] double getWeight() const noexcept;
    /// @}

    /// @name Polarity
    /// @{

    /// @brief Sets the polarity.
    /// @param[in] firstMotion  The pick's first motion.
    void setFirstMotion(MAssociate::Arrival::FirstMotion firstMotion) noexcept;
    /// @brief Gets the pick's first motion.
    [[nodiscard]] MAssociate::Arrival::FirstMotion getFirstMotion() const noexcept;
    /// @brief Sets the first motion's weight.  This is useful when the first 
    ///        motion is produced by an ML algorithm and has an accompanying
    ///        posterior probability.
    /// @param[in] weight  The first motion's weight.
    /// @throws std::invalid_argument if the weight is not positive.
    void setFirstMotionWeight(double weight);
    /// @result The first motion's weight.
    [[nodiscard]] double getFirstMotionWeight() const noexcept;
    /// @} 
private:
    class PickImpl;
    std::unique_ptr<PickImpl> pImpl;
};
std::ostream& operator<<(std::ostream &os,
                         const MAssociate::Pick &pick);
}
#endif
