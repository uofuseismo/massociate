#ifndef SFF_MINISEED_WAVEFORMIDENTIFIER_HPP
#define SFF_MINISEED_WAVEFORMIDENTIFIER_HPP
#include <memory>
namespace MAssociate
{
/*!
 * @class WaveformIdentifier "waveformIdentifier.hpp" "massociate/waveformIdentifier.hpp"
 * @brief Defines a station, network, channel, location name.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class WaveformIdentifier
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    WaveformIdentifier();
    /*!
     * @brief Copy constructor.
     * @param[in] waveformIdentifier  The waveform id class from which to
     *                                initialize.
     */
    WaveformIdentifier(const WaveformIdentifier  &waveformIdentifier);
    /*!
     * @brief Move constructor.
     * @param[in,out] waveformIdentifier  The  class to initialize from. On
     *                                    exit wawveformIdentifier's behavior is
     *                                    undefined.
     */
    WaveformIdentifier(WaveformIdentifier &&waveformIdentifier) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] waveformIdentifier  The waveform identifier to copy.
     * @result A deep copy of sncl.
     */
    WaveformIdentifier& operator=(const WaveformIdentifier &waveformIdentifier);
    /*!
     * @brief Move assignment operator.
     * @param[in] waveformIdentifier  The waveform identifier whose memory will
     *                                be moved to this.  On exit,
     *                                waveformIdentifier's is undefined.
     */
    WaveformIdentifier& operator=(WaveformIdentifier &&waveformIdentifier) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~WaveformIdentifier();
    /*!
     * @brief Clears the waveform identifier class.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Network
     * @{
     */
    /*!
     * @brief Sets the network name.
     * @param[in] name  The name of the network.
     */
    void setNetwork(const std::string &name) noexcept;
    /*!
     * @brief Gets the network name.
     * @result The name of the network.
     */
    [[nodiscard]] std::string getNetwork() const noexcept;
    /*! @} */

    /*! @name Station
     * @{
     */
    /*!
     * @brief Sets the station name.
     * @param[in] name  The name of the station.
     */
    void setStation(const std::string &name) noexcept;
    /*!
     * @brief Gets the station name.
     * @result The name of the station.
     */
    [[nodiscard]] std::string getStation() const noexcept;
    /*! @} */

    /*! @name Channel
     * @{
     */
    /*!
     * @brief Sets the channel name.
     * @param[in] name  The name of the channel.
     */
    void setChannel(const std::string &name) noexcept;
    /*!
     * @brief Gets the channel name.
     * @result The name of the channel.
     */
    [[nodiscard]] std::string getChannel() const noexcept;
    /*! @} */ 

    /*! @name Location Code
     * @{
     */
    /*!
     * @brief Sets the location code.
     * @param[in] name  The name of the location code.
     */
    void setLocationCode(const std::string &name) noexcept;
    /*!
     * @brief Gets the location code.
     * @result The name of the location code.
     */
    [[nodiscard]] std::string getLocationCode() const noexcept;
    /*! @} */

    /*!
     * @brief Convenience routine to determine if the network, station,
     *        channel location are all empty.
     * @result True indicates that there is no SNCL information.
     *         False indicates that at least one of the SNCL elements is
     *         defined.
     */
    [[nodiscard]] bool isEmpty() const noexcept;
private:
    class WaveformIdentifierImpl;
    std::unique_ptr<WaveformIdentifierImpl> pImpl;
};
bool operator==(const WaveformIdentifier &lhs, const WaveformIdentifier &rhs) noexcept;
bool operator!=(const WaveformIdentifier &lhs, const WaveformIdentifier &rhs) noexcept;
}
#endif
