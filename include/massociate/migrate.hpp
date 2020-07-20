#ifndef MASSOCIATE_MIGRATE_HPP
#define MASSOCIATE_MIGRATE_HPP
#include <memory>
namespace MAssociate
{
class Pick;
class MigrationParameters;
/*!
 * @class Migrate "migrate.hpp" "massociate/migrate.hpp"
 * @brief This class migrates the differential pick times.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license. 
 */
template<class T = float>
class Migrate
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor
     */
    Migrate();
    /*!
     * @brief Copy constructor.
     * @param[in] migrate  The migration class from which to initialize
     *                     this class.
     */
    Migrate(const Migrate &migrate);
    /*!
     * @brief Move constructor.
     * @param[in,out] migrate  The migration class from which to initialize
     *                         this class.  On exit, migrate's behavior is
     *                         undefined.
     */
    Migrate(Migrate &&migrate) noexcept;
    /*! @} */

    /*!
     * @brief Copy assignment
     * @param[in] migrate   The migrate class to copy to this.
     * @result A deep copy of the migrate class.
     * @note This will also copy travel time tables which can be very
     *       resource intensive.
     */
    Migrate& operator=(const Migrate &migrate);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] migrate   On input this is the migrate class whose
     *                          memory will be moved to this.  On exit,
     *                          migrate's behavior is undefined.
     * @return The memory from the migrate class moved to this.
     */
    Migrate& operator=(Migrate &&migrate) noexcept;

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor
     */
    ~Migrate();
    /*!
     * @brief Releases memory and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Step 1: Initialize
     * @{
     */
    /*!
     * @brief Initializes the class.
     * @param[in] parameters   The migration parameters.  The required
     *                         parameters must be set.
     * @throws std::invalid_argument if the required parameters are not set.
     * @note You will have to add the travel time tables in a second step.
     */
    void initialize(const MigrationParameters &parameters);
    /*!
     * @return True indicates that class is initialized.
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
     * @note There that there is no static correction added.
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

    /*! @name Step 4: Migrate
     * @{
     */
    /*!
     * @brief Migrates the picks.
     * @throws std::runtime_error if the class is not initialized or there are
     *         less than two picks.
     * @sa \c isInitialized(), \c getNumberOfPicks()
     */
    void migrate();
    /*! @} */
    /*!
     * @return True indicates that the migration image has been computed.
     */
    [[nodiscard]] bool haveImage() const noexcept;
    /*! @} */

    /*! @name Step 6: Get Optima
     * @{
     */
    /*!
     * @result The index of and value of the maximum in the migration image.
     * @throws std::runtime_error if the image is not yet computed.
     * @sa \c haveImage()
     */
    [[nodiscard]] std::pair<int, T> getImageMaximum() const;
    /*!
     * @result This is the travel time in seconds from each pick to the
     *         maximum in the migration image for each pick.  Note that the
     *         static corrections have been added into the maximum.
     * @throws std::runtime_error if the image is not yet coputed.
     * @sa \c getImageMaximum()
     */
    [[nodiscard]] std::vector<double> getTravelTimesToMaximum() const;
    /*!
     * @param[in] normalize   If true then normalize so that the contributions
     *                        sum to unity.
     * @result Each pick's contribution to the maximum.  This can be useful
     *         for tie breaking when a waveform identifier/phase pair 
     *         contributes twice but the associator must choose one or the
     *         other.
     */
    [[nodiscard]] std::vector<double> 
        getContributionToMaximum(const bool normalize = true) const;
    /*! @} */
private:
    class MigrateImpl;
    std::unique_ptr<MigrateImpl> pImpl;
};
}
#endif
