#ifndef MASSOCIATE_MIGRATION_PARAMETERS_HPP
#define MASSOCIATE_MIGRATION_PARAMETERS_HPP
#include <memory>
#include "massociate/enums.hpp"
namespace MAssociate
{
/*!
 * @class MigrationParameters "migrationParameters.hpp" "massociate/migrationParameters.hpp"
 * @brief Defines the parameters for the migration of differential
 *        arrival times.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class MigrationParameters
{
public:
    /*! @name Constructors
     * @{
     */
    MigrationParameters();
    /*!
     * @brief Copy constructor.
     * @param[in] parameters   The class from which to initialize this class.
     */
    MigrationParameters(const MigrationParameters &parameters);
    /*!
     * @brief Move constructor.
     * @param[in,out] parameters   The class from which to initialize this
     *                             class.  On exit, parameters's behavior
     *                             is undefined.
     */
    MigrationParameters(MigrationParameters &&parameters) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] parameters   The parameters class to copy to this.
     * @return A deep copy of the parameters class.
     */
    MigrationParameters& operator=(const MigrationParameters &parameters);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] parameters   The parameters class whose memory will be
     *                             moved to this.  On exit, parameters's
     *                             behavior is undefined.
     * @return The memory from parameters moved to this.
     */
    MigrationParameters& operator=(MigrationParameters &&parameters) noexcept;
    /*! @} */
    
    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor
     */
    ~MigrationParameters();
    /*!
     * @brief Releases memory and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Required Parameters
     * @{
     */
    /*!
     * @brief Sets the number of travel time tables that will be set on the
     *        migration.
     * @param[in] nTables  The number of travel time tables.
     *                     This must be at least 2. 
     * @throws std::invalid_argument if nTables is not at least 2.
     */
    void setNumberOfTravelTimeTables(int nTables);
    /*!
     * @result The number of travel time tables.
     * @throws std::runtime_error if \c haveNumberOfTables is false.
     */
    [[nodiscard]] int getNumberOfTravelTimeTables() const;
    /*!
     * @result True indicates that this variable has been set.
     */
    [[nodiscard]] bool haveNumberOfTravelTimeTables() const noexcept;

    /*!
     * @brief Sets the number of points in each travel time table.
     * @param[in] nPoints  The number of points in the table.
     * @throws std::invalid_argument if this is not positive.
     * @note This can change the value of the tile size.  It is therefore
     *       recommended to set this prior to setting the tile size.
     */
    void setNumberOfPointsInTravelTimeTable(int nPoints); 
    /*!
     * @result The number of point in each travel time table.
     * @throws std::runtime_error if \c haveNumberOfPointsInTravelTimeTable()
     *         is false.
     */
    [[nodiscard]] int getNumberOfPointsInTravelTimeTable() const;
    /*!
     * @result True indicates that this variable has been set.
     */
    [[nodiscard]] bool haveNumberOfPointsInTravelTimeTable() const noexcept;
    /*! @} */

    /*! @name Optional Parameters
     * @{
     */
    /*!
     * @brief This is a tuning parameter.  Ideally you want the number of 
     *        tables times the tile size to fit into the processor's cache
     *        size.  If this variable is too small then memory thrashing will
     *        slow down program execution.  If this variable is too large then
     *        cache misses will slow down program execution.
     * @param[in] tileSize  The tile size.  This must be positive.
     * @throws std::invalid_argument if tileSize is not positive.
     * @note Internally this will be set to the minimum of the number of
     *       points in the table and the specified tile size.
     * @sa \c setNumberOfPointsInTravelTimeTable()
     */
    void setTileSize(int tileSize);
    /*!
     * @result The tile size.
     */
    [[nodiscard]] int getTileSize() const noexcept;

    /*!
     * @brief Sets the analytic correlation function to migrate.
     */
    void setAnalyticCorrelationFunction(
        MAssociate::AnalyticCorrelationFunction function) noexcept;
    /*!
     * @return The analytic correlation function to migrate.
     */
    [[nodiscard]] MAssociate::AnalyticCorrelationFunction
        getAnalyticCorrelationFunction() const noexcept;
    /*! @} */

private:
    class MigrationParametersImpl;
    std::unique_ptr<MigrationParametersImpl> pImpl;
};
}
#endif
