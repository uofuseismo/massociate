#ifndef MASSOCIATE_ASSOCIATORPARAMETERS_HPP
#define MASSOCIATE_ASSOCIATORPARAMETERS_HPP
#include <memory>
#include "massociate/enums.hpp"
namespace MAssociate
{
class AssociatorParameters
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor
     */
    AssociatorParameters();
    /*!
     * @brief Copy constructor.
     * @param[in] parameters   The parameters class from which to initialize
     *                         this class.
     */
    AssociatorParameters(const AssociatorParameters &parameters);
    /*!
     * @brief Move constructor.
     * @param[in,out] parameters  The parameters class from which to initialize
     *                            this class.  On exit, parameters's behavior
     *                            is undefined.
     */
    AssociatorParameters(AssociatorParameters &&parameters) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] parameters   The parameters class to copy to this.
     * @return A deep copy of the associator parameters class.
     */
    AssociatorParameters& operator=(const AssociatorParameters &parameters);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] parameters   The parameters class whose memory will be
     *                             moved to this.  On exit, parameters's
     *                             behavior is undefined.
     * @return The memory from parameters moved to this.
     */
    AssociatorParameters& operator=(AssociatorParameters &&parameters) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor
     */
    ~AssociatorParameters();
    /*!
     * @brief Resets the class and clears all memory.
     */
    void clear() noexcept;
    /*! @} */

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
     * @brief Sets the minimum number of arrivals in a cluster required
     *        to nucleate an event.
     * @param[in] minArrivals  The minimum number of arrivals required to
     *                         nucleate an event.
     * @throws std::invalid_argument if minArrivals is not positive.
     * @note Technically four arrivals are required to locate an event
     *       so this number should be at least four.
     */
    void setMinimumNumberOfArrivalsToNucleate(int minArrivals);
    /*!
     *
     * @return
     */
    [[nodiscard]] int getMinimumNumberOfArrivalsToNucleate() const noexcept;

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


    /*! @name DBSCAN Clustering
     * @{
     */
    /*!
     * @brief DBSCAN adds elements to the cluster that are within epsilon of
     *        any point in the cluster.  This algorithm clusters on origin times
     *        so any point within epsilon seconds of any origin time in a cluster
     *        will be assimilated.
     * @param[in] epsilon   The maximum time between origin times in a cluster
     *                      measured in seconds.
     * @throws std::invalid_argument if this is negative.
     */
    void setDBSCANEpsilon(double epsilon);
    /*!
     * @result The maximum time (seconds0 between origin times in a cluster.
     */
    [[nodiscard]] double getDBSCANEpsilon() const noexcept;
    /*!
     * @brief Sets the number of origin times required to form a cluster.
     * @param[in] minimumClusterSize   The minimum number of origin times
     *                                 to form a cluster (seed an event).
     * @throws std::invalid_argument if this is not at least 2.
     * @note Since this will seed the event you'll like want this number to
     *       be at least 5 or 6.  Moreover, since DBSCAN does not perform
     *       causality checks this value should be greater than or equal
     *       to \c getMinimumNumberOfArrivalsToNucleate().
     */
    void setDBSCANMinimumClusterSize(int minimumClusterSize);
    /*!
     * @result The number of origin times required to form a cluster.
     */
    [[nodiscard]] int getDBSCANMinimumClusterSize() const noexcept;
    /*! @} */

    /*! @name PageRank
     * @{
     */
    /*!
     * @brief Sets the number of iterations in the  PageRank power iteration.
     * @param[in] nIterations   The number of iterations.
     * @throws std::invalid_argument if this is not positive.
     */
    void setPageRankNumberOfIterations(int nIterations);
    /*!
     * @result The number of page rank iterations.
     */
    [[nodiscard]] int getPageRankNumberOfIterations() const noexcept;
    /*!
     * @brief Sets the damping in the PageRank power iteration
     *        \f$
     *           \textbf{r}_{i+1}
     *         = d \mathcal{M} \textbf{r}_i + \frac{1 - d}{N} I
     *        \f$.
     *        where \$ d \f$ is the damping factor and \f$ \textbf{r}_i \f$
     *        are the page ranks at the i'th iteration.
     * @param[in] damping  The damping factor which must be in the range
     *                     \f$ [0,1 ] \f$.
     * @throws std::invalid_argument if the damping factor is out of range.
     */
    void setPageRankDampingFactor(double damping);
    /*!
     * @result The damping factor.
     */
    [[nodiscard]] double getPageRankDampingFactor() const noexcept;
    /*! @} */

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
     */
    void setTileSize(int tileSize);
    /*!
     * @result The tile size.
     */
    [[nodiscard]] int getTileSize() const noexcept;
private:
    class AssociatorParametersImpl;
    std::unique_ptr<AssociatorParametersImpl> pImpl;
};
}
#endif
