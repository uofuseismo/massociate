#ifndef MASSOCIATE_DIVIDED_RECTANGLES_HPP
#define MASSOCIATE_DIVIDED_RECTANGLES_HPP
#include <memory>
#include <massociate/optimizer.hpp>
#include <umps/logging/log.hpp>
namespace MAssociate
{
 class IMigrator;
 class Arrival;
 class Event;
}
namespace MAssociate
{
class DividedRectangles : public IOptimizer
{
public:
    /// @brief Constructor.
    DividedRectangles();
    explicit DividedRectangles(std::shared_ptr<UMPS::Logging::ILog> &logger);

    /// @brief Sets the migration engine.
    //void setMigrator(std::unique_ptr<IMigrator> &&migrator);
    /// @result Releases the migration engine.
    //[[nodiscard]] std::unique_ptr<IMigrator> releaseMigrator(); 
    /// @result True indicates the migration engine was set.
    //[[nodiscard]] bool haveMigrator() const noexcept;

    /// @result True indicates the arrivals were set.
    //[[nodiscard]] bool haveArrivals() const noexcept;
    void enableSearchDepth() noexcept;
    void disableSearchDepth() noexcept;
    [[nodiscard]] bool searchDepth() const noexcept;


    /// @brief To build an event, this number of arrivals is required to
    ///        contribute to the image maximum.
    //void setMinimumNumberOfArrivalsToBuildEvent(int nArrivals);
    /// @result After migrating at at a point, this number of arrivals is
    ///         required to contribute in order to build the event.
    //[[nodiscard]] int getMinimumNumberOfArrivalsToBuildEvent() const noexcept;

    void setRefinement(double refinement);
    [[nodiscard]] double getRefinement() const noexcept;

    /// @brief Sets the number of initial objective function evaluations to
    ///        perform in the depth search.
    /// @param[i] nEvaluations  The number of objective function evaluations.
    /// @throws std::invalid_argument if n is not positive.
    void setNumberOfInitialObjectiveFunctionEvaluations(int nEvaluations);
    /// @result The number of objective function evaluations.
    [[nodiscard]] int getNumberOfInitialObjectiveFunctionEvaluations() const noexcept;

    /// @brief Sets the number of objective function evaluations.  Loosely
    ///        speaking this is the number of iterations for DIRECT.
    /// @param[i] nEvaluations  The number of objective function evaluations.
    /// @throws std::invalid_argument if n is not positive.
    void setNumberOfObjectiveFunctionEvaluations(int nEvaluations);
    /// @result The number of objective function evaluations.
    [[nodiscard]] int getNumberOfObjectiveFunctionEvaluations() const noexcept;

    /// @brief Enables normalization of the simple box-boundaries.
    void enableNormalization() noexcept;
    /// @brief Disables normalization of the simple box-boundaries.
    void disableNormalization() noexcept;
    /// @result True indicates the simple box-boundaries will be normalized
    ///         prior to searching.
    /// @note By default normalization is set.
    [[nodiscard]] bool normalize() const noexcept;

    /// @brief Enables biasing towards local search.  This implements the
    ///        "A locally-biased form of the DIRECT algorithm" by 
    ///        Gablonsky and Kelley and is more efficient when the 
    ///        objective function does not have too many local minima.
    void enableLocallyBias() noexcept;
    /// @brief Disables biasing towards local search.  This implements
    ///        the classic "Lipshcitzian optimization without Lipschitz
    ///        constant" method of Jones et al.
    void disableLocallyBias() noexcept;
    /// @result Indicates whether or not to use the locally biased DIRECT
    ///         algorithm.
    [[nodiscard]] bool locallyBias() const noexcept;

    /// @brief Sets the depth at which the migration optimization will
    ///        be performed.
    /// @param[in] depth   The depth in meters.
    void setDepth(double depth);
    /// @result The depth at which the migration will be performed.
    [[nodiscard]] double getDepth() const noexcept;

    /// @brief Sets the search extent in X.
    /// @param[in] extentInX   The lower and upper extent to search in x in
    ///                        meters.  This should be a subset of the region's
    ///                        extent.
    void setExtentInX(const std::pair<double, double> &extentInX);
    /// @result The extent to search in x.  By default this is the region's
    ///         extent.
    [[nodiscard]] std::pair<double, double> getExtentInX() const;
    /// @brief Sets the search extent in Y.
    /// @param[in] extentInY   The lower and upper extent to search in y in
    ///                        meters.  This should be a subset of the region's
    ///                        extent.
    void setExtentInY(const std::pair<double, double> &extentInY);
    /// @result The extent to search in y.  By default this is the region's
    ///         extent.
    [[nodiscard]] std::pair<double, double> getExtentInY() const;
    /// @brief Sets the search extent in Z.
    /// @param[in] extentInZ   The lower and upper extent to search in z in
    ///                        meters.
    /// @throws std::invalid_argument if the search extent is not between
    ///         -8,500 and 800,000 meters.
    void setExtentInZ(const std::pair<double, double> &extentInZ);
    /// @result The extent to search in z.
    [[nodiscard]] std::pair<double, double> getExtentInZ() const;

    /// @brief Performs the particle swarm optimization.
    void optimize() override final;

    [[nodiscard]] double getOptimalValue() const override final;
    [[nodiscard]] bool haveOptimum() const noexcept override final;
    [[nodiscard]] std::tuple<double, double, double> getOptimalHypocenter() const override final;
    [[nodiscard]] double getOptimalOriginTime() const override final;
    [[nodiscard]] std::vector<std::pair<Arrival, double>> getContributingArrivals() const override final;

    /// @name Destructors
    /// @{

    /// @brief Destructor.
    ~DividedRectangles();
    /// @}

    DividedRectangles(const DividedRectangles &) = delete;
    DividedRectangles& operator=(const DividedRectangles &) = delete;
private:
    class DividedRectanglesImpl;
    std::unique_ptr<DividedRectanglesImpl> pImpl; 
};
}
#endif
