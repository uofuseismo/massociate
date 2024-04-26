#include <cmath>
#include <vector>
#include <limits>
#include <functional>
#include <pagmo/pagmo.hpp>
#include <uLocator/position/geographicRegion.hpp>
#include <uLocator/optimizers/originTime.hpp>
#include <umps/logging/standardOut.hpp>
#include <nlopt.hpp>
#include "massociate/dividedRectangles.hpp"
#include "massociate/migrator.hpp"
#include "massociate/arrival.hpp"
#include "massociate/event.hpp"

using namespace MAssociate;

namespace
{
struct FixedDepthDoubleDifferenceMigration
{
    struct Point3D
    {
        double x{0};
        double y{0};
        double value{0};
    }; 
    explicit FixedDepthDoubleDifferenceMigration(
        const nlopt::algorithm algorithm) :
        mOptimizer(algorithm, 2) 
    {
        mOptimizer.set_max_objective(mObjectiveFunction);
    }
    [[nodiscard]] double operator()(const unsigned int n,
                                    const double x[], double *gradient) const
    {
        double value{std::numeric_limits<double>::lowest()};
#ifndef NDEBUG
        assert(x != nullptr);
        assert(n == static_cast<unsigned int> (mParameters));
#endif
        try
        {
            if (gradient == nullptr)
            {
                value = mMigrationFunction(x[0], x[1]);
                if (mSaveHistory)
                {   
                    mHistory.push_back(Point3D {x[0], x[1], value});
                }
            }
            else
            {
                throw std::runtime_error("Cannot yet compute gradient"); 
            }
        }
        catch (const std::exception &e) 
        {
            auto errorMessage = "Problem for source at f.s. (x,y) = ("
                              + std::to_string(x[0]) + "," 
                              + std::to_string(x[1])
                              + ").  Failed with: "  + std::string {e.what()};
            std::cerr << errorMessage << std::endl;
        }
        return value;
    }   
    /// @brief Sets the search boundaries.
    void setSearchBoundaries(
        const std::vector<double> &lowerBoundaries,
        const std::vector<double> &upperBoundaries)
    {   
        if (static_cast<int> (lowerBoundaries.size()) != mParameters)
        {
            throw std::invalid_argument("lowerBoundaries size != 2");
        }
        if (lowerBoundaries.size() != upperBoundaries.size())
        {
            throw std::invalid_argument("upperBoundaries size != 2");
        }
        for (int i = 0; i < static_cast<int> (lowerBoundaries.size()); ++i)
        {
            if (lowerBoundaries[i] >= upperBoundaries[i])
            {
                throw std::invalid_argument("l[i] >= u[i]");
            }
        }
        mOptimizer.set_lower_bounds(lowerBoundaries);
        mOptimizer.set_upper_bounds(upperBoundaries);
        mLowerBounds = lowerBoundaries;
        mUpperBounds = upperBoundaries;
    } 
    /// @brief Sets the maximum number of objective function eavluations.
    void setNumberOfObjectiveFunctionEvaluations(const int nEvaluations)
    {   
        if (nEvaluations < 1)
        {
            throw std::invalid_argument("nEvaluations must be positive");
        }
        mOptimizer.set_maxeval(nEvaluations);
        if (mSaveHistory)
        {
            mHistory.reserve(nEvaluations);
        }
    }
    /// @brief Resets objective function counters
    void resetCounters()
    {
        mObjectiveFunctionEvaluations = 0;
        mGradientEvaluations = 0;
        if (mSaveHistory)
        {
            mHistory.clear();
            mHistory.reserve(mOptimizer.get_maxeval());
        }
    }
    nlopt::opt mOptimizer;
    MAssociate::IMigrator *mMigrator{nullptr}; 
    std::function<double (unsigned int, const double *, double *) >
        mObjectiveFunction
        {
            std::bind(&FixedDepthDoubleDifferenceMigration::operator(),
                      this,
                      std::placeholders::_1,
                      std::placeholders::_2,
                      std::placeholders::_3)
        };
    std::function< double(const double, const double) >
        mMigrationFunction;
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    std::vector<double> mLowerBounds;
    std::vector<double> mUpperBounds;
    mutable std::vector<Point3D> mHistory;
    int mParameters{2}; 
    int mObjectiveFunctionEvaluations{0};
    int mGradientEvaluations{0};
    bool mSaveHistory{false};
};

struct FreeDepthDoubleDifferenceMigration
{
    struct Point4D
    {
        double x{0};
        double y{0};
        double z{0};
        double value{0};
    };

    explicit FreeDepthDoubleDifferenceMigration(
        const nlopt::algorithm algorithm) :
        mOptimizer(algorithm, 3)
    {
        mOptimizer.set_max_objective(mObjectiveFunction);
    }
    [[nodiscard]] double operator()(int n, const double x[],
                                    double *gradient) const
    {
        double value{std::numeric_limits<double>::lowest()};
#ifndef NDEBUG
        assert(x != nullptr);
        assert(n == static_cast<unsigned int> (mParameters));
#endif
        try
        {
            if (gradient == nullptr)
            {
                value = mMigrationFunction(x[0], x[1], x[2]);
  //std::cout << x[0] << " " << x[1] << " " << x[2] << " " << value << std::endl;

                if (mSaveHistory)
                {
                    mHistory.push_back(Point4D {x[0], x[1], x[2], value});
                }
            }
            else
            {
                throw std::runtime_error("Cannot yet compute gradient");
            }
        }
        catch (const std::exception &e)
        {
            auto errorMessage = "Problem for source at (x,y,z) = ("
                              + std::to_string(x[0]) + ","
                              + std::to_string(x[1]) + ","
                              + std::to_string(x[2]) + ","
                              + ").  Failed with: "  + std::string {e.what()};
            std::cerr << errorMessage << std::endl;
        }
        return value;
    }
    /// @brief Sets the search boundaries.
    void setSearchBoundaries(
        const std::vector<double> &lowerBoundaries,
        const std::vector<double> &upperBoundaries)
    {   
        if (static_cast<int> (lowerBoundaries.size()) != mParameters)
        {
            throw std::invalid_argument("lowerBoundaries size != 3");
        }
        if (lowerBoundaries.size() != upperBoundaries.size())
        {
            throw std::invalid_argument("upperBoundaries size != 3");
        }
        for (int i = 0; i < static_cast<int> (lowerBoundaries.size()); ++i)
        {
            if (lowerBoundaries[i] >= upperBoundaries[i])
            {
                throw std::invalid_argument("l[i] >= u[i]");
            }
        }
        mOptimizer.set_lower_bounds(lowerBoundaries);
        mOptimizer.set_upper_bounds(upperBoundaries);
        mLowerBounds = lowerBoundaries;
        mUpperBounds = upperBoundaries;
    }   
    /// @brief Sets the maximum number of objective function eavluations.
    void setNumberOfObjectiveFunctionEvaluations(const int nEvaluations)
    {
        if (nEvaluations < 1)
        {
            throw std::invalid_argument("nEvaluations must be positive");
        }
        mOptimizer.set_maxeval(nEvaluations);
        if (mSaveHistory)
        {
            mHistory.reserve(nEvaluations);
        }
    }
    /// @brief Resets objective function counters
    void resetCounters()
    {
        mObjectiveFunctionEvaluations = 0;
        mGradientEvaluations = 0;
        if (mSaveHistory)
        {
            mHistory.clear();
            mHistory.reserve(mOptimizer.get_maxeval());
        }
    }

    nlopt::opt mOptimizer;
    MAssociate::IMigrator *mMigrator{nullptr};
    std::function<double (unsigned int, const double *, double *) >
        mObjectiveFunction
        {
            std::bind(&FreeDepthDoubleDifferenceMigration::operator(),
                      this,
                      std::placeholders::_1,
                      std::placeholders::_2,
                      std::placeholders::_3)
        };
    std::function< double(const double, const double, const double) >
        mMigrationFunction;
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    std::vector<double> mLowerBounds;
    std::vector<double> mUpperBounds;
    mutable std::vector<Point4D> mHistory;
    int mParameters{3}; 
    int mObjectiveFunctionEvaluations{0};
    int mGradientEvaluations{0};
    bool mSaveHistory{false};
};

}

class DividedRectangles::DividedRectanglesImpl
{
public:
    DividedRectanglesImpl(std::shared_ptr<UMPS::Logging::ILog> logger = nullptr) :
        mLogger(logger)
    {
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> ();
        }
    }
    [[nodiscard]] double migrateFixedDepth(const double x, const double y) const
    {
        return mMigrator->evaluate(x, y, mDepth);
    }
    [[nodiscard]] double migrateFreeDepth(const double x,
                                          const double y,
                                          const double z) const
    {
        return mMigrator->evaluate(x, y, z);
    }   
    std::function<double (const double, const double)>
    mFixedDepthMigrationFunction
    {
        std::bind(&DividedRectanglesImpl::migrateFixedDepth,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2)
    };
    std::function<double (const double, const double, const double)>
    mFreeDepthMigrationFunction
    {
        std::bind(&DividedRectanglesImpl::migrateFreeDepth,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2,
                  std::placeholders::_3)
    };  
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    IMigrator *mMigrator{nullptr};
    Event mEvent;
    std::pair<double, double> mExtentInX;
    std::pair<double, double> mExtentInY;
    // For my relocated/enhanced catalogs this represents 
    // the 1% percentile and 99% percentile.  Hence, this range
    // covers a lot of depths that we are likely to find.
    // The reason for using the relocated catalogs is we will be
    // using the updated velocity models and station/source corrections.
    std::pair<double, double> mExtentInZ{-1700, 22000};
    // For Yellowstone we would go from 0.9 km to 16 km but to 
    // be super-conservative I'll use -1 km.
    //std::pair<double, double> mExtentInZ{-1000, 16000};
    // Default search depths for Utah
    std::vector<double> mSearchDepths{-1000, 3000, 6400, 20000};
    // For YNP

    std::vector<std::pair<Arrival, double>> mContributingArrivals;
    double mOptimalMigrationValue{0};
    double mBestLatitude{0};
    double mBestLongitude{0};
    double mBestDepth{0};
    double mDepth{6000};
    double mRefinement{25000};
    int mNumberOfInitialObjectiveFunctionEvaluations{200};
    int mNumberOfObjectiveFunctionEvaluations{500};
    bool mNormalize{true};
    bool mLocallyBias{true};
    bool mHaveExtentInX{false};
    bool mHaveExtentInY{false};
    bool mHaveExtentInZ{true};
    bool mHaveOptimum{false};
    bool mSearchDepth{false};
};

/// Constructor
DividedRectangles::DividedRectangles() :
    IOptimizer (),
    pImpl(std::make_unique<DividedRectanglesImpl> ())
{
}

DividedRectangles::DividedRectangles(
    std::shared_ptr<UMPS::Logging::ILog> &logger) :
    IOptimizer(),
    pImpl(std::make_unique<DividedRectanglesImpl> (logger))
{
}


/*
/// Move constructor
DividedRectangles::DividedRectangles(ParticleSwarm &&direct) noexcept
{
    *this = std::move(direct);
}

/// Move assignment
DividedRectangles& DividedRectangles::operator=(ParticleSwarm &&direct) noexcept
{
    if (&pso == this){return *this;}
    pImpl = std::move(direct.pImpl);
    return *this;
}
*/

/// Search depth
void DividedRectangles::enableSearchDepth() noexcept
{
    pImpl->mSearchDepth = true;
}

void DividedRectangles::disableSearchDepth() noexcept
{
    pImpl->mSearchDepth = false;
}

/// Search in depth?
bool DividedRectangles::searchDepth() const noexcept
{
    return pImpl->mSearchDepth;
}

/// Refinement
void DividedRectangles::setRefinement(double refinement)
{
    if (refinement <= 0)
    {
        throw std::invalid_argument("Refinement must be positive");
    }
    pImpl->mRefinement = refinement;
}

double DividedRectangles::getRefinement() const noexcept
{
    return pImpl->mRefinement;
}

/// Number of object function evaluations
void DividedRectangles::setNumberOfInitialObjectiveFunctionEvaluations(
    const int n)
{
    if (n < 1)
    {   
        throw std::runtime_error(
           "Number of objective function evaluations must be positive");
    }   
    pImpl->mNumberOfInitialObjectiveFunctionEvaluations = n;
}

int DividedRectangles::getNumberOfInitialObjectiveFunctionEvaluations()
    const noexcept
{
    return pImpl->mNumberOfInitialObjectiveFunctionEvaluations;
}

/// Number of object function evaluations
void DividedRectangles::setNumberOfObjectiveFunctionEvaluations(const int n)
{
    if (n < 1)
    {
        throw std::runtime_error(
           "Number of objective function evaluations must be positive");
    }
    pImpl->mNumberOfObjectiveFunctionEvaluations = n;
}

int DividedRectangles::getNumberOfObjectiveFunctionEvaluations()
    const noexcept
{
    return pImpl->mNumberOfObjectiveFunctionEvaluations;
}

/// Locally bias?
void DividedRectangles::enableLocallyBias() noexcept
{
    pImpl->mLocallyBias = true;
}

void DividedRectangles::disableLocallyBias() noexcept
{
    pImpl->mLocallyBias = false;
}

bool DividedRectangles::locallyBias() const noexcept
{
    return pImpl->mLocallyBias;
}

/// Normalize?
void DividedRectangles::enableNormalization() noexcept
{
    pImpl->mNormalize = true;
}

void DividedRectangles::disableNormalization() noexcept
{
    pImpl->mNormalize = false;
}

bool DividedRectangles::normalize() const noexcept
{
    return pImpl->mNormalize;
}

/// The depth at which to perform the migration
void DividedRectangles::setDepth(const double depth)
{
    if (depth < -8500 || depth > 800000)
    {
        throw std::invalid_argument("Depth must be in range [-8500, 800,000]");
    }
    pImpl->mDepth = depth;
}

/// The depth
double DividedRectangles::getDepth() const noexcept
{
    return pImpl->mDepth;
}

/// Sets the extent
void DividedRectangles::setExtentInX(
    const std::pair<double, double> &extent)
{
    if (extent.second <= extent.first)
    {
        throw std::invalid_argument("extent.second <= extent.first in x");
    }
    pImpl->mExtentInX = extent;
    pImpl->mHaveExtentInX = true;
}


void DividedRectangles::setExtentInY(
    const std::pair<double, double> &extent)
{
    if (extent.second <= extent.first)
    {
        throw std::invalid_argument("extent.second <= extent.first in y");
    }
    pImpl->mExtentInY = extent;
    pImpl->mHaveExtentInY = true;
}

std::pair<double, double> DividedRectangles::getExtentInX() const
{
    if (!pImpl->mHaveExtentInX)
    {
        throw std::runtime_error("Extent in x not set");
    }
    return pImpl->mExtentInX;
}

std::pair<double, double> DividedRectangles::getExtentInY() const
{
    if (!pImpl->mHaveExtentInY)
    {
        throw std::runtime_error("Extent in y not set");
    }
    return pImpl->mExtentInY;
}

void DividedRectangles::setExtentInZ(
    const std::pair<double, double> &extent)
{
    if (extent.second <= extent.first)
    {   
        throw std::invalid_argument("extent.second <= extent.first in z");
    }
    if (extent.first < -8500)
    {
        throw std::invalid_argument("Minimum depth must be > -8500 m");
    }
    if (extent.second > 800000)
    {
        throw std::invalid_argument("Maximum depth must be < 800,000 m");
    }
    pImpl->mExtentInZ = extent;
    pImpl->mHaveExtentInZ = true;
}

std::pair<double, double> DividedRectangles::getExtentInZ() const
{
    if (!pImpl->mHaveExtentInZ)
    {
        throw std::runtime_error("Extent in z not set");
    }
    return pImpl->mExtentInZ;
}


/// Performs the optimization
void DividedRectangles::optimize()
{
    pImpl->mHaveOptimum = false;
    if (!haveMigrator())
    {
        throw std::runtime_error("Migration engine not set");
    }
    if (!haveArrivals())
    {
        throw std::runtime_error("Arrivals not set");
    }
    pImpl->mMigrator = IOptimizer::getMigratorHandle();
    // Figure out DIRECT algorithm
    auto directAlgorithm = nlopt::GN_DIRECT_L;
    if (normalize())
    {
        if (locallyBias())
        {
            pImpl->mLogger->debug("Using DIRECT locally biased algorithm");
            directAlgorithm = nlopt::GN_DIRECT_L;
        }
        else
        {
            pImpl->mLogger->debug("Using DIRECT classic algorithm");
            directAlgorithm = nlopt::GN_DIRECT;
        }
    }
    else
    {
        if (locallyBias())
        {
            pImpl->mLogger->debug(
                "Using DIRECT locally biased unscaled algorithm");
            directAlgorithm = nlopt::GN_DIRECT_L_NOSCAL;
        }
        else
        {
            pImpl->mLogger->debug("Using DIRECT unscaled algorithm");
            directAlgorithm = nlopt::GN_DIRECT_NOSCAL;
        }
    }
    // Ensure extent is set
    auto region = pImpl->mMigrator->getGeographicRegion();
    if (!pImpl->mHaveExtentInX)
    {
        setExtentInX(region->getExtentInX());
    }   
    if (!pImpl->mHaveExtentInY)
    {   
        setExtentInY(region->getExtentInY());
    }
    auto [x0, x1] = getExtentInX();
    auto [y0, y1] = getExtentInY();
    auto evaluationDepth = getDepth();
    double xOptimum{0};
    double yOptimum{0}; 
    double zOptimum{0};
    double fOptimum{0};
    // Make sure the migrator will not save the contribution history
    pImpl->mMigrator->disableSaveArrivalContributionList();
    // Do a series of 2D searches
    for (const auto &searchDepth : pImpl->mSearchDepths)
    {
        // Instantiate the NLOpt class 
        ::FixedDepthDoubleDifferenceMigration optimizer2D{directAlgorithm};
        optimizer2D.mMigrationFunction = pImpl->mFixedDepthMigrationFunction;
        optimizer2D.mLogger = pImpl->mLogger;
        // Set bounds 
        optimizer2D.setSearchBoundaries(std::vector<double> {x0, y0},
                                        std::vector<double> {x1, y1} );
        // Number of objective function evaluations
        optimizer2D.setNumberOfObjectiveFunctionEvaluations(
            getNumberOfInitialObjectiveFunctionEvaluations());

        std::vector<double> xLocation2D(2, 0);
        double fOptimumForDepth{0};
        try
        {
            setDepth(searchDepth);
            pImpl->mLogger->debug("Optimizing at "
                                + std::to_string(searchDepth)
                                + " depth ...");
            optimizer2D.mOptimizer.optimize(xLocation2D, fOptimumForDepth);
        }
        catch (const std::exception &e)
        {
            auto errorMessage = "2D DIRECT failed with: "
                              + std::string {e.what()};
            pImpl->mLogger->error(errorMessage);
            throw std::runtime_error(errorMessage);
        }
        if (fOptimumForDepth >= fOptimum)
        {
            xOptimum = xLocation2D.at(0);
            yOptimum = xLocation2D.at(1);
            zOptimum = getDepth();
            fOptimum = fOptimumForDepth;
        }
    }
    setDepth(evaluationDepth);
    pImpl->mLogger->debug("Initial search finished at: ("
                        + std::to_string(xOptimum) + ","
                        + std::to_string(yOptimum) + ","
                        + std::to_string(zOptimum) + ")");
    // Refine and optimize
    x0 = std::max(xOptimum - getRefinement(), getExtentInX().first);
    x1 = std::min(xOptimum + getRefinement(), getExtentInX().second);
    y0 = std::max(yOptimum - getRefinement(), getExtentInY().first);
    y1 = std::min(yOptimum + getRefinement(), getExtentInY().second);
    auto [z0, z1] = getExtentInZ();
    // Instantiate the NLOpt class 
    ::FreeDepthDoubleDifferenceMigration optimizer3D{directAlgorithm};
    optimizer3D.mMigrationFunction = pImpl->mFreeDepthMigrationFunction;
    optimizer3D.mLogger = pImpl->mLogger;
    // Set bounds
    optimizer3D.setSearchBoundaries(std::vector<double> {x0, y0, z0},
                                    std::vector<double> {x1, y1, z1} );
    // Number of objective function evaluations
    optimizer3D.setNumberOfObjectiveFunctionEvaluations(
        getNumberOfObjectiveFunctionEvaluations());
    std::vector<double> xLocation3D(3, 0);
    xLocation3D[0] = 0.5*(x0 + x1); //std::min(xOptimum, std::max(x0, xOptimum));
    xLocation3D[1] = 0.5*(y0 + y1); //std::min(yOptimum, std::max(y0, yOptimum)); 
    xLocation3D[2] = 0.5*(z0 + z1); //std::min(zOptimum, std::max(z0, zOptimum));
    double fOptimum3D{0};
    try
    {
        pImpl->mLogger->debug("Optimizing...");
        optimizer3D.mOptimizer.optimize(xLocation3D, fOptimum3D);
    }
    catch (const std::exception &e) 
    {
        auto errorMessage = "3D DIRECT failed with: "
                          + std::string {e.what()};
        pImpl->mLogger->error(errorMessage);
        throw std::runtime_error(errorMessage);
    }
    if (fOptimum3D >= fOptimum)
    {
        xOptimum = xLocation3D.at(0);
        yOptimum = xLocation3D.at(1);
        zOptimum = xLocation3D.at(2);
        fOptimum = fOptimum3D;
    }
    else
    {
        pImpl->mLogger->debug("2D DIRECT beat 3D DIRECT");
    }
    // Search through the known locations
    auto knownCandidateEventsScores
        = pImpl->mMigrator->evaluateAtKnownLocations();
    if (!knownCandidateEventsScores.empty())
    {
        auto bestKnownCandidateEventScore
            = std::max_element(knownCandidateEventsScores.begin(),
                               knownCandidateEventsScores.end());
        auto index
            = std::distance(knownCandidateEventsScores.begin(),
                            bestKnownCandidateEventScore);
        auto bestKnownCandidateEventOptimum
            = knownCandidateEventsScores.at(index);
        //std::cout << bestKnownCandidateEventOptimum << " " << fOptimum << std::endl;
        if (bestKnownCandidateEventOptimum > fOptimum)
        {
            pImpl->mLogger->info("Predetermined event location beat PSO");
            auto knownLocation = pImpl->mMigrator->getKnownSearchLocation(index);
            xOptimum = knownLocation->x();
            yOptimum = knownLocation->y();
            zOptimum = knownLocation->z();
            fOptimum = bestKnownCandidateEventOptimum;
        }
    }
    
    // Pick a winner and extract the hypocenter and contributing arrivals
    //auto optimumLocation = newPopulation.champion_x();
    //auto xOptimum = optimumLocation.at(0);
    //auto yOptimum = optimumLocation.at(1);
    //auto zOptimum = evaluationDepth;
    auto [latitude, longitude]
         = region->localToGeographicCoordinates(xOptimum, yOptimum);
//std::cout << latitude << " " << longitude << " " << zOptimum << " " << fOptimum << std::endl;
    pImpl->mBestLatitude = latitude;
    pImpl->mBestLongitude = longitude;
    pImpl->mBestDepth = zOptimum;
    try
    {
        pImpl->mMigrator->enableSaveArrivalContributionList();
        pImpl->mOptimalMigrationValue
            = pImpl->mMigrator->evaluate(xOptimum, yOptimum, zOptimum);
        pImpl->mContributingArrivals
            = pImpl->mMigrator->getContributingArrivals();
        pImpl->mHaveOptimum = true;
    }   
    catch (const std::exception &e) 
    {   
        pImpl->mLogger->debug("PSO could not build event: "
                            + std::string {e.what()});
        pImpl->mHaveOptimum = false;
    }
    pImpl->mMigrator = nullptr;
}

/// Have optimum?
bool DividedRectangles::haveOptimum() const noexcept
{
    return pImpl->mHaveOptimum;
}

/// Optimal value
double DividedRectangles::getOptimalValue() const
{
    if (!haveOptimum())
    {
        throw std::runtime_error("Optimal location not computed");
    }
    return pImpl->mOptimalMigrationValue;
}

/// Get the optimal hypocenter
std::tuple<double, double, double> DividedRectangles::getOptimalHypocenter() const
{
    if (!haveOptimum())
    {
        throw std::runtime_error("Optimal location not computed");
    }
    return std::tuple {pImpl->mBestLatitude,
                       pImpl->mBestLongitude,
                       pImpl->mBestDepth};
}

/// Optimal value
std::vector<std::pair<Arrival, double>>
    DividedRectangles::getContributingArrivals() const
{
    if (!haveOptimum())
    {   
        throw std::runtime_error("Optimal location not computed");
    }   
    return pImpl->mContributingArrivals;
}

/// Destructor
DividedRectangles::~DividedRectangles() = default; 

