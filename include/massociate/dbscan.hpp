#ifndef MASSOCIATE_DBSCAN_HPP
#define MASSOCIATE_DBSCAN_HPP
#include <memory>
#include <massociate/clusterer.hpp>
namespace MAssociate
{
/// @brief Interface to the DBSCAN clustering algorithm.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class DBSCAN : public IClusterer
{
public:
    /// @name Constructors
    /// @{
 
    /// @brief Constructor.
    DBSCAN();
    /// @}

    /// @name Initialization
    /// @{

    /// @brief Sets the DBSCAN parameters.
    /// @param[in] epsilon     The maximum difference between two points
    ///                        in a cluster.
    /// @param[in] minSamples  The minimum number of samples required to form
    ///                        a cluster.
    /// @throws std::invalid_argument if epsilon is negative or
    ///         minSamples is not positive.
    void initialize(double epsilon, int minSamples);
    /// @brief Determines if the class is initialized.
    /// @result True indicates that the class is initialized.
    bool isInitialized() const noexcept;
    /// @}

    /// @name Set Data to Cluster
    /// @{

    /// @brief Sets the data to cluster.
    /// @param[in] nObservations  The number of observations (rows).
    /// @param[in] nFeatures      The number of features (columns).
    /// @param[in] X              The feature matrix.  This is an
    ///                           [nObservations x nFeatures] matrix stored in
    ///                           row major format.
    /// @throws std::runtime_error if the class is not initialized.
    /// @throws std::invalid_argument if nObservations or nFeatures is less than
    ///         one or X is NULL.
    /// @sa \c isInitialized()
    void setData(int nObservations, int nFeatures, const std::vector<double> &X) override;
    /// @brief Sets the data to cluster.
    /// @param[in] nObservations  The number of observations (rows).
    /// @param[in] nFeatures      The number of features (columns).
    /// @param[in] X              The feature matrix.  This is an
    ///                           [nObservations x nFeatures] matrix stored in
    ///                           row major format.
    /// @param[in] weights        These are the weights of each observation.
    ///                           This is an array whose dimension
    ///                           [nObservations] and all values must be
    ///                           positive.
    /// @throws std::runtime_error if the class is not initialized.
    /// @throws std::invalid_argument if nObservations or nFeatures is less than
    ///         one, or X is NULL, or weights has an entry which is not postiive.
    /// @sa \c isInitialized()
    void setWeightedData(int nObservations, int nFeatures, const std::vector<double> &X,
                         const std::vector<double> &weights);
    /// @brief Determines if the data was set.
    /// @result True indicates that the data was set.
    [[nodiscard]] bool haveData() const noexcept;
    /// @}

    /// @name Clustering
    /// @{
 
    /// @brief Clusters the data using the DBSCAN algorithm.
    /// @throws std::runtime_error if no data was set or the class was not
    ///         initialized.
    /// @sa \c haveData(), \c isInitialized()
    void cluster() override final;
    /// @}

    /// @name Results
    /// @{

    /// @brief Gets the number of clusters.
    [[nodiscard]] int getNumberOfClusters() const noexcept override final;
    /// @brief Sets the data to cluster.
    /// @result The cluster label of each observation.
    /// @throws std::runtime_error if the clusters were not computed.
    /// @sa \c haveLabels()
    [[nodiscard]] std::vector<int> getLabels() const override final;
    /// @brief Determines if the labels were computed. 
    /// @result True indicates that DBSCAN was called and the labels
    ///         were extracted.
    [[nodiscard]] bool haveLabels() const noexcept override final;
    /// @}

    /// @name Destructors
    /// @{
 
    /// @brief Resets the class and releases all memory.
    void clear() noexcept;
    /// @brief Destructor.
    ~DBSCAN() override;
    /// @}

    /// @result The clustering algorithm.
    [[nodiscard]] IClusterer::Algorithm getAlgorithm() const noexcept override final;

    DBSCAN(const DBSCAN &) = delete;
    DBSCAN(DBSCAN &&) noexcept = delete;
    DBSCAN& operator=(const DBSCAN &) = delete;
    DBSCAN& operator=(DBSCAN &&) noexcept = delete;
private:
    class DBSCANImpl;
    std::unique_ptr<DBSCANImpl> pImpl;
};
}
#endif
