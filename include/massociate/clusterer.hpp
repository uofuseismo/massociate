#ifndef MASSOCIATE_CLUSTERER_HPP
#define MASSOCIATE_CLUSTERER_HPP
#include <vector>
namespace MAssociate
{
class IClusterer
{
public:
    enum class Algorithm
    {
        DBSCAN
    };
public:

    virtual void setData(int nObservations, int nFeatures, const std::vector<double> &X) = 0;

    virtual void cluster() = 0;

    /// @result The cluster labels for each example.
    [[nodiscard]] virtual std::vector<int> getLabels() const = 0;
    /// @result True indicates the cluster labels have been computed.
    [[nodiscard]] virtual bool haveLabels() const noexcept = 0;
    /// @result The number of clusters.
    [[nodiscard]] virtual int getNumberOfClusters() const = 0;

    [[nodiscard]] virtual Algorithm getAlgorithm() const noexcept = 0;

    virtual ~IClusterer() = default;
};
}
#endif
