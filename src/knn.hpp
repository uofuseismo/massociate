#ifndef KNN_HPP
#define KNN_HPP
#include <vector>
#include <string>
#include <cmath>
#ifndef NDEBUG
#include <cassert>
#endif
#include <oneapi/dal/algo/knn.hpp>
#include <oneapi/dal/table/homogen.hpp>
#include <oneapi/dal/table/row_accessor.hpp>
namespace
{
std::map<std::string, std::vector<std::string>>
    cluster(const std::vector<double> &x, const std::vector<double> &y,
            const std::vector<std::string> &stationNames,
            const int nNeighborsIn = 5)
{
    std::map<std::string, std::vector<std::string>> result;
#ifndef NDEBUG
    assert(x.size() == y.size());
#endif
    auto nStations = static_cast<int> (x.size());
    if (nStations < 1){return result;}
    int nNeighbors = std::min(nStations, nNeighborsIn);
    oneapi::dal::knn::descriptor
    <
        double,
        oneapi::dal::knn::method::brute_force,
        oneapi::dal::knn::task::search,
        oneapi::dal::minkowski_distance::descriptor<double>
    > descriptor(nNeighbors);
    std::vector<double> trainingFeaturesMatrix(2*nStations);
    for (int i = 0; i < nStations; ++i)
    {
        trainingFeaturesMatrix[2*i] = x[i];
        trainingFeaturesMatrix[2*i + 1] = y[i];
    }
    constexpr int64_t nFeatures{2};
    auto nTrainingObservations = static_cast<int64_t> (nStations);
    const auto XTrain
        = oneapi::dal::homogen_table::wrap<double>
         (trainingFeaturesMatrix.data(), nTrainingObservations, nFeatures);
    auto model = oneapi::dal::train(descriptor, XTrain).get_model();
    const auto trainingResult
         = oneapi::dal::infer(descriptor, XTrain, model);
    auto trainingIndices = trainingResult.get_indices();
    oneapi::dal::row_accessor<const double> accessor{trainingIndices};
    for (int64_t row = 0; row < nTrainingObservations; ++row)
    {
        const auto rowValues = accessor.pull({row, row + 1});
        std::vector<std::string> myNeighbors;
        for (int j = 0; j < nNeighbors; ++j)
        {
            auto stationIndex = static_cast<int> (std::round(rowValues[j]));
            myNeighbors.push_back(stationNames[stationIndex]);
        }
        result.insert(std::pair {stationNames[row], std::move(myNeighbors)});
    }
    return result;
}
}
#endif
