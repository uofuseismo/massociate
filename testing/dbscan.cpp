#include <cstdio>
#include <cstdlib>
#include <vector>
#include <mkl_lapacke.h>
#include "private/dbscan.hpp"
#include "private/pageRank.hpp"
#include <gtest/gtest.h>
namespace
{
/*
int getLargestEigenvalueAndEigenvector(
      const int nRows, double A[], const int lda,
      std::vector<double> *eigenValues, std::vector<double> *eigenVectors,
      const double abstol = std::numeric_limits<double>::epsilon()*100)
{
    const lapack_int il = nRows;
    const lapack_int iu = nRows;
    const double vl = 0;
    const double vu = 0;
    lapack_int ldz = iu - il + 1; //nRows;
    std::vector<double> evectors(nRows*ldz, 0);
    std::vector<double> evalues(nRows);
    std::vector<int> isuppz(2*std::max(iu-il+1,1));
    lapack_int mFound; // iu - il + 1
    auto info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'I', 'U',
                               nRows, A, lda, vl, vu, il, iu, abstol,
                               &mFound, evalues.data(), evectors.data(),
                               ldz, isuppz.data());
    if (info < 0)
    {
        std::runtime_error("Argument error: invalid argument: "
                         + std::to_string(info));
    }
    else if (info > 0)
    {
        std::cerr << "Algorithm failed" << std::endl;
        eigenValues->resize(0);
        eigenVectors->resize(0);
        return -1;
    }
    eigenValues->resize(mFound, 0);
    std::copy(evalues.begin(), evalues.begin() + mFound, eigenValues->begin());
    eigenVectors->resize(mFound*nRows, 0);
    for (int i=0; i<nRows; ++i)
    {
        auto row = eigenVectors->data() + i*mFound;
        for (int j=0; j<mFound; ++j)
        {
            row[j] = evectors[i*ldz+j];
        }
    } 
    return 0;
}
*/

TEST(dbscan, dbscan)
{
    int nObs = 6;
    int nFeatures = 2;
    std::vector<double> X(nObs*nFeatures);
    X[0] = 1;
    X[1] = 2;

    X[2] = 2;
    X[3] = 2;

    X[4] = 2;
    X[5] = 2;

    X[6] = 8;
    X[7] = 7;

    X[8] = 8;
    X[9] = 8;

    X[10] = 25;
    X[11] = 80;

    std::vector<int> labelsRef({0, 0, 0, 1, 1, -1});

    MAssociate::DBSCAN dbscan;
    double epsilon = 3;
    int minSamples = 2;
    EXPECT_NO_THROW(dbscan.initialize(epsilon, minSamples));
    EXPECT_NO_THROW(dbscan.setData(nObs, nFeatures, X.data()));
    EXPECT_NO_THROW(dbscan.cluster());
    auto nClusters = dbscan.getNumberOfClusters();
    EXPECT_EQ(nClusters, 2);
    //printf("nClusters: %d\n", nClusters);
    auto labels = dbscan.getLabels();
    EXPECT_EQ(labels.size(), labelsRef.size());
    for (int i=0; i<static_cast<int> (labels.size()); ++i)
    {
        EXPECT_EQ(labels[i], labelsRef[i]);
    }
    // Remove the outlier
    EXPECT_NO_THROW(dbscan.setData(nObs-1, nFeatures, X.data()));
    EXPECT_NO_THROW(dbscan.cluster());
    nClusters = dbscan.getNumberOfClusters();
    //printf("nClusters: %d\n", nClusters);
    EXPECT_EQ(nClusters, 2);
    labels = dbscan.getLabels();
    EXPECT_EQ(labels.size(), labelsRef.size()-1);
    for (int i=0; i<static_cast<int> (labels.size()); ++i)
    {
        EXPECT_EQ(labels[i], labelsRef[i]);
    }
}
TEST(pageRank, largestEigenvalue)
{
    const int nRows = 5;
    const int lda = nRows;
    std::vector<double> A({0.67, -0.20,  0.19, -1.06,  0.46,
                          -0.20,  3.82, -0.13,  1.06, -0.48,
                           0.19, -0.13,  3.27,  0.11,  1.10,
                          -1.06,  1.06,  0.11,  5.86, -0.98,
                           0.46, -0.48,  1.10, -0.98,  3.54});
    double eigenvalueRef = 6.934791784860098;
    std::vector<double> eigenvectorRef({0.18243522, -0.35606835, 0.10215486,
                                       -0.84049795, 0.35079952});
    std::vector<double> evalues, evectors;
    MAssociate::PageRankUtilities::getLargestEigenvalueAndEigenvector(
        nRows, A.data(), lda, &evalues, &evectors);
    EXPECT_EQ(1, static_cast<int> (evalues.size()));
    EXPECT_NEAR(evalues[0], eigenvalueRef, 1.e-5); 
    for (int i=0; i<5; ++i)
    {
        EXPECT_NEAR(evectors[i], eigenvectorRef[i], 1.e-5);
    } 
}

}
