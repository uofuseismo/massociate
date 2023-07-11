#include <iostream>
#include <limits>
#include <vector>
#include <algorithm>
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#include "private/pageRank.hpp"
using namespace MAssociate;

namespace
{

[[maybe_unused]]
void gemv(const int lda, const int nRows, const int nCols,
          const double *A, const double *x, double *y,
          const double alpha = 1,
          const double beta = 0,
          const bool lSymmetric = false)
{
    cblas_dgemv(CblasRowMajor, CblasNoTrans, nRows, nCols, alpha, A, lda,
                x, 1, beta, y, 1);
} 

[[maybe_unused]]
void gemv(const int lda, const int nRows, const int nCols,
          const float *A, const float *x, float *y,
          const float alpha = 1,
          const float beta = 0,
          const bool lSymmetric = false)
{
    cblas_sgemv(CblasRowMajor, CblasNoTrans, nRows, nCols, alpha, A, lda,
                x, 1, beta, y, 1);
}

template<class T>
std::vector<T> computeSymmetricAdjacency(const std::vector<std::pair<int, int>> &edges,
                                         const std::vector<T> &weights)
{
    auto nEdges = static_cast<int> (edges.size());
    std::vector<int> nodes(2*nEdges);
    for (int edge=0; edge<nEdges; ++edge)
    {
        nodes[2*edge]   = edges[edge].first;
        nodes[2*edge+1] = edges[edge].second;
    }
    std::stable_sort(nodes.begin(), nodes.end());
    auto ip = std::unique(nodes.begin(), nodes.end());
    nodes.resize(std::distance(nodes.begin(), ip));
    auto nRows = static_cast<int> (nodes.size());
    auto nCols = nRows;
    auto lda = nCols;
    // Fill upper triangle
    std::vector<T> M(lda*nCols, 0);
    for (int edge=0; edge<nEdges; ++edge)
    {
        auto row = std::lower_bound(nodes.begin(), nodes.end(),
                                    edges[edge].first);
        auto col = std::lower_bound(nodes.begin(), nodes.end(),
                                    edges[edge].second);
        int i = std::distance(nodes.begin(), row);
        int j = std::distance(nodes.begin(), col);
        auto idx = i*lda + j;
        M[idx] = M[idx] + weights[edge]; 
    }
    // Fill lower half
    for (int j=0; j<nCols; ++j)
    {
        for (int i=0; i<j-1; ++i)
        {
            M[j*lda + i] = M[i*lda + j];
        }
    }
    // Compute row sum by computing M*ones
    std::vector<T> ones(nCols, 1);
    std::vector<T> rowSum(nRows);
    gemv(lda, nRows, nCols, M.data(), ones.data(), rowSum.data(), 1, 0, true);
    // Normalize
    for (int i=0; i<nRows; ++i)
    {
        #pragma omp simd
        for (int j=0; j<nRows; ++j){M[i*lda + j] = M[i*lda + j]/rowSum[i];}
    } 
}

}


class PageRank::PageRankImpl
{
    std::vector<std::pair<int, int>> mEdges;
    std::vector<double> mWeights;
};

PageRank::PageRank() :
    pImpl(std::make_unique<PageRankImpl> ())
{
}

PageRank::~PageRank() = default;

int PageRankUtilities::getLargestEigenvalueAndEigenvector(
      const int nRows, double A[], const int lda,
      std::vector<double> *eigenValues, std::vector<double> *eigenVectors,
      const double abstol)
{
    const lapack_int il = nRows;
    const lapack_int iu = nRows;
    const double vl = 0;
    const double vu = 0;
    lapack_int ldz = iu - il + 1; //nRows;
    std::vector<double> evectors(nRows*ldz, 0);
    std::vector<double> evalues(nRows);
    lapack_int mFound; // iu - il + 1
    auto info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'I', 'U',
                               nRows, A, lda, vl, vu, il, iu, abstol,
                               &mFound, evalues.data(), evectors.data(),
                               ldz, NULL);
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
