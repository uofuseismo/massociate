#ifndef MASSOCIATE_PAGERANK_HPP
#define MASSOCIATE_PAGERANK_HPP
#include <memory>
namespace MAssociate
{
class PageRank
{
public:
    PageRank();
    ~PageRank();
    void initialize(int nIterations = 100, double damping = 0.85);
    bool isInitialized();
    void addEdge(const std::pair<int, int>, double weight = 1);
    void clearEdges() noexcept;
    void compute();
    std::vector<double> getPageRanks() const;
private:
    class PageRankImpl;
    std::unique_ptr<PageRankImpl> pImpl;
};
namespace PageRankUtilities
{
/*!
 * @brief Returns the largest eigenvalue/eigenvector pair.
 * @param[in] nRows  The number of rows (columns) in A.
 * @param[in,out] A  This is a row major matrix storing a symmetric matrix.
 *                   It has dimension [nRows x lda].  On exit, A may be
 *                   modified by LAPACK.
 * @param[in] lda    The leading dimension of A.  This probably should be nRows.
 * @param[out] eigenValues   On exit, this contains the largest 
 *                           eigenvalue(s).  The number of eigenvalues are
 *                           eigenValues.size().
 * @param[out] eigenVectors  A row major matrix of eigenvectors.  This has
 *                           dimension [nRows x eigenValues.size()].
 * @param[in] abstol  The absolute tolerance to which the eigenvalue(s)
 *                    must be computed.  It is likely unnecessary to find
 *                    eigenvalue to high precision. 
 * @result 0 indicates success.
 */
int getLargestEigenvalueAndEigenvector(
      int nRows, double A[], int lda,
      std::vector<double> *eigenValues, std::vector<double> *eigenVectors,
      double abstol = std::numeric_limits<double>::epsilon()*100);
}
}
#endif
