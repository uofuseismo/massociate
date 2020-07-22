#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <array>
#include <algorithm>
#include <mkl.h>
#include "bilinearInterpolation.hpp"

namespace
{
/// @brief Linearly interpolates f(x,y).
/// @param[in] x     The x location.
/// @param[in] y     The y location.
/// @param[in] x1    The start point of cell containing x.
/// @param[in] x2    The end point of the cell containing x.
/// @param[in] y1    The start point of the cell containing y.
/// @param[in] y2    The end point of the cell containing y.
/// @param[in] fq11  The value of f at (x1,y1).
/// @param[in] fq21  The value of f at (x2,y1). 
/// @param[in] fq12  The value of f at (x1,y2).
/// @parampin[ fq22  The value of f at (x2,y2).
/// @result The value of f linearly interpolated at (x,y).
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
#pragma omp declare simd
template<typename T>
T bilinearInterpolate(const T x, const T y,
                      const T x1, const T x2,
                      const T y1, const T y2,
                      const T fq11, const T fq21,
                      const T fq12, const T fq22)
{
    auto dx = x2 - x1;
    auto dy = y2 - y1;
    auto x_m_x1 = x - x1;
    auto x2_m_x = x2 - x;
    auto y_m_y1 = y - y1;
    auto y2_m_y = y2 - y;
    auto fxy = fq11*(x2_m_x*y2_m_y)
             + fq21*(x_m_x1*y2_m_y)
             + fq12*(x2_m_x*y_m_y1)
             + fq22*(x_m_x1*y_m_y1);
    fxy = fxy/(dx*dy);
    return fxy;
}
 
}

/// Double implementation
template<>
class BilinearInterpolation<double>::BilinearInterpolationImpl
{
public:
    void createTasksForNonUniformGrid()
    {
        auto nx = static_cast<MKL_INT> (mXPoints.size());
        MKL_INT status;
        status = dfdNewTask1D(&mTaskX,
                              nx, mXPoints.data(), DF_QUASI_UNIFORM_PARTITION,
                              1, mXPoints.data(), DF_NO_HINT);
        if (status != DF_STATUS_OK)
        {
            throw std::runtime_error("Failed to create x data fitting task");
        }
        mHaveTaskX = true;

        auto ny = static_cast<MKL_INT> (mYPoints.size());
        status = dfdNewTask1D(&mTaskY,
                              ny, mYPoints.data(), DF_QUASI_UNIFORM_PARTITION,
                              1, mYPoints.data(), DF_NO_HINT);
        if (status != DF_STATUS_OK)
        {
            throw std::runtime_error("Failed to create y data fitting task");
        }
        mHaveTaskY = true;
    }
    void createTasksForUniformGrid()
    {
        auto nx = static_cast<MKL_INT> (mNumberOfXPoints);
        auto status = dfdNewTask1D(&mTaskX,
                                   nx, mXEndPoints.data(), DF_UNIFORM_PARTITION,
                                   1, mXEndPoints.data(), DF_NO_HINT);
        if (status != DF_STATUS_OK)
        {
            throw std::runtime_error("Failed to create x data fitting task");
        }
        mHaveTaskX = true;

        auto ny = static_cast<MKL_INT> (mNumberOfYPoints);
        status = dfdNewTask1D(&mTaskY,
                              ny, mYEndPoints.data(), DF_UNIFORM_PARTITION,
                              1, mYEndPoints.data(), DF_NO_HINT);
        if (status != DF_STATUS_OK)
        {
            throw std::runtime_error("Failed to create y data fitting task");
        }
        mHaveTaskY = true;
    }
    void searchPoints(const MKL_INT nq, const double xq[], const double yq[])
    {
        MKL_INT status = DF_STATUS_OK;
        // Search x
        mCellX.resize(nq, 0);
        if (std::is_sorted(xq, xq + nq))
        {
            status = dfdSearchCells1D(mTaskX, DF_METHOD_STD,
                                      nq, xq, DF_SORTED_DATA,
                                      DF_NO_APRIORI_INFO, mCellX.data());
        }
        else
        {
            status = dfdSearchCells1D(mTaskX, DF_METHOD_STD,
                                      nq, xq, DF_NO_HINT,
                                      DF_NO_APRIORI_INFO, mCellX.data());
        }
        if (status != DF_STATUS_OK)
        {
            throw std::runtime_error("Failed to search xq");
        }
        // Search y
        mCellY.resize(nq, 0);
        if (std::is_sorted(yq, yq + nq))
        {
            status = dfdSearchCells1D(mTaskY, DF_METHOD_STD,
                                      nq, yq, DF_SORTED_DATA,
                                      DF_NO_APRIORI_INFO, mCellY.data());
        }
        else
        {
            status = dfdSearchCells1D(mTaskY, DF_METHOD_STD,
                                      nq, yq, DF_NO_HINT,
                                      DF_NO_APRIORI_INFO, mCellY.data());
        }
        if (status != DF_STATUS_OK)
        {
            throw std::runtime_error("Failed to search yq");
        }
    }

    DFTaskPtr mTaskX;
    DFTaskPtr mTaskY;
    std::vector<double> mInterpolatedF;
    std::vector<double> mXPoints;
    std::vector<double> mYPoints;
    std::vector<double> mF;
    std::vector<double> mFInterpolated;
    std::vector<MKL_INT> mCellX;
    std::vector<MKL_INT> mCellY;
    std::array<double, 2> mXEndPoints;
    std::array<double, 2> mYEndPoints;
    double mDeltaX = 0;
    double mDeltaY = 0;
    int mNumberOfXPoints = 0;
    int mNumberOfYPoints = 0;
    bool mUniformGrid = false;
    bool mHaveTaskX = false;
    bool mHaveTaskY = false;
    bool mHaveInterpolatedFunction = false;
    bool mInitialized = false;
};

/// C'tor
template<class T>
BilinearInterpolation<T>::BilinearInterpolation() :
    pImpl(std::make_unique<BilinearInterpolationImpl> ())
{
}

/// Copy c'tor
template<class T>
BilinearInterpolation<T>::BilinearInterpolation(
    const BilinearInterpolation &bilin)
{
    *this = bilin;
}

/// Move c'tor
template<class T>
BilinearInterpolation<T>::BilinearInterpolation(
    BilinearInterpolation &&bilin) noexcept
{
    *this = std::move(bilin);
}

/// Copy assignment
template<class T>
BilinearInterpolation<T>&
BilinearInterpolation<T>::operator=(const BilinearInterpolation &bilin)
{
    pImpl = std::make_unique<BilinearInterpolationImpl> ();
    if (bilin.isInitialized())
    {
        // The tasks points to memory so it is safer to simply reinitialize
        auto nx = bilin.pImpl->mNumberOfXPoints;
        auto ny = bilin.pImpl->mNumberOfYPoints;
        auto nxy = nx*ny;
        if (bilin.pImpl->mUniformGrid)
        {
            std::pair<T, T> xLimits(bilin.pImpl->mXEndPoints[0],
                                    bilin.pImpl->mYEndPoints[1]);
            std::pair<T, T> yLimits(bilin.pImpl->mYEndPoints[0],
                                    bilin.pImpl->mYEndPoints[1]);
            initialize(xLimits, yLimits, nx, ny, bilin.pImpl->mF.data()); 
        }
        else
        {
            initialize(nx, bilin.pImpl->mXPoints.data(),
                       ny, bilin.pImpl->mYPoints.data(),
                       nxy, bilin.pImpl->mF.data());
        }
        // If result was computed then get it
        pImpl->mInterpolatedF = bilin.pImpl->mInterpolatedF;
        pImpl->mHaveInterpolatedFunction
            = bilin.pImpl->mHaveInterpolatedFunction; 
    }
    return *this;
}

/// Move assignment 
template<class T>
BilinearInterpolation<T>&
BilinearInterpolation<T>::operator=(BilinearInterpolation &&bilin) noexcept
{
    if (&bilin == this){return *this;}
    pImpl = std::move(bilin.pImpl);
    return *this;
}

/// Destructor
template<class T>
BilinearInterpolation<T>::~BilinearInterpolation()
{
    clear();
}

/// Initialize non-uniform grid
template<class T>
template<typename U>
void BilinearInterpolation<T>::initialize(const int nx, const U x[],
                                          const int ny, const U y[],
                                          const int nxy, const U f[])
{
    clear();
    if (nx < 2){throw std::invalid_argument("nx must be at least 2");}
    if (ny < 2){throw std::invalid_argument("ny must be at least 2");}
    if (nxy != nx*ny){throw std::invalid_argument("nx*ny != nxy");}
    if (x == nullptr){throw std::invalid_argument("x is NULL");}
    if (y == nullptr){throw std::invalid_argument("y is NULL");}
    if (f == nullptr){throw std::invalid_argument("f is NULL");}
    if (!std::is_sorted(x, x + nx))
    {
        throw std::invalid_argument("x must be sorted in increasing order");
    }
    if (!std::is_sorted(y, y + ny))
    {
        throw std::invalid_argument("y must be sorted in increasing order");
    }
    // Copy
    pImpl->mNumberOfXPoints = nx;
    pImpl->mXPoints.resize(nx);
    std::copy(x, x + nx, pImpl->mXPoints.begin());
    pImpl->mXEndPoints[0] = x[0];
    pImpl->mXEndPoints[1] = x[nx-1];

    pImpl->mNumberOfYPoints = ny;
    pImpl->mYPoints.resize(ny);
    std::copy(y, y + ny, pImpl->mYPoints.begin());
    pImpl->mYEndPoints[0] = y[0];
    pImpl->mYEndPoints[1] = y[ny-1];

    pImpl->mF.resize(nxy);
    std::copy(f, f + nxy, pImpl->mF.begin());
    // Create the tasks
    pImpl->createTasksForNonUniformGrid();
    pImpl->mUniformGrid = false;
    pImpl->mInitialized = true;
}

/// Initialize uniform grid
template<class T>
template<typename U>
void BilinearInterpolation<T>::initialize(const std::pair<U, U> &xLimits,
                                          const std::pair<U, U> &yLimits,
                                          const int nx, const int ny,
                                          const U f[])
{
    clear();
    if (nx < 2){throw std::invalid_argument("nx must be at least 2");}
    if (ny < 2){throw std::invalid_argument("ny must be at least 2");}
    if (f == nullptr){throw std::invalid_argument("f is NULL");}
    if (xLimits.first >= xLimits.second)
    {
        throw std::invalid_argument("xLimits.first must be < xLimits.second");
    }
    if (yLimits.first >= yLimits.second)
    {
        throw std::invalid_argument("yLimits.first must be < xLimits.second");
    }
    // Copy
    auto nxy = nx*ny;
    pImpl->mNumberOfXPoints = nx;
    pImpl->mNumberOfYPoints = ny;
    pImpl->mF.resize(nxy);
    std::copy(f, f + nxy, pImpl->mF.begin());
    pImpl->mXEndPoints[0] = xLimits.first;
    pImpl->mXEndPoints[1] = xLimits.second;
    pImpl->mYEndPoints[0] = yLimits.first;
    pImpl->mYEndPoints[1] = yLimits.second;
    pImpl->mDeltaX = (xLimits.second - xLimits.first)/(nx - 1);
    pImpl->mDeltaY = (yLimits.second - yLimits.first)/(ny - 1);
    // Create the tasks
    pImpl->createTasksForUniformGrid();
    pImpl->mUniformGrid = true;
    pImpl->mInitialized = true;
}

/// Reset the class
template<class T>
void BilinearInterpolation<T>::clear() noexcept
{
    if (pImpl->mHaveTaskX){dfDeleteTask(&pImpl->mTaskX);}
    if (pImpl->mHaveTaskY){dfDeleteTask(&pImpl->mTaskY);}
    std::fill(pImpl->mXEndPoints.begin(), pImpl->mXEndPoints.end(), 0);
    std::fill(pImpl->mYEndPoints.begin(), pImpl->mYEndPoints.end(), 0);
    pImpl->mInterpolatedF.clear();
    pImpl->mXPoints.clear();
    pImpl->mYPoints.clear();
    pImpl->mF.clear(); 
    pImpl->mDeltaX = 0;
    pImpl->mDeltaY = 0;
    pImpl->mNumberOfXPoints = 0;
    pImpl->mNumberOfYPoints = 0;
    pImpl->mUniformGrid = false;
    pImpl->mHaveInterpolatedFunction = false;
    pImpl->mHaveTaskX = false;
    pImpl->mHaveTaskY = false;
    pImpl->mInitialized = false;
}

/// Determines if the class is initialized
template<class T>
bool BilinearInterpolation<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Interpolates at points (xq,yq)
template<class T>
void BilinearInterpolation<T>::interpolate(const int nq,
                                           const T xq[],
                                           const T yq[])
{
    pImpl->mHaveInterpolatedFunction = false;
    if (!isInitialized())
    {
        throw std::runtime_error("Class not initialized");
    }
    if (nq < 1){return;}
    std::vector<T> fout(nq, 0);
    // Check bounds points
    auto [xmin, xmax] = std::minmax_element(xq, xq + nq);
    if (*xmin < pImpl->mXEndPoints[0] || *xmax > pImpl->mXEndPoints[1])
    {
        throw std::invalid_argument("[xmin,xmax] = ["
                                  + std::to_string(*xmin) + ","
                                  + std::to_string(*xmax) 
                                  + "] must be in range ["
                                  + std::to_string(pImpl->mXEndPoints[0]) + ","
                                  + std::to_string(pImpl->mXEndPoints[1])
                                  + "]"); 
    }
    auto [ymin, ymax] = std::minmax_element(yq, yq + nq);
    if (*ymin < pImpl->mYEndPoints[0] || *ymax > pImpl->mYEndPoints[1])
    {
        throw std::invalid_argument("[ymin,ymax] = ["
                                  + std::to_string(*ymin) + ","
                                  + std::to_string(*ymax)
                                  + "] must be in range ["
                                  + std::to_string(pImpl->mYEndPoints[0]) + ","
                                  + std::to_string(pImpl->mYEndPoints[1])
                                  + "]");
    }
    // Figure out the cells at which to interpolate
    pImpl->searchPoints(nq, xq, yq);
    // Set constants and pointers 
    auto nx = pImpl->mNumberOfXPoints;
    const auto cellX = pImpl->mCellX.data();
    const auto cellY = pImpl->mCellY.data();
    const auto f = pImpl->mF.data();
    pImpl->mInterpolatedF.resize(nq, 0);
    auto fq = pImpl->mInterpolatedF.data();
    // Interpolate on a non-uniform grid
    if (!pImpl->mUniformGrid)
    {
        const auto x = pImpl->mXPoints.data();
        const auto y = pImpl->mYPoints.data();
        #pragma omp simd
        for (int i=0; i<nq; ++i)
        {
            // Get cell indices
            auto ix = cellX[i] - 1;
            auto iy = cellY[i] - 1;
            // Get cell positions
            auto x1 = x[ix];
            auto x2 = x[ix+1];
            auto y1 = y[iy];
            auto y2 = y[iy+1];
#ifndef DNDEBUG
            assert( (xq[i] >= x1) && (xq[i] <= x2) );
            assert( (yq[i] >= y1) && (yq[i] <= y2) );
#endif
            // Get function values
            auto indx = iy*nx + ix;
            auto fq11 = f[indx];
            auto fq21 = f[indx + 1];
            auto fq12 = f[indx + nx];
            auto fq22 = f[indx + nx + 1]; 
            // Interpolate
            fq[i] = bilinearInterpolate(xq[i], yq[i], x1, x2, y1, y2,
                                        fq11, fq21, fq12, fq22);
        }
    }        
    else
    {
        // Interpolate on a uniform grid
        auto x0 = pImpl->mXEndPoints[0];
        auto y0 = pImpl->mYEndPoints[0];
        auto dx = pImpl->mDeltaX;
        auto dy = pImpl->mDeltaY;
        #pragma omp simd
        for (int i=0; i<nq; ++i)
        {
            // Get cell indices
            auto ix = cellX[i] - 1;
            auto iy = cellY[i] - 1;
            // Get cell positions
            auto x1 = x0 + dx*ix;
            auto x2 = x0 + dx*(ix + 1);
            auto y1 = y0 + dy*iy;
            auto y2 = y0 + dy*(iy + 1);
#ifndef DNDEBUG
            assert( (xq[i] >= x1) && (xq[i] <= x2) );
            assert( (yq[i] >= y1) && (yq[i] <= y2) );
#endif
            // Get function values
            auto indx = iy*nx + ix;
            auto fq11 = f[indx];
            auto fq21 = f[indx + 1];
            auto fq12 = f[indx + nx];
            auto fq22 = f[indx + nx + 1];
            // Interpolate
            fq[i] = bilinearInterpolate(xq[i], yq[i], x1, x2, y1, y2,
                                        fq11, fq21, fq12, fq22);
        }
    }
    pImpl->mHaveInterpolatedFunction = true;
} 

template<class T>
bool BilinearInterpolation<T>::haveInterpolatedFunction() const noexcept
{
    return pImpl->mHaveInterpolatedFunction;
}

template<class T>
std::vector<T> BilinearInterpolation<T>::getInterpolatedFunction() const
{
    if (!haveInterpolatedFunction())
    {
        throw std::runtime_error("interpolate not yet called");
    }
    return pImpl->mInterpolatedF;
}

template<class T>
const T *BilinearInterpolation<T>::getInterpolatedFunctionPointer() const
{
    if (!haveInterpolatedFunction())
    {
        throw std::runtime_error("interpolate not yet called");
    }
    return pImpl->mInterpolatedF.data();
}

template<class T>
int BilinearInterpolation<T>::getInterpolatedFunctionSize() const noexcept
{
    return static_cast<int> (pImpl->mInterpolatedF.size());
}

///--------------------------------------------------------------------------///
///                         Template Instantiation                           ///
///--------------------------------------------------------------------------///
template class BilinearInterpolation<double>;
template void BilinearInterpolation<double>::initialize(int nx,
                                                        const double x[],
                                                        int ny,
                                                        const double y[],
                                                        int nxy,
                                                        const double f[]);
template void BilinearInterpolation<double>::initialize(int nx,
                                                        const float x[],
                                                        int ny,
                                                        const float y[],
                                                        int nxy,
                                                        const float f[]);
template void BilinearInterpolation<double>::initialize(
    const std::pair<double, double> &xLimits,
    const std::pair<double, double> &yLimits,
    int nx, int ny, const double f[]);
template void BilinearInterpolation<double>::initialize(
    const std::pair<float, float> &xLimits,
    const std::pair<float, float> &yLimits,
    int nx, int ny, const float f[]);
/*
template void BilinearInterpolation<double>::interpolate(int nq,
                                                         const double xq[],
                                                         const double yq[]);
*/
