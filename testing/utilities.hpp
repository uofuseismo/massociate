#ifndef UTILITIES_HPP
#define UTILITIES_HPP
#include <random>
#include <vector>
#include <cmath>
namespace
{

/// Generates n uniform random numbers on the interval [lower,upper]
[[maybe_unused]]
std::vector<int> generateUniformRandomNumbers(
    const int n, const int lower, const int upper, const int seed)
{
    std::vector<int> result(n, 0);
    if (seed > 0)
    {
        std::mt19937 rng(static_cast<unsigned> (seed));
        std::uniform_int_distribution<int> dist(lower, upper);
        for (int i=0; i<n; ++i)
        {
            result[i] = dist(rng);
        }
    }
    else
    {
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<int> dist(lower, upper);
        for (int i=0; i<n; ++i)
        {
            result[i] = dist(rng);
        }
    }
    return result;
}

/// Makes a travel time table in a homogeneous medium
template<typename T>
[[maybe_unused]]
std::vector<T> makeTravelTimeTable(const T xr, const T yr, const T zr,
                                   const T dx, const T dy, const T dz,
                                   const int nx, const int ny, const int nz,
                                   const T vel,
                                   const T x0, const T y0, const T z0)
{
    std::vector<T> ttimes(nx*ny*nz, -1);
    auto slow = static_cast<T> (1./vel);
    #pragma omp simd collapse(3)
    for (int iz=0; iz<nz; ++iz)
    {
        for (int iy=0; iy<ny; ++iy)
        {
            for (int ix=0; ix<nx; ++ix)
            {
                auto x = x0 + ix*dx;
                auto y = y0 + iy*dy;
                auto z = z0 + iz*dz;
                auto diffx = x - xr;
                auto diffy = y - yr;
                auto diffz = z - zr;
                auto d = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
                ttimes[iz*nx*ny + iy*nx + ix] = slow*d;
            }
        }
    }
    return ttimes;
}

/// Defines the 3D points
template<typename T>
[[maybe_unused]]
void makePoints(const T x0, const T y0, const T z0,
                const T dx, const T dy, const T dz,
                const int nx, const int ny, const int nz,
                std::vector<T> *x,
                std::vector<T> *y,
                std::vector<T> *z)
{
    auto npts = nx*ny*nz;
    x->resize(npts, 0);
    y->resize(npts, 0);
    z->resize(npts, 0);
    auto xPtr = x->data();
    auto yPtr = y->data();
    auto zPtr = z->data();
    #pragma omp simd collapse(3)
    for (int iz=0; iz<nz; ++iz)
    {
        for (int iy=0; iy<ny; ++iy)
        {
            for (int ix=0; ix<nx; ++ix)
            {
                auto indx = iz*nx*ny + iy*nx + ix;
                xPtr[indx] = x0 + ix*dx;
                yPtr[indx] = y0 + iy*dy;
                zPtr[indx] = z0 + iz*dz;
            }
        }
    } 
}

/// Computes a travel time in a homogeneous medium
template<typename T>
[[maybe_unused]]
T computeTravelTime(const T xs, const T ys, const T zs,
                    const T xr, const T yr, const T zr,
                    const T vel)
{
    auto dx = xr - xs;
    auto dy = yr - ys;
    auto dz = zr - zs;
    return sqrt(dx*dx + dy*dy + dz*dz)/vel;
}

[[maybe_unused]]
double infinityNorm(const int n, const double x[], const double y[])
{
    double emax = 0;
    #pragma omp simd reduction(max: emax)
    for (int i=0; i<n; ++i)
    {
        emax = std::max(emax, std::abs(x[i] - y[i]));
    }
    return emax;
}

[[maybe_unused]]
double infinityNorm(const int n, const double x[], const float y[])
{
    double emax = 0;
    #pragma omp simd reduction(max: emax)
    for (int i=0; i<n; ++i)
    {
        emax = std::max(emax, std::abs(static_cast<double> (x[i] - y[i])));
    }
    return emax;
}

}
#endif
