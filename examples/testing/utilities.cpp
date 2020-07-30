#include <cmath>
#include <iostream>
#include <vector>
#include "bilinearInterpolation.hpp"
#include <gtest/gtest.h>

namespace
{

double fxy(const double x, const double y)
{
    return 2*x + 0.1*y;
}

template<typename T>
void makeGrid(const int nx, const int ny,
              const T x0, const T dx,
              const T y0, const T dy,
              std::vector<T> &x,
              std::vector<T> &y,
              std::vector<T> &f,
              bool phase1 = true)
{
    auto nxy = nx*ny;
    if (phase1)
    {
        x.resize(nx, 0);
        y.resize(ny, 0);
    }
    else
    {
        x.resize(nxy, 0);
        y.resize(nxy, 0);
    }
    f.resize(nxy, 0);
    for (int iy=0; iy<ny; ++iy)
    {
        for (int ix=0; ix<nx; ++ix)
        {
            auto indx = iy*nx + ix;
            auto ixuse = indx;
            auto iyuse = indx;
            if (phase1)
            {
                ixuse = ix;
                iyuse = iy;
            }
            x[ixuse] = x0 + ix*dx;
            y[iyuse] = y0 + iy*dy;
            f[indx] = fxy(x[ixuse], y[iyuse]);
        }
    }
}

TEST(Utilities, BilinearInterpolation)
{
    BilinearInterpolation bint;
    int nx = 5;
    int ny = 6;
    double dx = 0.2;
    double dy = 0.5;
    double x0 = 0;
    double y0 = 0;
    double x1 = x0 + (nx - 1)*dx;
    double y1 = y0 + (ny - 1)*dy;
    int nxInt = 15;
    int nyInt = 18;
    double dxInt = (x1 - x0)/(nxInt - 1);
    double dyInt = (y1 - y0)/(nyInt - 1);
    std::vector<double> x, y, f;
    makeGrid(nx, ny, x0, dx, y0, dy, x, y, f, true);
    std::vector<double> xq, yq, fqRef;
    makeGrid(nxInt, nyInt, x0, dxInt, y0, dyInt, xq, yq, fqRef, false);
    bint.initialize(x.size(), x.data(), y.size(), y.data(),
                    f.size(), f.data());
    bint.interpolate(xq.size(), xq.data(), yq.data());
    auto fq = bint.getInterpolatedFunction();    
    double errorMax = 0;
    for (int i=0; i<static_cast<int> (fqRef.size()); ++i)
    {
        errorMax = std::max(errorMax, std::abs(fqRef[i] - fq[i]));
    }
    EXPECT_NEAR(errorMax, 0, 1.e-14);
    //std::cout << errorMax << std::endl;

    std::pair<double, double> xLimits({x0, x1});
    std::pair<double, double> yLimits({y0, y1});
    bint.initialize(xLimits, yLimits, nx, ny, f.data());
    bint.interpolate(xq.size(), xq.data(), yq.data());
    errorMax = 0;
    for (int i=0; i<static_cast<int> (fqRef.size()); ++i)
    {
        errorMax = std::max(errorMax, std::abs(fqRef[i] - fq[i]));
    }
    EXPECT_NEAR(errorMax, 0, 1.e-14);
    //std::cout << errorMax << std::endl;
}

}
