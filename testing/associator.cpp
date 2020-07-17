#include <string>
#include <cmath>
#include <vector>
#include "utilities.hpp"
#include "massociate/associatorParameters.hpp"
#include "massociate/associator.hpp"
#include "massociate/mesh/cartesian/points3d.hpp"
#include "massociate/arrival.hpp"
#include "massociate/pick.hpp"
#include "massociate/waveformIdentifier.hpp"
#include <gtest/gtest.h>

namespace
{

TEST(MAssociate, AssociatorParameters)
{
    int minArrivals = 8;
    int dbscanClusterSize = 12;
    double dbscanEpsilon = 0.33;
    int pageRankIterations = 25;
    double pageRankDamping = 0.6;
    int nTables= 8;
    MAssociate::AssociatorParameters parameters;

    EXPECT_NO_THROW(
        parameters.setMinimumNumberOfArrivalsToNucleate(minArrivals));
    EXPECT_EQ(parameters.getMinimumNumberOfArrivalsToNucleate(), minArrivals);

    EXPECT_NO_THROW(parameters.setNumberOfTravelTimeTables(nTables));
    EXPECT_EQ(parameters.getNumberOfTravelTimeTables(), nTables);

    EXPECT_NO_THROW(parameters.setDBSCANEpsilon(dbscanEpsilon));
    EXPECT_NEAR(parameters.getDBSCANEpsilon(), dbscanEpsilon, 1.e-14); 

    EXPECT_NO_THROW(parameters.setDBSCANMinimumClusterSize(dbscanClusterSize));
    EXPECT_EQ(parameters.getDBSCANMinimumClusterSize(), dbscanClusterSize);

    EXPECT_NO_THROW(parameters.setPageRankDampingFactor(pageRankDamping));
    EXPECT_NEAR(parameters.getPageRankDampingFactor(), pageRankDamping, 1.e-14);

    EXPECT_NO_THROW(
        parameters.setPageRankNumberOfIterations(pageRankIterations));
    EXPECT_EQ(parameters.getPageRankNumberOfIterations(), pageRankIterations);

    MAssociate::AssociatorParameters pCopy(parameters);
    EXPECT_EQ(pCopy.getMinimumNumberOfArrivalsToNucleate(), minArrivals);
    EXPECT_EQ(pCopy.getNumberOfTravelTimeTables(), nTables);
    EXPECT_NEAR(pCopy.getDBSCANEpsilon(), dbscanEpsilon, 1.e-14);
    EXPECT_EQ(pCopy.getDBSCANMinimumClusterSize(), dbscanClusterSize);
    EXPECT_NEAR(pCopy.getPageRankDampingFactor(), pageRankDamping, 1.e-14);
    EXPECT_EQ(pCopy.getPageRankNumberOfIterations(), pageRankIterations);

}

TEST(MAssociate, Associator)
{
    int minArrivals = 8;
    int nrec = 8;
    double x0 = 0;
    double y0 = 0;
    double z0 = 0;
    double x1 = 50*1.e3;
    double y1 = 50*1.e3;
    double z1 = 15*1.e3;
    double dx = 1000;
    double dy = 1000;
    double dz = 1000;
    auto nx = static_cast<int> ( std::round((x1 - x0)/dx) ) + 1;
    auto ny = static_cast<int> ( std::round((y1 - y0)/dy) ) + 1;
    auto nz = static_cast<int> ( std::round((z1 - z0)/dz) ) + 1;
    double vp = 4000;
    double vs = 2400;
    MAssociate::WaveformIdentifier sncl;
    sncl.setNetwork("FK");
    sncl.setChannel("HHZ");
    sncl.setLocationCode("00");
    std::vector<double> srcX( { x0 +   (dx*nx)/2,
                                x0 + 4*(dx*nx)/5,
                                x0 +   (dx*nx)/2} );
    std::vector<double> srcY( { y0 +   3*(dy*ny)/5, 
                                y0 +     (dy*ny)/3,
                                y0 +   3*(dy*ny)/5} );
    std::vector<double> srcZ( { z0 +   (dz*nz)/2,
                                z0 + 3*(dz*nz)/5,
                                z0 +   (dz*nz)/2} );
    std::vector<double> originTimes({5, 5, 10});
    std::vector<int> iSrcX(3), iSrcY(3), iSrcZ(3), iSrc(3);
    for (int i=0; i<3; ++i)
    {
        iSrcX[i] = static_cast<int> ( std::round((srcX[i] - x0)/dx) );
        iSrcY[i] = static_cast<int> ( std::round((srcY[i] - y0)/dy) );
        iSrcZ[i] = static_cast<int> ( std::round((srcZ[i] - z0)/dz) );
        iSrc[i] = iSrcZ[i]*nx*ny + iSrcY[i]*nx + iSrcX[i];
        srcX[i] = x0 + dx*iSrcX[i];
        srcY[i] = y0 + dy*iSrcY[i];
        srcZ[i] = z0 + dz*iSrcZ[i];
    }
    // Randomly put some receivers at the surface
    auto iRecX = generateUniformRandomNumbers(nrec, 0, nx - 1, 6493);
    auto iRecY = generateUniformRandomNumbers(nrec, 0, ny - 1, 3475);
    std::vector<int> iRecZ(iRecX.size(), 0); 
    // Create the geometry
    std::vector<double> xPoints, yPoints, zPoints;
    makePoints(x0, y0, z0, dx, dy, dz, nx, ny, nz,
               &xPoints, &yPoints, &zPoints);
    MAssociate::Mesh::Cartesian::Points3D<float> geometry;
    geometry.setNumberOfPoints(xPoints.size());
    geometry.setXPositions(xPoints.size(), xPoints.data());
    geometry.setYPositions(yPoints.size(), yPoints.data());
    geometry.setZPositions(zPoints.size(), zPoints.data());   
    // Initialize association engine
    MAssociate::AssociatorParameters parameters;
    auto nTables = 2*nrec;
    MAssociate::Associator<float> associator;
    parameters.setNumberOfTravelTimeTables(nTables);
    associator.initialize(parameters, geometry);
    // Create travel time tables
    for (int i=0; i<nrec; ++i)
    {
        auto xr = x0 + iRecX[i]*dx;
        auto yr = y0 + iRecY[i]*dy;
        auto zr = z0 + iRecZ[i]*dz;
        auto pTable = makeTravelTimeTable(xr, yr, zr, dx, dy, dz, nx, ny, nz,
                                          vp, x0, y0, z0);
        auto sTable = makeTravelTimeTable(xr, yr, zr, dx, dy, dz, nx, ny, nz,
                                          vs, x0, y0, z0);
    } 
    // Create a list of picks
    int nPicks = 0;
    std::vector<MAssociate::Pick> picks;
    for (int isrc=0; isrc<static_cast<int> (iSrc.size()); ++isrc)
    {
        for (int i=0; i<nrec; ++i)
        {
            std::string station = "T" + std::to_string(i + 1);
            sncl.setStation(station);
            auto xr = x0 + iRecX[i]*dx;
            auto yr = y0 + iRecY[i]*dy;
            auto zr = z0 + iRecZ[i]*dz;
            double pTime = computeTravelTime(xr, yr, zr,
                                             srcX[isrc], srcY[isrc], srcZ[isrc],
                                             vp) + originTimes[isrc];
            double sTime = computeTravelTime(xr, yr, zr,
                                             srcX[isrc], srcY[isrc], srcZ[isrc],
                                             vs) + originTimes[isrc];
            MAssociate::Pick pick;
            pick.setIdentifier(nPicks);
            pick.setWaveformIdentifier(sncl);
            pick.setTime(pTime);
            pick.setPhaseName("P");
            pick.setStandardDeviation(0.2/std::sqrt(12.));
            picks.push_back(pick);
            nPicks = nPicks + 1;
            if (i%2 == isrc%2)
            {
                pick.setIdentifier(nPicks);
                pick.setTime(sTime);
                pick.setPhaseName("S");
                pick.setStandardDeviation(0.4/std::sqrt(12.));
                picks.push_back(pick);
                nPicks = nPicks + 1;
            }
        }
    }
    // Add picks
 
/*
    std::vector<double> xSources({1000, 5000, 5000});
    std::vector<double> ySources({2000, 6000, 6000});
    std::vector<double> zSources({5000
*/
}

}
