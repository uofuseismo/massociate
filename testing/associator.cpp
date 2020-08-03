#include <string>
#include <cmath>
#include <vector>
#include "utilities.hpp"
#include "massociate/associatorParameters.hpp"
#include "massociate/associator.hpp"
#include "massociate/mesh/cartesian/points3d.hpp"
#include "massociate/arrival.hpp"
#include "massociate/event.hpp"
#include "massociate/pick.hpp"
#include "massociate/waveformIdentifier.hpp"
#include "private/weightedStatistics.hpp"
#include <gtest/gtest.h>

namespace
{

TEST(MAssociate, AssociatorParameters)
{
    int minArrivals = 8;
    int dbscanClusterSize = 12;
    double dbscanEpsilon = 0.33;
    int pageRankIterations = 25;
    int evidInterval = 3;
    double pageRankDamping = 0.6;
    auto otObjFn = MAssociate::OriginTimeObjectiveFunction::L1;
    int nTables = 8;
    int tileSize = 128;
    auto function = MAssociate::AnalyticCorrelationFunction::GAUSSIAN;
    MAssociate::AssociatorParameters parameters;
 
    EXPECT_EQ(parameters.getEventIdentifierInterval(), 1);

    EXPECT_NO_THROW(
        parameters.setMinimumNumberOfArrivalsToNucleate(minArrivals));
    EXPECT_EQ(parameters.getMinimumNumberOfArrivalsToNucleate(), minArrivals);

    EXPECT_NO_THROW(parameters.setNumberOfTravelTimeTables(nTables));
    EXPECT_EQ(parameters.getNumberOfTravelTimeTables(), nTables);
 
    EXPECT_NO_THROW(parameters.setTileSize(tileSize));
    EXPECT_EQ(parameters.getTileSize(), tileSize);

    parameters.setAnalyticCorrelationFunction(function);
    EXPECT_EQ(parameters.getAnalyticCorrelationFunction(), function);

    parameters.setOriginTimeObjectiveFunction(otObjFn);
    EXPECT_EQ(parameters.getOriginTimeObjectiveFunction(), otObjFn);

    EXPECT_NO_THROW(parameters.setDBSCANEpsilon(dbscanEpsilon));
    EXPECT_NEAR(parameters.getDBSCANEpsilon(), dbscanEpsilon, 1.e-14); 

    EXPECT_NO_THROW(parameters.setDBSCANMinimumClusterSize(dbscanClusterSize));
    EXPECT_EQ(parameters.getDBSCANMinimumClusterSize(), dbscanClusterSize);

    EXPECT_NO_THROW(parameters.setPageRankDampingFactor(pageRankDamping));
    EXPECT_NEAR(parameters.getPageRankDampingFactor(), pageRankDamping, 1.e-14);

    EXPECT_NO_THROW(
        parameters.setPageRankNumberOfIterations(pageRankIterations));
    EXPECT_EQ(parameters.getPageRankNumberOfIterations(), pageRankIterations);

    EXPECT_NO_THROW(parameters.setEventIdentifierInterval(evidInterval));
    EXPECT_EQ(parameters.getEventIdentifierInterval(), evidInterval);

    MAssociate::AssociatorParameters pCopy(parameters);
    EXPECT_EQ(pCopy.getMinimumNumberOfArrivalsToNucleate(), minArrivals);
    EXPECT_EQ(pCopy.getNumberOfTravelTimeTables(), nTables);
    EXPECT_EQ(pCopy.getAnalyticCorrelationFunction(), function);
    EXPECT_EQ(pCopy.getOriginTimeObjectiveFunction(), otObjFn);
    EXPECT_EQ(pCopy.getTileSize(), tileSize);
    EXPECT_NEAR(pCopy.getDBSCANEpsilon(), dbscanEpsilon, 1.e-14);
    EXPECT_EQ(pCopy.getDBSCANMinimumClusterSize(), dbscanClusterSize);
    EXPECT_NEAR(pCopy.getPageRankDampingFactor(), pageRankDamping, 1.e-14);
    EXPECT_EQ(pCopy.getPageRankNumberOfIterations(), pageRankIterations);
    EXPECT_EQ(pCopy.getEventIdentifierInterval(), evidInterval);
}

TEST(MAssociator, weightedStatistics)
{
    // https://rdrr.io/cran/spatstat/man/weighted.median.html
    // library(spatstat)
    // x <- c(1.1, 5.3, 3.7, 2.1, 7.0, 9.9)
    // y <- c(1.1, 0.4, 2.1, 3.5, 1.2, 0.8)
    // weighted.median(x, y)
    std::vector<double> x({1.1, 5.3, 3.7, 2.1, 7.0, 9.9});
    std::vector<double> wts({1.1, 0.4, 2.1, 3.5, 1.2, 0.8});
    auto wm = weightedMedian(x.size(), x.data(), wts.data());
    EXPECT_NEAR(wm, 2.085714, 1.e-6);
    // xpost <- c(0.1,0.35,0.05,0.1,0.15,0.05,0.2)
    // weighted.median(xpost, xpost) 
    std::vector<double> xpost({0.1,0.35,0.05,0.1,0.15,0.05,0.2});
    std::vector<double> wpost(xpost);
    wm = weightedMedian(xpost.size(), xpost.data(), wpost.data());
    EXPECT_NEAR(wm, 0.1625, 1.e-6);
    // Do an even length example
    wm = weightedMedian(xpost.size() - 1, xpost.data(), wpost.data());
    EXPECT_NEAR(wm, 0.1333333, 1.e-6);
    // And test an if statement conditional
    std::vector<double> x2({1, 2});
    std::vector<double> w2({1, 4});
    wm = weightedMedian(x2.size(), x2.data(), w2.data());
    EXPECT_NEAR(wm, 1.375, 1.e-6);
    // And an edge case of n=1
    EXPECT_NEAR(x[0], weightedMedian(1, x.data(), wts.data()), 1.e-6); 
    //------------------------------------------------------------------------//
    // Regular median
    wm = median(x.size(), x.data());
    EXPECT_NEAR(wm, 4.5, 1.e-6);
    wm = median(x.size()-1, x.data());
    EXPECT_NEAR(wm, 3.7, 1.e-6);
    EXPECT_NEAR(x[0], median(1, x.data()), 1.e-12);
    //------------------------------------------------------------------------//
    // Weighted mean and mean
    wm = weightedMean(xpost.size(), xpost.data(), xpost.data());
    EXPECT_NEAR(wm, 0.21, 1.e-6);
    EXPECT_NEAR(weightedMean(1, x.data(), x.data()), x[0], 1.e-12);
    wm = mean(xpost.size(), xpost.data());
    EXPECT_NEAR(wm, 0.1428571, 1.e-6);
    EXPECT_NEAR(mean(1, x.data()), x[0], 1.e-12); 
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
    parameters.setMinimumNumberOfArrivalsToNucleate(minArrivals);
    EXPECT_NO_THROW(associator.initialize(parameters, geometry));
    EXPECT_TRUE(associator.isInitialized());
    // Create travel time tables
    for (int i=0; i<nrec; ++i)
    {
        std::string station = "T" + std::to_string(i + 1);
        auto xr = x0 + iRecX[i]*dx;
        auto yr = y0 + iRecY[i]*dy;
        auto zr = z0 + iRecZ[i]*dz;
        auto pTable = makeTravelTimeTable(xr, yr, zr, dx, dy, dz, nx, ny, nz,
                                          vp, x0, y0, z0);
        auto sTable = makeTravelTimeTable(xr, yr, zr, dx, dy, dz, nx, ny, nz,
                                          vs, x0, y0, z0);
        associator.setTravelTimeTable(sncl.getNetwork(), station, "P",
                                      pTable.size(), pTable.data());
        associator.setTravelTimeTable(sncl.getNetwork(), station, "S",
                                      sTable.size(), sTable.data());
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
            pick.setIdentifier(nPicks + 100*isrc);
            pick.setWaveformIdentifier(sncl);
            pick.setTime(pTime);
            pick.setPhaseName("P");
            pick.setStandardDeviation(0.2/std::sqrt(12.));
            picks.push_back(pick);
            nPicks = nPicks + 1;
            if (i%2 == isrc%2)
            {
                pick.setIdentifier(nPicks + 100*isrc);
                pick.setTime(sTime);
                pick.setPhaseName("S");
                pick.setStandardDeviation(0.4/std::sqrt(12.));
                picks.push_back(pick);
                nPicks = nPicks + 1;
            }
        }
    }
    // Add picks
    for (const auto &pick : picks)
    {
        EXPECT_NO_THROW(associator.addPick(pick));
    }
    // Associate
    associator.associate(); 
/*
    std::vector<double> xSources({1000, 5000, 5000});
    std::vector<double> ySources({2000, 6000, 6000});
    std::vector<double> zSources({5000
*/
    // Check the associations
    EXPECT_EQ(associator.getNumberOfEvents(), 3);
    auto events = associator.getEvents();
    std::sort(events.begin(), events.end(), 
              [](const MAssociate::Event &a, const MAssociate::Event &b)
              {
                 return a.getOriginTime() < b.getOriginTime();
              });
    for (const auto &event : events)
    {
    }
    // Test clear
    associator.clearPicks();
    EXPECT_EQ(associator.getNumberOfPicks(), 0);
    associator.clearEvents();    
    EXPECT_EQ(associator.getNumberOfEvents(), 0);
}

}
