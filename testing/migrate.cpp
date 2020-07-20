#include <vector>
#include <cmath>
#include "massociate/migrationParameters.hpp"
#include "massociate/migrate.hpp"
#include "massociate/pick.hpp"
#include "massociate/waveformIdentifier.hpp"
#include "private/analyticSignals.hpp"
#include "utilities.hpp"
#include <gtest/gtest.h>

namespace
{
/// @brief Simple direct correlation.
/// @param[in] g1   The first signal.
/// @param[in] g2   The second signal.
/// @result \f$ g1 \star g2 \f$ where the convention is the max of the 
///         correlation will shift g1 so that it best aligns with g2.
/// @throws std::invalid_argument if g1.size() != g2.size().
template<typename T>
std::vector<T> simpleDirectCorrelation(const std::vector<T> &g1,
                                       const std::vector<T> &g2,
                                       double dt = 1);

TEST(MAssociate, MigrationParameters)
{
    MAssociate::MigrationParameters parameters;
    int nTables = 5;
    int nPoints = 110;
    int tileSize = 90;
    int dbscanClusterSize = 5;
    double dbscanEpsilon = 0.3;
    auto function = MAssociate::AnalyticCorrelationFunction::BOXCAR;

    EXPECT_NO_THROW(parameters.setNumberOfPointsInTravelTimeTable(nPoints));
    EXPECT_EQ(parameters.getNumberOfPointsInTravelTimeTable(), nPoints);

    EXPECT_NO_THROW(parameters.setNumberOfTravelTimeTables(nTables));
    EXPECT_EQ(parameters.getNumberOfTravelTimeTables(), nTables);

    EXPECT_NO_THROW(parameters.setTileSize(tileSize));
    EXPECT_EQ(parameters.getTileSize(), tileSize);
    EXPECT_NO_THROW(parameters.setTileSize(nPoints + 1));
    EXPECT_EQ(parameters.getTileSize(), nPoints);

    parameters.setAnalyticCorrelationFunction(function);
    EXPECT_EQ(parameters.getAnalyticCorrelationFunction(), function);
    EXPECT_NO_THROW(parameters.setDBSCANEpsilon(dbscanEpsilon));
    EXPECT_NEAR(parameters.getDBSCANEpsilon(), dbscanEpsilon, 1.e-12);
    EXPECT_NO_THROW(parameters.setDBSCANMinimumClusterSize(dbscanClusterSize));
    EXPECT_EQ(parameters.getDBSCANMinimumClusterSize(), dbscanClusterSize);

    MAssociate::MigrationParameters pCopy(parameters);
    EXPECT_EQ(pCopy.getNumberOfPointsInTravelTimeTable(), nPoints);
    EXPECT_EQ(pCopy.getNumberOfTravelTimeTables(), nTables);
    EXPECT_EQ(pCopy.getTileSize(), nPoints);
    EXPECT_EQ(pCopy.getAnalyticCorrelationFunction(), function);
    EXPECT_NEAR(pCopy.getDBSCANEpsilon(), dbscanEpsilon, 1.e-12);
    EXPECT_EQ(pCopy.getDBSCANMinimumClusterSize(), dbscanClusterSize);

    parameters.clear();
}

TEST(MAssociate, AnalyticGaussianCorrelation)
{
    double tmin = 0;
    double tmax = 5;
    double dt = 0.01;
    double p1 = 1.3;  // P pick time
    double p2 = 2.1;  // S pick time
    double sd1 = 0.08; // P wave standard deviation
    double sd2 = 0.20; // S wave stadnard deviation
    auto npts = static_cast<int> ( std::round((tmax - tmin)/dt) ) + 1;
    // Create Gaussian signals to correlate
    std::vector<double> times(npts, 0);
    std::vector<double> g1(npts, 0);
    std::vector<double> g2(npts, 0);
    auto xnorm1 = 1./(sd1*sqrt(2*M_PI));
    auto xnorm2 = 1./(sd2*sqrt(2*M_PI));
    for (int i=0; i<npts; ++i)
    {
        times[i] = tmin + i*dt;
        g1[i] = xnorm1*exp( -(pow(times[i] - p1, 2))/(2*sd1*sd1) );
        g2[i] = xnorm2*exp( -(pow(times[i] - p2, 2))/(2*sd2*sd2) );
    }
    // Compute correlation times
    int xcLen = 2*npts - 1;
    auto xc = simpleDirectCorrelation(g1, g2);
    // Correlation is defined in samples but approximates an integral hence
    // we need a normalization.
    for (int i=0; i<xc.size(); ++i){xc[i] = dt*xc[i];}
    std::vector<double> xcAnalytic(xcLen, 0);
    std::vector<double> xcTimes(xcLen, 0);
    auto xnormAmp = analyticGaussianCorrelationAmplitudeNormalization(sd1, sd2);
    auto xnormExp = analyticGaussianCorrelationExponentNormalization(sd1, sd2);
//std::cout << xnormAmp << std::endl;
//std::cout << xnormExp << std::endl;
    for (int i=0; i<xcLen; ++i)
    {
        xcTimes[i] =-tmax + i*dt;
        xcAnalytic[i] = analyticGaussianCorrelation(xcTimes[i], p1, p2,
                                                    xnormAmp, xnormExp);
        EXPECT_NEAR(xcAnalytic[i], xc[i], 1.e-7);
        //std::cout << xcTimes[i] << "," << xc[i] << "," <<  xcAnalytic[i] << std::endl;
    }
    auto error = infinityNorm(xcLen, xc.data(), xcAnalytic.data());
    EXPECT_NEAR(error, 0, 1.e-7);
}

TEST(MAssociate, AnalyticBoxcar)
{
    double tmin = 0;
    double tmax = 5;
    double dt = 0.01;
    double p1 = 0.5;  // P pick time
    double p2 = 2.1;  // S pick time
    double w1 = 2*0.08; // P wave width
    double w2 = 2*0.20;// S wave width
    auto npts = static_cast<int> ( std::round((tmax - tmin)/dt) ) + 1;
    // P pick leads and S pick and have different widths
    std::vector<double> times(npts, 0);
    std::vector<double> b1(npts, 0);
    std::vector<double> b2(npts, 0);
    for (int i=0; i<npts; ++i)
    {
        times[i] = tmin + i*dt;
        if (times[i] >= p1 - w1/2 && times[i] < p1 + w1/2){b1[i] = 1/w1;}
        if (times[i] >= p2 - w2/2 && times[i] < p2 + w2/2){b2[i] = 1/w2;}
    }
    int xcLen = 2*npts - 1;
    auto xc = simpleDirectCorrelation(b1, b2, dt);
    std::vector<double> xcAnalytic(xcLen, 0);
double xsum =0;
    for (int i=0; i<xcLen; ++i)
    {
        double deltaT =-tmax + i*dt;
        xcAnalytic[i] = analyticBoxcarCorrelation(deltaT, p1, p2, w1, w2);
xsum = xsum + xcAnalytic[i]*dt;
/*
if (xc[i] > 1.e-8 || xcAnalytic[i] > 1.e-8)
{
   std::cout << times[i] << " " << xc[i] << " " << xcAnalytic[i] << " " << xc[i]/xcAnalytic[i] << std::endl;
}
*/
    }
//std::cout << xsum << std::endl;
    auto error = infinityNorm(xc.size(), xc.data(), xcAnalytic.data());
    EXPECT_NEAR(error, 0, 1.e-10);
    //std::cout << "Error: " << error << std::endl;
    // Put first pick after second pick same width
    p1 = 3.1;
    p2 = 2;
    w1 = 0.1;
    w2 = 0.1;
    std::fill(b1.begin(), b1.end(), 0);
    std::fill(b2.begin(), b2.end(), 0);
    for (int i=0; i<npts; ++i)
    {
        if (times[i] >= p1 - w1/2 && times[i] < p1 + w1/2){b1[i] = 1/w1;}
        if (times[i] >= p2 - w2/2 && times[i] < p2 + w2/2){b2[i] = 1/w2;}
    }
    xc = simpleDirectCorrelation(b1, b2, dt);
    for (int i=0; i<xcLen; ++i)
    {
        double deltaT =-tmax + i*dt;
        xcAnalytic[i] = analyticBoxcarCorrelation(deltaT, p1, p2, w1, w2); 
/*
if (xc[i] > 1.e-8 || xcAnalytic[i] > 1.e-8)
{
   std::cout << times[i] << " " << xc[i] << " " << xcAnalytic[i] << " " << xc[i]/xcAnalytic[i] << std::endl;
}
*/
    }
    error = infinityNorm(xc.size(), xc.data(), xcAnalytic.data());
    EXPECT_NEAR(error, 0, 1.e-10);
    //std::cout << "Error: " << error << std::endl;
    // Another thing is first pick after second pick but different widths
    p2 = 1.9;
    p1 = 1.6;
    w1 = 0.1;
    w2 = 0.2;
    std::fill(b1.begin(), b1.end(), 0);
    std::fill(b2.begin(), b2.end(), 0);
    for (int i=0; i<npts; ++i)
    {
        if (times[i] >= p1 - w1/2 && times[i] < p1 + w1/2){b1[i] = 1/w1;}
        if (times[i] >= p2 - w2/2 && times[i] < p2 + w2/2){b2[i] = 1/w2;}
    }
    xc = simpleDirectCorrelation(b1, b2, dt);
    for (int i=0; i<xcLen; ++i)
    {
        double deltaT =-tmax + i*dt;
        xcAnalytic[i] = analyticBoxcarCorrelation(deltaT, p1, p2, w1, w2);
    }
    error = infinityNorm(xc.size(), xc.data(), xcAnalytic.data());
    EXPECT_NEAR(error, 0, 1.e-10);
    //std::cout << "Error: " << error << std::endl;
}

TEST(MAssociate, Migrate)
{
    const double sqrt12 = std::sqrt(12.);
    MAssociate::Migrate<float> migrate;
    MAssociate::Migrate<double> migrateGauss;
    std::string network = "FK";
    std::string channel = "HHZ";
    std::string location = "01";
    MAssociate::WaveformIdentifier sncl;
    MAssociate::MigrationParameters parameters;
    sncl.setNetwork(network);
    sncl.setChannel(channel);
    sncl.setLocationCode(location);
    int nx = 41;
    int ny = 43;
    int nz = 10;
    double x0 = 0;
    double y0 = 0;
    double z0 = 0;
    double dx = 1000;
    double dy = 1000;
    double dz = 1000;
    double vp = 4000;
    double vs = 2200;
    double pStaticCorrection =-0.05;
    double sStaticCorrection = 0.05;
    int nrec = 11;
    int isrcx = 15;
    int isrcy = 37;
    int isrcz = 7;
    int isrc = isrcz*nx*ny + isrcy*nx + isrcx;
    auto xs = x0 + dx*isrcx;
    auto ys = y0 + dy*isrcy;
    auto zs = z0 + dz*isrcz;
    auto ot = 25;
    auto irecx = generateUniformRandomNumbers(nrec, 0, nx - 1, 4993);
    auto irecy = generateUniformRandomNumbers(nrec, 0, ny - 1, 4995);
    std::vector<int> irecz(nrec, 0); // Put all these at surface
    // Create the migration parameters
    parameters.setNumberOfTravelTimeTables(2*nrec);
    parameters.setNumberOfPointsInTravelTimeTable(nx*ny*nz);
    parameters.setTileSize(512);
    parameters.setAnalyticCorrelationFunction(
         MAssociate::AnalyticCorrelationFunction::BOXCAR);
    EXPECT_NO_THROW(migrate.initialize(parameters));
    // -> Gaussian
    parameters.setAnalyticCorrelationFunction(
         MAssociate::AnalyticCorrelationFunction::GAUSSIAN);
    EXPECT_NO_THROW(migrateGauss.initialize(parameters));
    // Create the travel time tables 
    for (int i=0; i<nrec; ++i)
    {
        EXPECT_FALSE(migrate.haveAllTravelTimeTables());
        std::string station = "T" + std::to_string(i + 1);
        auto xr = x0 + irecx[i]*dx;
        auto yr = y0 + irecy[i]*dy;
        auto zr = z0 + irecz[i]*dz;
        auto ptable = makeTravelTimeTable(xr, yr, zr, dx, dy, dz, nx, ny, nz,
                                          vp, x0, y0, z0);
        auto stable = makeTravelTimeTable(xr, yr, zr, dx, dy, dz, nx, ny, nz,
                                          vs, x0, y0, z0);
        migrate.setTravelTimeTable(network, station, "P",
                                   ptable.size(), ptable.data());
        migrate.setTravelTimeTable(network, station, "S",
                                   stable.size(), stable.data());
        migrateGauss.setTravelTimeTable(network, station, "P",
                                        ptable.size(), ptable.data());
        migrateGauss.setTravelTimeTable(network, station, "S",
                                        stable.size(), stable.data());
    }
    EXPECT_TRUE(migrate.haveAllTravelTimeTables());
    // Figure out max dt
    auto maxDT = migrate.getMaximumDifferentialTravelTime();
    float dtMaxRef = 0;
    for (int i=0; i<nrec; ++i)
    {
        std::string station = "T" + std::to_string(i + 1);
        std::vector<std::string> phases({"P", "S"});
        for (int ip=0; ip<2; ++ip)
        {
            auto t1 = migrate.getTravelTimeTable(network, station, phases[ip]);
            for (int j=i+1; j<nrec; ++j)
            {
                std::string station = "T" + std::to_string(j + 1);
                for (int jp=0; jp<2; ++jp)
                {
                    auto t2 = migrate.getTravelTimeTable(network, station,
                                                         phases[jp]);
                    #pragma omp simd reduction(min:dtMaxRef)
                    for (int k=0; k<static_cast<int> (t2.size()); ++k)
                    {
                        dtMaxRef = std::max(dtMaxRef, std::abs(t2[k] - t1[k]));
                    }
                }
            }
        }
    }
    EXPECT_NEAR(dtMaxRef, maxDT, 1.e-4);
    // Verify I can get the tables back
    // Create the travel time tables and a list of picks
    for (int i=0; i<nrec; ++i)
    {
        std::string station = "T" + std::to_string(i + 1);
        auto xr = x0 + irecx[i]*dx;
        auto yr = y0 + irecy[i]*dy;
        auto zr = z0 + irecz[i]*dz;
        auto pTableRef = makeTravelTimeTable(xr, yr, zr, dx, dy, dz, nx, ny, nz,
                                             vp, x0, y0, z0);
        auto sTableRef = makeTravelTimeTable(xr, yr, zr, dx, dy, dz, nx, ny, nz,
                                             vs, x0, y0, z0);
        // Get the tables
        auto pTable = migrate.getTravelTimeTable(network, station, "P");
        auto sTable = migrate.getTravelTimeTable(network, station, "S");
        EXPECT_EQ(pTableRef.size(), pTable.size());
        EXPECT_EQ(sTableRef.size(), sTable.size());
        EXPECT_LT(infinityNorm(pTable.size(), pTableRef.data(), pTable.data()),
                  1.e-6);
        EXPECT_LT(infinityNorm(sTable.size(), sTableRef.data(), sTable.data()), 
                  1.e-6);
    }
    // Add the picks
    std::vector<MAssociate::Pick> picks;
    int nPicks = 0;
    for (int i=0; i<nrec; ++i)
    {
        std::string station = "T" + std::to_string(i + 1);
        sncl.setStation(station);
        auto xr = x0 + irecx[i]*dx;
        auto yr = y0 + irecy[i]*dy;
        auto zr = z0 + irecz[i]*dz;
        double pTime = computeTravelTime(xr, yr, zr, xs, ys, zs, vp)
                     + ot + pStaticCorrection;
        double sTime = computeTravelTime(xr, yr, zr, xs, ys, zs, vs)
                     + ot + sStaticCorrection;
        MAssociate::Pick pick;
        pick.setIdentifier(nPicks);
        pick.setWaveformIdentifier(sncl);
        pick.setTime(pTime);
        pick.setPhaseName("P");
        pick.setStandardDeviation(0.2/sqrt12);
        pick.setStaticCorrection(pStaticCorrection);
        EXPECT_NO_THROW(migrate.addPick(pick));
        picks.push_back(pick);
        nPicks = nPicks + 1;        
        if (i%2 == 0)
        {
            pick.setIdentifier(nPicks);
            pick.setTime(sTime);
            pick.setPhaseName("S");
            pick.setStandardDeviation(0.4/sqrt12);
            pick.setStaticCorrection(sStaticCorrection);
            EXPECT_NO_THROW(migrate.addPick(pick));
            picks.push_back(pick);
            nPicks = nPicks + 1;
        }
    }
    EXPECT_EQ(migrate.getNumberOfPicks(), nPicks);
    // There is no reason to pop any picks because of causality in this
    // example so do a straight sum of the migration contribution for all pairs
    double sumGaussian = 0;
    double sumBoxcar = 0;
    for (int ia=0; ia<static_cast<int> (picks.size()); ++ia)
    {
        EXPECT_NO_THROW(migrateGauss.addPick(picks[ia]));
        for (int ja=ia+1; ja<static_cast<int> (picks.size()); ++ja)
        {
            double std1 = picks[ia].getStandardDeviation();
            double std2 = picks[ja].getStandardDeviation();
            // Perfect data
            double deltaT = 0;
            double p1 = 0;
            double p2 = 0;
            auto normalizeAmp
                = analyticGaussianCorrelationAmplitudeNormalization(std1, std2);
            auto normalizeExp
                = analyticGaussianCorrelationExponentNormalization(std1, std2);
            auto boxcarWeight1 = std::sqrt(12.)*std1;
            auto boxcarWeight2 = std::sqrt(12.)*std2;
            sumGaussian = sumGaussian
                        + analyticGaussianCorrelation(deltaT, p1, p2,
                                                      normalizeAmp,
                                                      normalizeExp);
            sumBoxcar = sumBoxcar
                      + analyticBoxcarCorrelation(deltaT, p1, p2,
                                                  boxcarWeight1, boxcarWeight2);
        }
    }
    // Migrate
    EXPECT_NO_THROW(migrate.migrate());
    auto imageMax = migrate.getImageMaximum();
    auto ttimesToMax = migrate.getTravelTimesToMaximum();
    EXPECT_EQ(imageMax.first, isrc);
    EXPECT_NEAR(imageMax.second, sumBoxcar, 1.e-3); // 477.5*1.e-6 is mach eps
    for (int ia=0; ia<static_cast<int> (picks.size()); ++ia)
    {
        auto waveid = picks[ia].getWaveformIdentifier();
        // Note that there is no static correction 
        auto tEst = migrate.getTravelTime(waveid.getNetwork(),
                                          waveid.getStation(),
                                          picks[ia].getPhaseName(),
                                          imageMax.first);
        tEst = tEst + picks[ia].getStaticCorrection(); // Fix the static correction
        auto tObs = picks[ia].getTime();
        EXPECT_NEAR(tObs - ot, tEst, 1.e-4);
        EXPECT_NEAR(ttimesToMax[ia], tObs - ot, 1.e-4);
    }
    // Migrate Gaussian
    EXPECT_NO_THROW(migrateGauss.migrate());
    imageMax = migrateGauss.getImageMaximum();
    EXPECT_EQ(imageMax.first, isrc);
    EXPECT_NEAR(imageMax.second, sumGaussian, 1.e-3); // 509.329*1.e-6 -> 1.e-3
    // Clear the picks
    migrate.clearPicks();
    EXPECT_EQ(migrate.getNumberOfPicks(), 0);
}


template<typename T>
std::vector<T> simpleDirectCorrelation(const std::vector<T> &g1,
                                       const std::vector<T> &g2,
                                       const double dt)
{
    if (g1.size() != g2.size())
    {
        throw std::invalid_argument("Only does g1.size() == g2.size()");
    }
    auto npts = static_cast<int> (g1.size());
    auto xclen = 2*npts - 1;
    // Cross-correlate manually - use same convention as IPP which asks:
    // What is the shift that moves g1 to best align to g2?
    std::vector<T> xc(xclen, 0);
    std::vector<T> awork(xclen, 0);
    for (int i=0; i<npts; i++)
    {
        auto k = npts - 1 - i;
        awork[i] = g2[k];
    }
    for (int i=1; i<=xclen; i++)
    {
        auto kmin = std::max(0, i-npts);
        auto kmax = std::min(i - 1, npts - 1);
        for (int k=kmin; k<=kmax; k++)
        {
            xc[xclen-i] = xc[xclen-i] + g1[k]*awork[i-k-1];
        }
    }
    if (dt != 1)
    {
        for (int i=0; i<static_cast<int> (xc.size()); ++i){xc[i] = dt*xc[i];}
    }
    return xc;
}

}
