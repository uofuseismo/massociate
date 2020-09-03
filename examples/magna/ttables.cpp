#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <iostream>
#include <vector>
#include <exception>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Constants.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include "bilinearInterpolation.hpp"
#include "hypoStation.hpp"
#include "h5io.hpp"
#include "blackList.hpp"

struct TravelTimeInterpolator
{
    double getTime(const double offset, const double depth)
    {
        bilin.interpolate(1, &offset, &depth);
        auto ptr = bilin.getInterpolatedFunctionPointer();
        return ptr[0];
    }
    std::vector<double> getTimes(const std::vector<double> &offsets,
                                 const std::vector<double> &depths)
    {
        if (offsets.size() != depths.size())
        {
            throw std::invalid_argument("offsets.size() != depths.size()");
        }
        bilin.interpolate(offsets.size(), offsets.data(), depths.data());
        return bilin.getInterpolatedFunction();
    }
    BilinearInterpolation<double> bilin;
    std::string phaseLabel;
    std::vector<double> travelTimes;
    std::vector<double> depths;
    std::vector<double> offsets;
};

struct Receiver
{
    std::string network;
    std::string station;
    double latitude = 0;
    double longitude = 0;
    double depth = 0;
};

bool operator==(const Receiver &a, const Receiver &b)
{
    return (a.network == b.network && a.station == b.station);
}

std::vector<Receiver> getReceiverListFromCSV(const std::string &fileName,
                                             const bool haveHeader = true)
{
    std::ifstream infl(fileName);
    if (!infl){throw std::runtime_error("Could not open: " + fileName);}
    std::vector<Receiver> receivers;
    std::string line;
    std::vector<std::string> split;
    int lineNumber = 0;
    while (std::getline(infl, line)) 
    {
        lineNumber = lineNumber + 1;
        if (lineNumber == 1 && haveHeader){continue;}
        //UU,ALT,ENE,ENN,ENZ,01,40.59028,-111.6375,2635.0
        boost::split(split, line, boost::is_any_of(",\n\t"),
                     boost::token_compress_on);
        Receiver receiver;
        receiver.network = split[0];
        receiver.station = split[1];
        receiver.latitude = std::stod(split[6]);
        receiver.longitude = std::stod(split[7]);
        receiver.depth = std::stod(split[8]);
        bool lnew = true;
        for (const auto &r : receivers)
        {
            if (receiver == r)
            {
                lnew = false;
                break;
            }
        }
        if (lnew){receivers.push_back(receiver);} 
    }
    infl.close();
    std::cout << "Number of receivers: " << receivers.size() << std::endl; 
    return receivers;
}

std::vector<Receiver> getReceiverListFromNodalFile(const std::string &fileName)
{
    std::ifstream infl(fileName);
    if (!infl){throw std::runtime_error("Could not open: " + fileName);}
    std::vector<Receiver> receivers;
    std::vector<std::string> split;
    std::string line;
    while (std::getline(infl, line))
    {
        // -112.02516440 40.75396926 1272.3 001 001 5462 UU001
        boost::split(split, line, boost::is_any_of(" ,\n\t"),
                     boost::token_compress_on);
        Receiver receiver;
        receiver.network = "UU";
        receiver.station = split[3];
        receiver.latitude = std::stod(split[1]);
        receiver.longitude = std::stod(split[0]);
        receiver.depth = std::round(std::stod(split[2]));
//std::cout << receiver.network << "," << receiver.station << "," << receiver.latitude << "," << receiver.longitude << "," << receiver.depth << std::endl;
        if (isBlackListed(receiver.network, receiver.station)){continue;}
        bool lnew = true;
        for (const auto &r : receivers)
        {
            if (receiver == r)
            {
                lnew = false;
                break;
            }
        }
        if (lnew){receivers.push_back(receiver);}
    }
    infl.close();
    std::cout << "Number of nodes: " << receivers.size() << std::endl;
    return receivers;
}

/// @brief Creates a travel time table interpolation structure from the
///        travel time tables computed by GrowClust. 
/// @param[in] fileName   The file name with the travel time table.
/// @param[in] phaseLabel The phase label (e.g., P or S).
/// @param[out] interp    This is the travel time with a travel time
///                       bilinear interpolation method.
void createGrowClustTravelTimeTable(const std::string &fileName,
                                    const std::string &phaseLabel,
                                    TravelTimeInterpolator *interp)
{
    std::ifstream infl(fileName);
    if (!infl)
    {
        throw std::runtime_error("Could not open: " + fileName);
    }
    int lineNumber = 0;
    std::string line;
    std::vector<std::string> split;
    int nDepths = 0;
    int nDistances = 0;
    std::vector<double> depths;
    std::vector<double> offsets;
    std::vector<double> travelTimes;
    while (std::getline(infl, line))
    {
        lineNumber = lineNumber + 1;
        if (lineNumber == 1){continue;} // Header
        // Remove leading blanks and parse
        boost::algorithm::trim_left(line);
        boost::split(split, line, boost::is_any_of("\n\t "),
                     boost::token_compress_on);
        // Read sizes
        if (lineNumber == 2)
        {
            nDistances = std::stoi(split[0]);
            nDepths = std::stoi(split[1]);
            offsets.reserve(nDistances);
            depths.resize(nDepths);
            travelTimes.resize(nDistances*nDepths, -1);
            continue;
        }
        // Read depths
        if (lineNumber == 3)
        {
            for (int id=0; id<nDepths; ++id)
            {
                depths[id] = std::stod(split[id]);
            }
            continue;
       }
       offsets.push_back(std::stod(split[0]));
       if (static_cast<int> (split.size() - 1) != nDepths)
       {
           throw std::runtime_error("Invalid size");
       }
       for (int i=0; i<nDepths; ++i)
       {
           auto indx = i*nDistances + offsets.size() - 1;
           travelTimes.at(indx) = std::stod(split[i+1]);
       }
    }
    infl.close();
    auto [tmin, tmax]
        = std::minmax_element(travelTimes.begin(), travelTimes.end());
    if (*tmin < 0){throw std::runtime_error("Failed to unpack ttimes");}
    // Create linearly interpolation
    interp->phaseLabel = phaseLabel;
    interp->offsets = offsets;
    interp->depths = depths;
    interp->travelTimes = travelTimes;
    interp->bilin.initialize(offsets.size(), offsets.data(),
                             depths.size(), depths.data(),
                             travelTimes.size(), travelTimes.data());
    // Check travel time interpolator works:
    //std::cout << travelTimes[10*nDistances+12] << " "  << interp->getTime(offsets[12], depths[10]) << std::endl;
}

/// @brief Create a grid in spherical coordinates.
void createGrid(const double lat0, const double lat1, const int nLat,
                const double lon0, const double lon1, const int nLon,
                const double z0,   const double z1,   const int nDep,
                std::vector<double> *lons,
                std::vector<double> *lats,
                std::vector<double> *depths)
{
    if (nDep < 2 || nLon < 2 || nLat < 2)
    {
        throw std::invalid_argument("need more than 2 grid points in x, y, z");
    }
    int nPts = nLat*nLon*nDep;
    std::cout << "Number of points in grid: " << nPts << std::endl;
    auto dLat = (lat1 - lat0)/(nLat - 1);
    auto dLon = (lon1 - lon0)/(nLon - 1);
    auto dDep = (z1 - z0)/(nDep - 1);
    lons->resize(nPts, 0);
    lats->resize(nPts, 0);
    depths->resize(nPts, -10);
    auto lonPtr = lons->data();
    auto latPtr = lats->data();
    auto depPtr = depths->data();
    // x
    #pragma omp simd collapse(3)
    for (int iLon=0; iLon<nLon; ++iLon)
    {
        // y 
        for (int iLat=0; iLat<nLat; ++iLat)
        {
            // z
            for (int id=0; id<nDep; ++id) 
            {
                auto indx = iLon*nDep*nLat + iLat*nDep + id; 
                lonPtr[indx] = lon0 + iLon*dLon;
                latPtr[indx] = lat0 + iLat*dLat;
                depPtr[indx] = z0   + id*dDep;
            }
        }
    }
    auto dmin = std::min_element(depths->begin(), depths->end());
    if (*dmin <=-10){throw std::runtime_error("algorithm failure");}
}

/// @brief Creates the candidate source points in the grid search.
/// @param[out] lons     The longitudes in degrees.
/// @param[out] lats     The latitudes in degrees.
/// @param[out] depths   The depths in kilometers.
void createSourcePoints(std::vector<double> *lats,
                        std::vector<double> *lons,
                        std::vector<double> *depths)
{
    double kmPerDeg = 111.195;
    double deltaCoarse = 4; // 5 km spacing
    double lat0Coarse = 40.5;
    double lat1Coarse = 41.0;
    double lon0Coarse =-112.24;
    double lon1Coarse =-111.76;
    double z0Coarse = 1.5; // 1.5 km is sea-level
    double z1Coarse = 21.5;
    int nLatCoarse = static_cast<int> (std::round( (lat1Coarse - lat0Coarse)
                                                  /deltaCoarse*kmPerDeg) );
    int nLonCoarse = static_cast<int> (std::round( (lon1Coarse - lon0Coarse)
                                                  /deltaCoarse*kmPerDeg) );
    int nDepCoarse = static_cast<int> (std::round( (z1Coarse - z0Coarse)
                                                  /deltaCoarse ));
 
    double deltaFine = 1.25; // 1 km spacing
    double lon0Fine =-112.2;
    double lon1Fine =-111.8;
    double lat0Fine = 40.6;
    double lat1Fine = 40.85;
    double z0Fine = 0;
    double z1Fine = 18;
    int nLatFine = static_cast<int> (std::round( (lat1Fine - lat0Fine)
                                                /deltaFine*kmPerDeg) );
    int nLonFine = static_cast<int> (std::round( (lon1Fine - lon0Fine)
                                                /deltaFine*kmPerDeg) );
    int nDepFine = static_cast<int> (std::round( (z1Fine - z0Fine)
                                                /deltaFine ));

    std::vector<double> latsCoarse, lonsCoarse, depthsCoarse;
    createGrid(lat0Coarse, lat1Coarse, nLatCoarse,
               lon0Coarse, lon1Coarse, nLonCoarse,
               z0Coarse, z1Coarse, nDepCoarse,
               &lonsCoarse, &latsCoarse, &depthsCoarse);
    // Add in Bingham
    lonsCoarse.push_back(-112.1453);
    latsCoarse.push_back(40.5162);
    depthsCoarse.push_back(0);
    //std::cout << depthsCoarse.size() << std::endl;
    // Create the fine grid
    std::vector<double> latsFine, lonsFine, depthsFine;
    createGrid(lat0Fine, lat1Fine, nLatFine,
               lon0Fine, lon1Fine, nLonFine,
               z0Fine,   z1Fine, nDepFine,
               &lonsFine, &latsFine, &depthsFine);
    lonsCoarse.insert(std::end(lonsCoarse),
                      std::begin(lonsFine), std::end(lonsFine));
    latsCoarse.insert(std::end(latsCoarse),
                      std::begin(latsFine), std::end(latsFine));
    depthsCoarse.insert(std::end(depthsCoarse),
                        std::begin(depthsFine), std::end(depthsFine));
    //std::cout << depthsCoarse.size() << std::endl;
    *lons = lonsCoarse;
    *lats = latsCoarse;
    *depths = depthsCoarse; 
}

// @brief Compute the source-to-receiver distances.
/// @param[in] receiverLatitude   The receiver's latitude in degrees.
/// @param[in] receiverLongitude  The receiver's longitude in degrees.
/// @param[in] sourceLatitudes    The source latitudes in degrees.
/// @param[in] sourceLongitudes   The source longitudes in degrees.
/// @param[out] distances         The source receiver distances in kilometers.
void computeDistances(const double receiverLatitude,
                      const double receiverLongitude,
                      const std::vector<double> &sourceLatitudes,
                      const std::vector<double> &sourceLongitudes,
                      std::vector<double> *distances)
{
    GeographicLib::Geodesic geodesic{GeographicLib::Constants::WGS84_a(),
                                     GeographicLib::Constants::WGS84_f()};
    distances->resize(sourceLatitudes.size());
    auto dPtr = distances->data();
    for (int i=0; i<static_cast<int> (sourceLatitudes.size()); ++i)
    {
        geodesic.Inverse(sourceLatitudes[i], sourceLongitudes[i],
                         receiverLatitude,   receiverLongitude,
                         dPtr[i]);
        dPtr[i] = dPtr[i]*1.e-3; // km
    }
}


int main()
{
    struct TravelTimeInterpolator pTimes;
    struct TravelTimeInterpolator sTimes;
/*
    auto receivers = getReceiverListFromCSV("../magna/magna_3c_stations.csv", true);
    auto nodes = getReceiverListFromNodalFile("../magna/Magna-Locs_2020-SN.txt");
    receivers.insert(std::end(receivers), std::begin(nodes), std::end(nodes));
*/
    HypoStation stations;
    stations.read("../magna/locate/magna.sta");
    std::cout << "Total number of receivers: " << stations.networks.size() << std::endl;
//    std::cout << "Total number of receivers: " << receivers.size() << std::endl;
    createGrowClustTravelTimeTable("../magna/TT.wasatch.pg", "P", &pTimes);
    createGrowClustTravelTimeTable("../magna/TT.wasatch.sg", "S", &sTimes);
    std::vector<double> lons, lats, depths, distances;
    createSourcePoints(&lats, &lons, &depths);
    H5IO h5io;
    h5io.openFileForWriting("../magna/travelTimeTables.h5");
    h5io.setGeometry(lats.size(), lats.data(), lons.data(), depths.data());
std::ofstream locFile("locations.txt");
std::ofstream staFile("stations.txt");

    for (int irec=0; irec<static_cast<int> (stations.networks.size()); ++irec)
    {
        computeDistances(stations.latitudes[irec], stations.longitudes[irec],
                         lats, lons, &distances);
        auto dmax = std::max_element(distances.begin(), distances.end());
        auto pTravelTimes = pTimes.getTimes(distances, depths);
        auto sTravelTimes = sTimes.getTimes(distances, depths);
        auto pmax = std::max_element(pTravelTimes.begin(), pTravelTimes.end());
        auto smax = std::max_element(sTravelTimes.begin(), sTravelTimes.end());
        if (irec == 0)
        {
           for (int j=0; j<static_cast<int> (lats.size()); ++j)
           {
               locFile << lons[j] << " " << lats[j] << " " << depths[j] << " " << distances[j]
                       << " " << pTravelTimes[j] << std::endl;
           }
        }
        staFile << stations.longitudes[irec] << " " << stations.latitudes[irec] << std::endl;
       std::cout << stations.stations[irec] << " " << stations.latitudes[irec] << " " 
                 << stations.longitudes[irec]
                 << " " << *dmax << " " << *pmax << " " << *smax << std::endl;
        h5io.addTravelTimeTable(stations.networks[irec],
                                stations.stations[irec],
                                "P",
                                pTravelTimes.size(), pTravelTimes.data());
        h5io.addTravelTimeTable(stations.networks[irec],
                                stations.stations[irec],
                                "S",
                                sTravelTimes.size(), sTravelTimes.data());
    }
/*
    for (int irec=0; irec<static_cast<int> (receivers.size()); ++irec)
    {
        computeDistances(receivers[irec].latitude, receivers[irec].longitude,
                         lats, lons, &distances); 
        auto dmax = std::max_element(distances.begin(), distances.end());
        auto pTravelTimes = pTimes.getTimes(distances, depths); 
        auto sTravelTimes = sTimes.getTimes(distances, depths);
        auto pmax = std::max_element(pTravelTimes.begin(), pTravelTimes.end());
        auto smax = std::max_element(sTravelTimes.begin(), sTravelTimes.end());
if (irec == 0)
{
for (int j=0; j<lats.size(); ++j)
{
   locFile << lons[j] << " " << lats[j] << " " << distances[j] << " " << pTravelTimes[j] << std::endl;
}
}
staFile << receivers[irec].longitude << " " << receivers[irec].latitude << std::endl;
std::cout << receivers[irec].station << " " << receivers[irec].latitude << " " << receivers[irec].longitude << " " << *dmax << " " << *pmax << " " << *smax << std::endl;
        h5io.addTravelTimeTable(receivers[irec].network,
                                receivers[irec].station,
                                "P",
                                pTravelTimes.size(), pTravelTimes.data());
        h5io.addTravelTimeTable(receivers[irec].network,
                                receivers[irec].station,
                                "S",
                                sTravelTimes.size(), sTravelTimes.data());
    }
*/
} 
