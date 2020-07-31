#ifndef HYPOSTATION_HPP
#define HYPOSTATION_HPP
#include <vector>
#include <string>
#include <fstream>
/*!
 * @brief Reads the station file (with static corrections) used by HypoDD.
 * @note It is assumed the first correction is the P correction in seconds
 *       and the second correction is the S correction in seconds.
 */
struct HypoStation
{
    /// Reads the hypo station file
    void read(const std::string &fileName)
    {
        clear();
        auto ifl = std::ifstream(fileName);
        if (!ifl)
        {
            throw std::invalid_argument("file = "
                                      + fileName + " does not exist");
        }
        std::string line;
        std::vector<std::string> split;
        while (std::getline(ifl, line))
        {
            //ALT   UU  ENE  40 35.4168 111 38.2500 2635 .0  P -0.22  0.00  0.00  0.00 0 0.00 01
            boost::split(split, line, boost::is_any_of(" \n\t"),
                         boost::token_compress_on);
            // Have lines that effectively repeat - differ only by channel
            bool lskip = false;
            for (int i=0; i<static_cast<int> (stations.size()); ++i)
            {
                if (split[0] == stations[i] &&
                    split[1] == networks[i])
                {
                    lskip = true;
                    break;
                }
            }
            if (lskip){continue;}
            stations.push_back(split[0]);
            networks.push_back(split[1]);
            double lat = std::stod(split[3])
                       + std::stod(split[4])/60.;
            double lon = std::stod(split[5])
                       + std::stod(split[6])/60.;
            latitudes.push_back(lat);
            longitudes.push_back(-lon); // Works positive west
            elevations.push_back(std::stod(split[7]));
            pCorrections.push_back(std::stod(split[10]));
            sCorrections.push_back(std::stod(split[11]));
        }
        ifl.close();
    }
    /// Gets the station elevation
    void getPosition(const std::string &network,
                     const std::string &station,
                     double *latitude, double *longitude,
                     double *elevation) const
    {
        auto ncorr = static_cast<int> (networks.size());
        for (int i=0; i<ncorr; ++i)
        {
            if (network == networks[i] && station == stations[i])
            {
                *latitude = latitudes[i];
                *longitude = longitudes[i];
                *elevation = elevations[i];
                return; 
            }
        }
        throw std::runtime_error("cant find " + network + "." + station);
    }
    /// Gets the static correction for the network/station/phase
    double getCorrection(const std::string &network,
                         const std::string &station,
                         const std::string &phase) const noexcept
    {
        if (phase == "P")
        {
            return getPCorrection(network, station);
        }
        return getSCorrection(network, station);
    }
    /// Gets the P static correction for the network/station
    double getPCorrection(const std::string &network,
                          const std::string &station) const noexcept
    {
        auto ncorr = static_cast<int> (networks.size());
        for (int i=0; i<ncorr; ++i)
        {
            if (network == networks[i] && station == stations[i])
            {
                return pCorrections[i];
            }
        }
        std::cerr << "Couldnt find: " << network << "." << station <<std::endl;
        return 0;
    }
    /// Gets the S static correction for the network/station
    double getSCorrection(const std::string &network,
                          const std::string &station) const noexcept
    {
        auto ncorr = static_cast<int> (networks.size());
        for (int i=0; i<ncorr; ++i)
        {
            if (network == networks[i] && station == stations[i])
            {
                return sCorrections[i];
            }
        }
        std::cerr << "Couldnt find: " << network << "." << station << std::endl;
        return 0;
    }
    /// Clears the class
    void clear() noexcept
    {
        networks.clear();
        stations.clear();
        pCorrections.clear();
        sCorrections.clear();
        latitudes.clear();
        longitudes.clear();
        elevations.clear();
    }
    std::vector<std::string> networks;
    std::vector<std::string> stations;
    std::vector<double> pCorrections;
    std::vector<double> sCorrections; 
    std::vector<double> latitudes;
    std::vector<double> longitudes;
    std::vector<double> elevations;
};
#endif
