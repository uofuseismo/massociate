#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <massociate/associator.hpp>
#include <massociate/associatorParameters.hpp>
#include <massociate/mesh/spherical/points3d.hpp>
#include <massociate/waveformIdentifier.hpp>
#include <massociate/pick.hpp>
#include <massociate/arrival.hpp>
#include <massociate/event.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <sff/utilities/time.hpp>
#include "hypoStation.hpp"
#include "blackList.hpp"
#include "h5io.hpp"
#ifdef BUILD_MPI
#include <mpi.h>
#endif

/*
bool operator==(const Pick &a, const Pick &b)
{
    if (a.network   != b.network){return false;}
    if (a.station   != b.station){return false;}
    if (a.component != b.component){return false;}
    if (a.phase     != b.phase){return false;}
    if (std::abs(a.error - b.error) > 1.e-8){return false;}
    if (a.time.getEpochalTime() != b.time.getEpochalTime()){return false;}
    return true;
}
*/

bool parsePick(const std::string &line,
               MAssociate::Pick *pick,
               bool heuristicWeight = true,
               double pModelingError = 0.1,
               double sModelingError = 0.2)
{
    // http://alomax.free.fr/nlloc/index.html
    //CAPU   UU   ENZ  ? S      ? 20200317 1346  0.7750 GAU  1.00e-01 -1.00e+00 -1.00e+00 -1.00e+00
    //FTT    UU   ENZ  ? S      ? 20200317 1346  0.9950 GAU  1.00e-01 -1.00e+00 -1.00e+00 -1.00e+00
    //CTU    UU   HHZ  ? S      ? 20200317 1346  1.9100 GAU  1.00e-01 -1.00e+00 -1.00e+00 -1.00e+00
    //GMU    UU   EHZ  ? S      ? 20200317 1346  2.0750 GAU  1.00e-01 -1.00e+00 -1.00e+00 -1.00e+00
    //MHD    UU   ENZ  ? P      D 20200317 1346  2.0927 GAU  1.00e-01 -1.00e+00 -1.00e+00 -1.00e+00
    char cStation[7], cInstrument[5], cComponent[5], cPhaseOnset[2],
         cPhase[7], cFirstMotion[2], cError[4];
    std::fill(cStation,     cStation + 7,     '\0');
    std::fill(cInstrument,  cInstrument + 5,  '\0');
    std::fill(cComponent,   cComponent,       '\0');
    std::fill(cPhaseOnset,  cPhaseOnset + 2,  '\0');
    std::fill(cPhase,       cPhase+ 7,        '\0');
    std::fill(cFirstMotion, cFirstMotion + 2, '\0');
    std::fill(cError,       cError + 4,       '\0');
    int year, month, day, hour, minute;
    double seconds, errorMagnitude, codaDuration, amplitude, period, priorWeight;
    sscanf(line.c_str(),
         "%6s %4s %4s %1s %6s %1s %4d%2d%2d %2d%2d %lf %3s %lf %lf %lf %lf %lf",
         cStation, cInstrument, cComponent, cPhaseOnset, cPhase, cFirstMotion,
         &year, &month, &day, &hour, &minute, &seconds,
         cError,
         &errorMagnitude, &codaDuration, &amplitude, &period, &priorWeight);
    std::string station(cStation);
    std::string network(cInstrument);
    std::string component(cComponent);
    std::string phase(cPhase);

    boost::algorithm::trim(station);
    boost::algorithm::trim(network);
    boost::algorithm::trim(component);
    boost::algorithm::trim(phase);
    if (isBlackListed(network, station)){return true;} // Might be blacklisted

    MAssociate::WaveformIdentifier waveid;
    waveid.setNetwork(network);
    waveid.setStation(station);
    waveid.setChannel(component);
    waveid.setLocationCode("01");

    auto polarity = MAssociate::Polarity::UNKNOWN;
    if (cFirstMotion[0] == 'D'){polarity = MAssociate::Polarity::DILATATIONAL;}
    if (cFirstMotion[0] == 'C'){polarity = MAssociate::Polarity::COMPRESSIONAL;}

    SFF::Utilities::Time time;
    time.setYear(year);
    time.setMonth(month);
    time.setDayOfMonth(day);
    time.setHour(hour);
    time.setMinute(minute);
    auto second = static_cast<int> (seconds);
    auto microSecond = static_cast<int> (std::round((seconds - second)*1.e6));
    time.setSecond(second);
    time.setMicroSecond(microSecond);
    //std::cout << line << std::endl;
    //std::cout << waveid << " " << phase << " " << static_cast<int> (polarity) << " " << time << std::endl;
    //getchar();

    double width = 0.1;
    if (heuristicWeight)
    {
        if (phase == "P" && polarity == MAssociate::Polarity::UNKNOWN)
        {
            width = 0.2;
        }
        else
        {
            width = 0.35; 
        }
    }
    else
    {
        if (phase == "P")
        {
            width = errorMagnitude + pModelingError;
        }
        else
        {
            width = errorMagnitude + sModelingError;
        }
    } 
    auto stddev = width/std::sqrt(12.0);  // Uniform distribution

    pick->setWaveformIdentifier(waveid);
    pick->setPhaseName(phase);
    pick->setPolarity(polarity);
    pick->setStandardDeviation(stddev);
    pick->setTime(time.getEpochalTime());
    pick->setPolarityWeight(codaDuration);
    return false;
}

struct StaticCorrection
{
    std::string network;
    std::string station;
    std::string phase;
    double correction;
};

struct StaticCorrections
{
    void read(const std::string &csvFile, const bool header = true)
    {
        corrections.clear();
        std::ifstream infl(csvFile);
        if (!infl){throw std::runtime_error("Could not open: " + csvFile);}
        int lineNumber = 0;
        std::string line;
        std::vector<std::string> split;
        while (std::getline(infl, line)) 
        {
            lineNumber = lineNumber + 1;
            if (header && lineNumber == 1){continue;}
            boost::split(split, line, boost::is_any_of(",\n\t"),
                         boost::token_compress_on);
            StaticCorrection correction;
            correction.network = split[0];
            correction.station = split[1];
            correction.phase = split[2];
            correction.correction = std::stod(split[3]); 
            corrections.push_back(correction);
        }
        infl.close();
    }
    double getCorrection(const std::string &network,
                         const std::string &station,
                         const std::string &phase) const
    {
        for (const auto &correction : corrections)
        {
            if (network == correction.network &&
                station == correction.station &&
                phase   == correction.phase)
            {
                return correction.correction;
            }
         }
         return 0;
    }
    std::vector<StaticCorrection> corrections; 
};

struct PickList
{
    void read(const std::string &fileName,
              const bool heuristicObservation = false,
              const double ptol = 0.2,
              const double stol = 0.4)
    {
        if (heuristicObservation){std::cout << "TODO: RERUN BUT SAVE PICK WIDTH" << std::endl;}
        picks.clear();
        std::ifstream infl(fileName);
        if (!infl){throw std::runtime_error("Could not open: " + fileName);}
        std::string line;
        picks.reserve(2*32768);
        int pickID = 0;
        MAssociate::Pick pick;
        while (std::getline(infl, line))
        {
            pickID = pickID + 1;
            auto blackListed = parsePick(line, &pick, heuristicObservation);
            if (blackListed){continue;} 
            pick.setIdentifier(pickID);
            picks.push_back(pick);
//if (picks.size() == 25065){std::cout << "FIXME: BREAKING EARLY" << std::endl; break;}
        }
        // Attempt to pop duplciates
        auto unknown = MAssociate::Polarity::UNKNOWN;
        std::vector<bool> removePick(picks.size(), false);
        MAssociate::WaveformIdentifier waveid, waveid2;
        std::string network, station, channel, phase;
        for (int ip=0; ip<static_cast<int> (picks.size()); ++ip)
        {
            waveid = picks[ip].getWaveformIdentifier();
            network = waveid.getNetwork();
            station = waveid.getStation();
            channel = waveid.getChannel();
            phase = picks[ip].getPhaseName();
            double tol = ptol;
            if (phase == "S"){tol = stol;}
            auto time = picks[ip].getTime(); 
            if (removePick[ip]){continue;} // Already popped
            if (isCollocated(waveid.getNetwork(), waveid.getStation()))
            {
                for (int jp=ip+1; jp<static_cast<int> (picks.size()); ++jp)
                {
                    auto time2 = picks[jp].getTime();
                    if (std::abs(time - time2) > tol){continue;}
                    waveid2 = picks[jp].getWaveformIdentifier();
                    if (network == waveid2.getNetwork() &&
                        station == waveid2.getStation() &&
                        phase == picks[jp].getPhaseName())
                    {
                        auto channel2 = waveid2.getChannel();
                        // Prefer "better" pick
                        if (picks[ip].getPolarity() == unknown)
                        {
                            if (picks[jp].getPolarity() != unknown)
                            {
                                removePick[ip] = true;
                                removePick[jp] = false;
                                continue;
                            }
                        }
                        if (picks[jp].getPolarity() == unknown)
                        {
                            if (picks[ip].getPolarity() != unknown)
                            {
                                removePick[ip] = false;
                                removePick[jp] = true;
                                continue;
                            } 
                        }
                        // Otherwise prefer high gain to strong motion
                        if (channel[1] == 'H' && channel2[1] == 'N')
                        {
                            removePick[ip] = false;
                            removePick[jp] = true;
                            continue;
                        }
                        else if (channel[1] == 'N' && channel2[1] == 'H')
                        {
                            removePick[ip] = true;
                            removePick[jp] = false;
                            continue;
                        }
                        else
                        {
                            removePick[ip] = true;
                            removePick[jp] = false;
                            std::cout << "Possible duplicate? " << waveid << " " << waveid2 << " " << channel << " " <<  channel2 << std::endl;
                        }
                    }
                }
            }
        }
        int nRemove = std::count(removePick.begin(), removePick.end(), true);
        if (nRemove > 0)
        {
            std::cout << "Removing " + std::to_string(nRemove)
                      << " picks" << std::endl;
            auto temp = picks; 
            picks.clear();
            picks.reserve(temp.size() - nRemove);
            for (int ip=0; ip<static_cast<int> (temp.size()); ++ip)
            {
                if (removePick[ip])
                {
                     //std::cout << "Removing: "
                     //          << temp[ip].getWaveformIdentifier() << std::endl;
                     continue;
                }
                picks.push_back(temp[ip]);
            }
        }
        std::cout << "Number of picks: "
                  << std::to_string(picks.size()) << std::endl;
    }
    void setStatics(const HypoStation &corrections)
    {
        for (int ip=0; ip<static_cast<int> (picks.size()); ++ip)
        {
            auto phase = picks[ip].getPhaseName();
            auto waveid = picks[ip].getWaveformIdentifier();
            auto network = waveid.getNetwork();
            auto station = waveid.getStation();
            auto corr = corrections.getCorrection(network, station, phase);
//std::cout << "P " << corr << " " << station << std::endl;
            if (corr != 0){picks[ip].setStaticCorrection(corr);} 
        }
    }
    void setPStatics(const StaticCorrections &corrections)
    {
        for (int ip=0; ip<static_cast<int> (picks.size()); ++ip)
        {
            auto phase = picks[ip].getPhaseName();
            if (phase != "P"){continue;}
            auto waveid = picks[ip].getWaveformIdentifier();
            auto network = waveid.getNetwork();
            auto station = waveid.getStation();
            auto corr = corrections.getCorrection(network, station, phase);
//std::cout << "P " << corr << " " << station << std::endl;
            if (corr != 0){picks[ip].setStaticCorrection(corr);}
        }
    }
    void setSStatics(const StaticCorrections &corrections)
    {
        for (int ip=0; ip<static_cast<int> (picks.size()); ++ip)
        {
            auto phase = picks[ip].getPhaseName();
            if (phase != "S"){continue;}
            auto waveid = picks[ip].getWaveformIdentifier();
            auto network = waveid.getNetwork();
            auto station = waveid.getStation();
            auto corr = corrections.getCorrection(network, station, phase);
//std::cout << "S " << corr << " " << station << std::endl;
            if (corr != 0){picks[ip].setStaticCorrection(corr);}
        }
    }
    std::vector<MAssociate::Pick> picks;
};


int main(int argc, char *argv[])
{
#ifdef BUILD_MPI
    MPI_Init(&argc, &argv);
#endif
    int nprocs = 1;
    int myid = 0;
#ifdef BUILD_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif
    // Parameters - one day to be read from a config file or something
    std::string travelTimeTableName = "../magna/travelTimeTables.h5";
    //std::string pickFile = "../magna/mainshock.txt";
    //std::string pickFile = "../magna/split_event.txt";
    //std::string pickFile = "../magna/associated_picks.txt";
    int nDays = 44; // 44
//nDays = 40;
    bool heuristicWidth = false;
    std::string pStaticsFile = "/Users/bbaker/trainMagna/p_statics.csv";
    std::string sStaticsFile = "/Users/bbaker/trainMagna/s_statics.csv";
    //std::string hypoStationFile = "../magna/locate/magna.sta"; 
    //std::string hypoStationFile = "../magna/locate/magna.updated.all.sta"; // re-optimized statics
    std::string hypoStationFile = "../magna/locate/magna.updated.all.iter2.sta"; // final iteration
    int dbscanMinClusterSize = 8; // No causality so err on the larger side
    int minNumberOfPArrivals = 2;
    double dbscanEpsilon = 0.35; // (seconds)
    // Load the static corrections
    if (myid == 0){std::cout << "Loading statics..." << std::endl;}
    //StaticCorrections pStatics;
    //StaticCorrections sStatics;
    //pStatics.read(pStaticsFile);
    //sStatics.read(sStaticsFile);
    HypoStation stations;
    stations.read(hypoStationFile);
    // Create the associator parameters
    MAssociate::AssociatorParameters parameters; 
    parameters.setAnalyticCorrelationFunction(
        MAssociate::AnalyticCorrelationFunction::BOXCAR);
    parameters.setOriginTimeObjectiveFunction(
        MAssociate::OriginTimeObjectiveFunction::L1);
    parameters.setDBSCANMinimumClusterSize(dbscanMinClusterSize);
    parameters.setDBSCANEpsilon(dbscanEpsilon);
    // Load the travel time tables
    MAssociate::Mesh::Spherical::Points3D<float> points;
    H5IO h5io;
    std::vector<std::string> tableNames;
    int nTables = 0;
    for (int rank=0; rank<nprocs; ++rank)
    {
        if (rank == myid)
        {
            h5io.openFileForReading(travelTimeTableName);
            h5io.getGeometry(&points); 
            tableNames = h5io.getTravelTimeTableNames(); 
            nTables = static_cast<int> (tableNames.size());
        }
#ifdef BUILD_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
    // Create the associator
    MAssociate::Associator<float> associator;
    parameters.setNumberOfTravelTimeTables(nTables);
    associator.initialize(parameters, points);
    // Add the travel time tables to the associator
    if (myid == 0){std::cout << "Loading travel time tables..." << std::endl;}
    for (int rank=0; rank<nprocs; ++rank)
    {
        if (rank == myid)
        {
            for (const auto &tableName : tableNames)
            {
                std::vector<std::string> nsp; //network/station/phase
                boost::split(nsp, tableName, boost::is_any_of("."));
                std::vector<double> ttimes;
                h5io.getTravelTimeTable(nsp[0], nsp[1], nsp[2], &ttimes);
                associator.setTravelTimeTable(nsp[0], nsp[1], nsp[2],
                                              ttimes.size(), ttimes.data());
            }
        }
#ifdef BUILD_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
    h5io.close();
    uint64_t evidOffset = 10000;
    // Load the picks
    int iday0 = 40;
//nDays = 40;
    for (int jday=iday0; jday<nDays; jday=jday+nprocs) //++iday)
    {
        int iday = jday + myid;
        if (iday >=nDays){break;}
        // Load the picks
        //std::string pickFile = "../magna/associated_picks."
        std::string pickFile = "/data/machineLearning/utahMLPicker/magnaCatalog/associate/associated_picks."
                             + std::to_string(iday+1) + ".txt";
        std::cout << "Loading picks on process " << myid << std::endl;
        PickList picks;
        picks.read(pickFile, heuristicWidth);
        //std::cout << "Attaching static corrections to arrivals..." << std::endl;
        picks.setStatics(stations);
        //picks.setPStatics(pStatics);
        //picks.setSStatics(sStatics);

        associator.clearPicks();
        associator.clearEvents();
        // Attach the picks to the associator
        std::cout << "Setting picks on process " << myid << std::endl; 
        for (const auto &pick : picks.picks)
        {
            associator.addPick(pick); 
        }
        // Associate (this takes awhile)
        std::cout << "Associating day = " << iday
                  << " on process " << myid << std::endl;
        associator.associate();
  
        // Get the events 
        auto events = associator.getEvents();
        // Put them in increasing order of origin time
        std::sort(events.begin(), events.end(),
                  [](const MAssociate::Event &a, const MAssociate::Event &b)
                  {
                      return a.getOriginTime() < b.getOriginTime();
                  });
        std::setprecision(6);
        std::string arrivalFileName = "../magna/prelimArrivals."
                                    + std::to_string(iday+1) + ".csv";
        std::string catalogFileName = "../magna/prelimLocs."
                                    + std::to_string(iday+1) + ".csv";
        std::ofstream arrivalFile(arrivalFileName); //"../magna/prelimArrivals.csv");
        std::ofstream catalogFile(catalogFileName); //"..//magna/prelimLocs.csv");
        arrivalFile << "evid,network,station,channel,location_code,phase,arrival_time,first_motion,first_motion_weight,origin_time,event_latitude,event_longitude,event_depth" << std::endl;
        for (const auto &event : events)
        {
            //if (event.getNumberOfPArrivals() < minNumberOfPArrivals){continue;}
            SFF::Utilities::Time ot(event.getOriginTime());
            auto lon = event.getXPosition();
            auto lat = event.getYPosition();
            auto depth = event.getZPosition();
            //std::cout << ot << "," << lat << "," << lon << "," << depth << std::endl;
            auto arrivals = event.getArrivals();
            for (const auto &arrival : arrivals)
            {
                auto evid = event.getIdentifier() + evidOffset*(iday+1);
                SFF::Utilities::Time arrivalTime(arrival.getTime());
                auto fm = static_cast<int> (arrival.getPolarity());
                auto waveid = arrival.getWaveformIdentifier();
                auto channel = waveid.getChannel();
                if (arrival.getPhaseName() == "S"){channel[2] = 'N';}
                arrivalFile << evid << ","
                            << waveid.getNetwork() << ","
                            << waveid.getStation() << ","
                            << channel << ","
                            << waveid.getLocationCode() << ","
                            << arrival.getPhaseName() << ","
                            << arrivalTime << "," << fm << "," 
                            << arrival.getPolarityWeight() << ","
                            << ot << "," << lat << "," << lon << "," << depth << std::endl; 
            }
            catalogFile << ot << "," << lat << "," << lon << "," << depth << std::endl;
        }
        arrivalFile.close();
        catalogFile.close();
    } // Loop on days
    std::cout << "Process " << myid << " has finished" << std::endl;
#ifdef BUILD_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
