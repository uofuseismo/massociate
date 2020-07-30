#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <array>
#include <hdf5.h>
#include <massociate/mesh/spherical/points3d.hpp>
#include "h5io.hpp"

namespace
{

std::vector<std::string> getGroupMemberNames(const hid_t gid)
{
    std::array<char, 512> groupName;
    std::fill(groupName.begin(), groupName.end(), '\0');
    auto lenos = H5Iget_name(gid, groupName.data(), groupName.size());
    hsize_t nobj;
    H5Gget_num_objs(gid, &nobj);
    std::array<char, 512> memberName;
    std::vector<std::string> names;
    names.reserve(nobj);
    for (hsize_t i=0; i<nobj; ++i)
    {
        std::fill(memberName.begin(), memberName.end(), '\0');
        lenos = H5Gget_objname_by_idx(gid, i, memberName.data(),
                                      memberName.size());
        std::string name(memberName.data(), lenos);
        names.push_back(name);
    }
    return names;
}

std::string networkStationPhaseToName(const std::string &network,
                                      const std::string &station,
                                      const std::string &phase)
{
    return network + "." + station + "." + phase;
}

}

class H5IO::H5IOImpl
{
public:
    void closeFile()
    {
        mTableNames.clear();
        if (mFileOpenForWriting){H5Fclose(mFile);}
        if (mFileOpenForReading){H5Fclose(mFile);}
        mFileOpenForWriting = false;
        mFileOpenForReading = false;
    }
    ~H5IOImpl()
    {
        closeFile();
    }
    std::vector<std::string> mTableNames;
    hid_t mFile;
    bool mFileOpenForWriting = false;
    bool mFileOpenForReading = false;
};

/// C'tor
H5IO::H5IO() :
    pImpl(std::make_unique<H5IOImpl> ())
{
}

/// Destructor
H5IO::~H5IO() = default;

/// Open file for writing
void H5IO::openFileForWriting(const std::string &fileName)
{
    pImpl->closeFile();
    pImpl->mFile = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC,
                             H5P_DEFAULT, H5P_DEFAULT);
    pImpl->mFileOpenForWriting = true;
}

/// Open file for reading
void H5IO::openFileForReading(const std::string &fileName)
{
    pImpl->closeFile();
    pImpl->mFile = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    pImpl->mFileOpenForReading = true;
    // Get table names
    auto gid = H5Gopen2(pImpl->mFile, "/TravelTimeTables", H5P_DEFAULT);
    pImpl->mTableNames = getGroupMemberNames(gid);
    H5Gclose(gid);
}

/// Sets the geometry
void H5IO::setGeometry(const int n,
                       const double latitudes[],
                       const double longitudes[],
                       const double depths[])
{
    if (!pImpl->mFileOpenForWriting)
    {
        throw std::runtime_error("File not open for writing");
    }
    if (n < 1 || latitudes == nullptr || longitudes == nullptr ||
        depths == nullptr)
    {
        throw std::invalid_argument("n must be positive and no null ptrs");
    }
    // Check if the group exists
    hid_t gid;
    if (!H5Lexists(pImpl->mFile, "/Geometry", H5P_DEFAULT))
    {
        gid = H5Gcreate2(pImpl->mFile, "/Geometry",
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
        gid = H5Gopen2(pImpl->mFile, "/Geometry", H5P_DEFAULT);
    }

    hsize_t dims[1];
    dims[0] = n;
    auto memorySpace = H5Screate_simple(1, dims, NULL);
    auto dataSet = H5Dcreate2(gid, "latitudes",
                              H5T_NATIVE_DOUBLE, memorySpace,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    auto status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, latitudes);
    if (status < 0)
    {
        std::cerr << "Failed writing lats " << std::endl;
    }
    H5Dclose(dataSet);

    dataSet = H5Dcreate2(gid, "longitudes",
                         H5T_NATIVE_DOUBLE, memorySpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, longitudes);
    if (status < 0)
    {
        std::cerr << "Failed writing lons " << std::endl; 
    }
    H5Dclose(dataSet);

    dataSet = H5Dcreate2(gid, "depths",
                         H5T_NATIVE_DOUBLE, memorySpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, depths);
    if (status < 0)
    {
        std::cerr << "Failed writing depths " << std::endl;
    }

    H5Dclose(dataSet);
    H5Sclose(memorySpace);
    H5Gclose(gid);
}

/// Get travel time table names
[[maybe_unused]]
std::vector<std::string> H5IO::getTravelTimeTableNames() const noexcept
{
    return pImpl->mTableNames;
}

/// Determines if the travel time table exists
bool H5IO::haveTravelTimeTable(const std::string &network,
                               const std::string &station,
                               const std::string &phase) const noexcept
{
    auto name = networkStationPhaseToName(network, station, phase);
    for (const auto &tableName : pImpl->mTableNames)
    {
        if (name == tableName){return true;} 
    }
    return false;
}

/// Get the geometry
template<typename U>
void H5IO::getGeometry(
    MAssociate::Mesh::Spherical::Points3D<U> *points) const
{
    if (!pImpl->mFileOpenForReading)
    {
        throw std::runtime_error("File not open for reading");
    }
    auto gid = H5Gopen2(pImpl->mFile, "/Geometry", H5P_DEFAULT);
    if (!H5Lexists(gid, "latitudes", H5P_DEFAULT))
    {
        throw std::runtime_error("latitudes doesnt exist");
    }
    if (!H5Lexists(gid, "longitudes", H5P_DEFAULT))
    {
        throw std::runtime_error("longitudes doesnt exist");
    }
    if (!H5Lexists(gid, "depths", H5P_DEFAULT))
    {
        throw std::runtime_error("depths doesnt exist");
    }
    // Open data, get dimensions, and read
    auto dataSet = H5Dopen(gid, "latitudes", H5P_DEFAULT);
    auto memorySpace = H5Dget_space(dataSet);
    hsize_t dims[1];
    auto ndims = H5Sget_simple_extent_dims(memorySpace, dims, NULL);
    if (ndims != 1){throw std::runtime_error("only arrays");}
    H5Sclose(memorySpace);

    std::vector<double> latitudes(dims[0]);
    H5Dread(dataSet, H5T_NATIVE_DOUBLE,
            H5S_ALL, H5S_ALL, H5P_DEFAULT, latitudes.data());
    H5Dclose(dataSet);

    dataSet = H5Dopen(gid, "longitudes", H5P_DEFAULT);
    std::vector<double> longitudes(dims[0]);
    H5Dread(dataSet, H5T_NATIVE_DOUBLE,
            H5S_ALL, H5S_ALL, H5P_DEFAULT, longitudes.data());
    H5Dclose(dataSet);

    dataSet = H5Dopen(gid, "depths", H5P_DEFAULT);
    std::vector<double> depths(dims[0]);
    H5Dread(dataSet, H5T_NATIVE_DOUBLE,
            H5S_ALL, H5S_ALL, H5P_DEFAULT, depths.data());
    H5Dclose(dataSet);

    H5Gclose(gid);
    // Now create points
    points->setNumberOfPoints(depths.size());
    points->setLatitudes(latitudes.size(), latitudes.data());
    points->setLongitudes(longitudes.size(), longitudes.data());
    points->setDepths(depths.size(), depths.data());
}

/// Reads a travel time table
template<typename U>
void H5IO::getTravelTimeTable(const std::string &network,
                              const std::string &station,
                              const std::string &phase,
                              std::vector<U> *ttimes) const
{
    if (!pImpl->mFileOpenForReading)
    {
        throw std::runtime_error("file not open for reading");
    }
    if (!haveTravelTimeTable(network, station, phase))
    {
        throw std::runtime_error("no corresponding table");
    }
    auto gid = H5Gopen(pImpl->mFile, "/TravelTimeTables", H5P_DEFAULT);
    auto dataSetName = networkStationPhaseToName(network, station, phase);
    auto dataSet = H5Dopen(gid, dataSetName.c_str(), H5P_DEFAULT);
    auto memorySpace = H5Dget_space(dataSet);
    hsize_t dims[1];
    auto ndims = H5Sget_simple_extent_dims(memorySpace, dims, NULL);
    if (ndims != 1){throw std::runtime_error("only arrays");}
    H5Sclose(memorySpace);
    
    std::vector<double> work(dims[0], -1);
    H5Dread(dataSet, H5T_NATIVE_DOUBLE,
            H5S_ALL, H5S_ALL, H5P_DEFAULT, work.data());
    ttimes->resize(work.size());
    std::copy(work.begin(), work.end(), ttimes->data());
    H5Dclose(dataSet);
    H5Gclose(gid);
}

/// Write a travel time table
void H5IO::addTravelTimeTable(const std::string &network,
                              const std::string &station,
                              const std::string &phase,
                              const int n, const double ttimes[])
{
    if (!pImpl->mFileOpenForWriting)
    {
        throw std::runtime_error("File not open for writing");
    }
    if (n < 1 || ttimes == nullptr)
    {
        throw std::invalid_argument("n must be positive and ttimes not null");
    }
    // Check if the group exists
    hid_t gid;
    if (!H5Lexists(pImpl->mFile, "/TravelTimeTables", H5P_DEFAULT))
    {
        gid = H5Gcreate2(pImpl->mFile, "/TravelTimeTables",
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    } 
    else
    {
        gid = H5Gopen2(pImpl->mFile, "/TravelTimeTables", H5P_DEFAULT);
    }
    auto dataSetName = networkStationPhaseToName(network, station, phase);
    // Check if the dataset exists
    if (H5Lexists(pImpl->mFile, dataSetName.c_str(), H5P_DEFAULT))
    {
        std::cerr << "Dataset: " << dataSetName
                  << " already exists" << std::endl;
    }
    hsize_t dims[1];
    dims[0] = n;
    auto memorySpace = H5Screate_simple(1, dims, NULL);
    auto dataSet = H5Dcreate2(gid, dataSetName.c_str(),
                              H5T_NATIVE_DOUBLE, memorySpace,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    std::vector<double> work(ttimes, ttimes+n);
    auto status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, work.data());
    if (status < 0)
    {
        std::cerr << "Failed writing dataset " << dataSetName << std::endl;
    }
    H5Dclose(dataSet);
    H5Sclose(memorySpace);
    H5Gclose(gid);
}

void H5IO::close() noexcept
{
    pImpl->closeFile();
}

///--------------------------------------------------------------------------///
///                        Template Instantiation                            ///
///--------------------------------------------------------------------------///
template void H5IO::getGeometry(MAssociate::Mesh::Spherical::Points3D<double> *points) const;
template void H5IO::getGeometry(MAssociate::Mesh::Spherical::Points3D<float> *points) const;
template void H5IO::getTravelTimeTable(const std::string &network,
                                       const std::string &station,
                                       const std::string &phase,
                                       std::vector<double> *ttimes) const;
template void H5IO::getTravelTimeTable(const std::string &network,
                                       const std::string &station,
                                       const std::string &phase,
                                       std::vector<float> *ttimes) const;
