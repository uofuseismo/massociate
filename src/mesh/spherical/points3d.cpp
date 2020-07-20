#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include "massociate/mesh/spherical/points3d.hpp"
#include "massociate/mesh/cartesian/points3d.hpp"

using namespace MAssociate::Mesh;
using namespace MAssociate::Mesh::Spherical;

template<class T>
class Points3D<T>::Points3DImpl
{
public:
    //std::map<std::string, std::vector<T>> mScalarNodalFields;
    std::vector<T> mLatitudes;
    std::vector<T> mLongitudes;
    std::vector<T> mDepths;
    int mScalarFieldSize = 0;
    const MAssociate::Geometry mGeometry
       = MAssociate::Geometry::SPHERICAL_POINTS_3D;
};

/// Constructor
template<class T>
Points3D<T>::Points3D() :
    pImpl(std::make_unique<Points3DImpl>())
{
}

/// Copy c'tor
template<class T>
Points3D<T>::Points3D(const Points3D<T> &points)
{
    *this = points;
}

/// Copy assignment
template<class T>
Points3D<T>& Points3D<T>::operator=(const Points3D &points)
{
    if (&points == this){return *this;}
    pImpl = std::make_unique<Points3DImpl> (*points.pImpl);
    return *this;
}

/// Deep copy
template<class T>
std::unique_ptr<Points3D<T>>
Points3D<T>::cloneSphericalPoints3D() const
{
    auto result = std::make_unique<Points3D<T>> (*this);
    return result;
}

/// Destructor
template<class T>
Points3D<T>::~Points3D() = default;

/// Clears the class
template<class T>
void Points3D<T>::clear() noexcept
{
    pImpl->mLatitudes.clear();
    pImpl->mLongitudes.clear();
    pImpl->mDepths.clear();
    pImpl->mScalarFieldSize = 0;
}

/// Sets the number of points
template<class T>
void Points3D<T>::setNumberOfPoints(const int nPoints)
{
    if (nPoints < 1)
    {
        throw std::invalid_argument("Number of points must be positive");
    }
    clear();
    pImpl->mScalarFieldSize = nPoints;
}

/// Gets the number of points
template<class T>
int Points3D<T>::getNumberOfPoints() const noexcept
{
    return getScalarFieldSize();
}

/// Sets the latitudes
template<class T>
template<typename U>
void Points3D<T>::setLatitudes(int nPoints, const U *latitudes)
{
    auto n = getNumberOfPoints();
    if (nPoints != n)
    {
        throw std::invalid_argument("nPoints = " + std::to_string(nPoints)
                                        + " must = " + std::to_string(n));
    }
    if (latitudes == nullptr)
    {
        throw std::invalid_argument("latitudes cannot be NULL");
    }
    const auto[min, max] = std::minmax_element(latitudes, latitudes + n);
    if (*min < -90 || *max > 90)
    {
        throw std::invalid_argument("All latitudes must be in range [90,90]");
    }
    pImpl->mLatitudes.resize(n);
    std::copy(latitudes, latitudes + n, pImpl->mLatitudes.begin());
}

template<class T>
bool Points3D<T>::haveLatitudes() const noexcept
{
    return !pImpl->mLatitudes.empty();
}

/// Sets the longitudes
template<class T>
template<typename U>
void Points3D<T>::setLongitudes(int nPoints, const U *longitudes)
{
    auto n = getNumberOfPoints();
    if (nPoints != n)
    {
        throw std::invalid_argument("nPoints = " + std::to_string(nPoints)
                                        + " must = " + std::to_string(n));
    }
    if (longitudes == nullptr)
    {
        throw std::invalid_argument("longitudes cannot be NULL");
    }
    const auto[min, max] = std::minmax_element(longitudes, longitudes + n);
    if (*min < -549 || *max >= 540)
    {
        throw std::invalid_argument(
            "All longitudes must be in range [-540,540)");
    }
    pImpl->mLongitudes.resize(n);
    std::copy(longitudes, longitudes + n, pImpl->mLongitudes.begin());
}

template<class T>
bool Points3D<T>::haveLongitudes() const noexcept
{
    return !pImpl->mLongitudes.empty();
}

/// Sets the depths
template<class T>
template<typename U>
void Points3D<T>::setDepths(int nPoints, const U *depths)
{
    auto n = getNumberOfPoints();
    if (nPoints != n)
    {
        throw std::invalid_argument("nPoints = " + std::to_string(nPoints)
                                        + " must = " + std::to_string(n));
    }
    if (depths == nullptr)
    {
        throw std::invalid_argument("depths cannot be NULL");
    }
    pImpl->mDepths.resize(n);
    std::copy(depths, depths+n, pImpl->mDepths.begin());
}

template<class T>
bool Points3D<T>::haveDepths() const noexcept
{
    return !pImpl->mDepths.empty();
}

/// Gets the scalar field size
template<class T>
int Points3D<T>::getScalarFieldSize() const noexcept
{
    return pImpl->mScalarFieldSize;
}

/// Get longitude corresponding to index
template<class T>
T Points3D<T>::getLongitude(const size_t index) const
{
    if (!haveLongitudes()){throw std::runtime_error("Longitudes not set");}
    return pImpl->mLongitudes.at(index);
}

/// Get latitude corresponding to index
template<class T>
T Points3D<T>::getLatitude(const size_t index) const
{
    if (!haveLatitudes()){throw std::runtime_error("Latitudes not set");}
    return pImpl->mLatitudes.at(index);
}

/// Get depth corresponding to index
template<class T>
T Points3D<T>::getDepth(const size_t index) const
{
    if (!haveDepths()){throw std::runtime_error("Depths not set");}
    return pImpl->mDepths.at(index);
}

/// Get geometry
template<class T>
MAssociate::Geometry Points3D<T>::getGeometry() const noexcept
{
    return pImpl->mGeometry;
}

///--------------------------------------------------------------------------///
///                        Template Instantiation                            ///
///--------------------------------------------------------------------------///
template class MAssociate::Mesh::Spherical::Points3D<double>;
template class MAssociate::Mesh::Spherical::Points3D<float>;
 
template void MAssociate::Mesh::Spherical::Points3D<double>::setLatitudes(const int n, const double x[]);
template void MAssociate::Mesh::Spherical::Points3D<double>::setLatitudes(const int n, const float x[]);
template void MAssociate::Mesh::Spherical::Points3D<double>::setLongitudes(const int n, const double x[]);
template void MAssociate::Mesh::Spherical::Points3D<double>::setLongitudes(const int n, const float x[]);
template void MAssociate::Mesh::Spherical::Points3D<double>::setDepths(const int n, const double x[]);
template void MAssociate::Mesh::Spherical::Points3D<double>::setDepths(const int n, const float x[]);

template void MAssociate::Mesh::Spherical::Points3D<float>::setLatitudes(const int n, const double x[]);
template void MAssociate::Mesh::Spherical::Points3D<float>::setLatitudes(const int n, const float x[]);
template void MAssociate::Mesh::Spherical::Points3D<float>::setLongitudes(const int n, const double x[]);
template void MAssociate::Mesh::Spherical::Points3D<float>::setLongitudes(const int n, const float x[]);
template void MAssociate::Mesh::Spherical::Points3D<float>::setDepths(const int n, const double x[]);
template void MAssociate::Mesh::Spherical::Points3D<float>::setDepths(const int n, const float x[]);

