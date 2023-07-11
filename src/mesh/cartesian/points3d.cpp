#include <string>
#include <stdexcept>
#include <vector>
#include "massociate/mesh/cartesian/points3d.hpp"
#include "massociate/mesh/spherical/points3d.hpp"

using namespace MAssociate::Mesh::Cartesian;

template<class T>
class Points3D<T>::Points3DImpl
{
public:
    std::vector<T> mX;
    std::vector<T> mY;
    std::vector<T> mZ;
    int mScalarFieldSize = 0;
    const MAssociate::Geometry mGeometry
        = MAssociate::Geometry::CARTESIAN_POINTS_3D;
}; 

/// C'tor 
template<class T>
Points3D<T>::Points3D() :
    pImpl(std::make_unique<Points3DImpl> ())
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
Points3D<T>::cloneCartesianPoints3D() const
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
    pImpl->mX.clear();
    pImpl->mY.clear();
    pImpl->mZ.clear();
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

/// Sets/gets X
template<class T>
template<typename U>
void Points3D<T>::setXPositions(int nPoints, const U *x)
{
    auto n = getNumberOfPoints();
    if (nPoints != n)
    {
        throw std::invalid_argument("nPoints = " + std::to_string(nPoints)
                                  + " must = " + std::to_string(n));
    }
    if (x == nullptr)
    {
        throw std::invalid_argument("x cannot be NULL");
    }
    pImpl->mX.resize(n);
    std::copy(x, x+n, pImpl->mX.begin());
}

template<class T>
bool Points3D<T>::haveXPositions() const noexcept
{
    return !pImpl->mX.empty();
}

template<class T>
T Points3D<T>::getXPosition(const size_t index) const
{
    if (!haveXPositions()){throw std::runtime_error("x positions not set");}
    return pImpl->mX.at(index);
}

/// Sets/gets Y
template<class T>
template<typename U>
void Points3D<T>::setYPositions(int nPoints, const U *y)
{
    auto n = getNumberOfPoints();
    if (nPoints != n)
    {
        throw std::invalid_argument("nPoints = " + std::to_string(nPoints)
                                  + " must = " + std::to_string(n));
    }
    if (y == nullptr)
    {
        throw std::invalid_argument("x cannot be NULL");
    }
    pImpl->mY.resize(n);
    std::copy(y, y+n, pImpl->mY.begin());
}

template<class T>
bool Points3D<T>::haveYPositions() const noexcept
{
    return !pImpl->mY.empty();
}

template<class T>
T Points3D<T>::getYPosition(const size_t index) const
{
    if (!haveYPositions()){throw std::runtime_error("y positions not set");}
    return pImpl->mY.at(index);
}

/// Sets/gets Z
template<class T>
template<typename U>
void Points3D<T>::setZPositions(int nPoints, const U *z)
{
    auto n = getNumberOfPoints();
    if (nPoints != n)
    {
        throw std::invalid_argument("nPoints = " + std::to_string(nPoints)
                                  + " must = " + std::to_string(n));
    }
    if (z == nullptr)
    {
        throw std::invalid_argument("x cannot be NULL");
    }
    pImpl->mZ.resize(n);
    std::copy(z, z+n, pImpl->mZ.begin());
}

template<class T>
bool Points3D<T>::haveZPositions() const noexcept
{
    return !pImpl->mZ.empty();
}

template<class T>
T Points3D<T>::getZPosition(const size_t index) const
{
    if (!haveZPositions()){throw std::runtime_error("z positions not set");}
    return pImpl->mZ.at(index);
}

/// Gets the number of points
template<class T>
int Points3D<T>::getNumberOfPoints() const noexcept
{
    return getScalarFieldSize();
}

/// Gets the scalar field size
template<class T>
int Points3D<T>::getScalarFieldSize() const noexcept
{
    return pImpl->mScalarFieldSize;
}

/// Geometry
template<class T>
MAssociate::Geometry Points3D<T>::getGeometry() const noexcept
{
    return pImpl->mGeometry;
}

///--------------------------------------------------------------------------///
///                        Template Instantiation                            ///
///--------------------------------------------------------------------------///
template class MAssociate::Mesh::Cartesian::Points3D<double>;
template class MAssociate::Mesh::Cartesian::Points3D<float>;

template void MAssociate::Mesh::Cartesian::Points3D<double>::setXPositions(const int n, const double x[]);
template void MAssociate::Mesh::Cartesian::Points3D<double>::setXPositions(const int n, const float x[]);
template void MAssociate::Mesh::Cartesian::Points3D<double>::setYPositions(const int n, const double x[]);
template void MAssociate::Mesh::Cartesian::Points3D<double>::setYPositions(const int n, const float x[]);
template void MAssociate::Mesh::Cartesian::Points3D<double>::setZPositions(const int n, const double x[]);
template void MAssociate::Mesh::Cartesian::Points3D<double>::setZPositions(const int n, const float x[]);

template void MAssociate::Mesh::Cartesian::Points3D<float>::setXPositions(const int n, const double x[]);
template void MAssociate::Mesh::Cartesian::Points3D<float>::setXPositions(const int n, const float x[]);
template void MAssociate::Mesh::Cartesian::Points3D<float>::setYPositions(const int n, const double x[]);
template void MAssociate::Mesh::Cartesian::Points3D<float>::setYPositions(const int n, const float x[]);
template void MAssociate::Mesh::Cartesian::Points3D<float>::setZPositions(const int n, const double x[]);
template void MAssociate::Mesh::Cartesian::Points3D<float>::setZPositions(const int n, const float x[]);
