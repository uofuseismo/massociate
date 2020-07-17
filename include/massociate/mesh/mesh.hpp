#ifndef MASSOCIATE_MESH_IMESH_HPP
#define MASSOCIATE_MESH_IMESH_HPP
#include <memory>
#include "massociate/enums.hpp"
namespace MAssociate::Mesh
{
namespace Spherical
{
  template<class T> class Points3D;
}
namespace Cartesian
{
  template<class T> class Points3D;
}
}
namespace MAssociate::Mesh
{
/*!
 * @brief An abstract base class that defines a mesh geometry.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
template<class T = float>
class IMesh
{
public:
    /*!
     * @brief Default destructor.
     */
    virtual ~IMesh() = default;
    /*!
     * @brief Resets the class.
     */
    virtual void clear() noexcept = 0;
    /*!
     * @brief Ensures the abstract base class can copy itself.
     */
    virtual std::unique_ptr<Spherical::Points3D<T>> cloneSphericalPoints3D() const{return nullptr;}
    virtual std::unique_ptr<Cartesian::Points3D<T>> cloneCartesianPoints3D() const{return nullptr;}
    /*!
     * @result The number of points in the scalar field.
     */
    [[nodiscard]] virtual int getScalarFieldSize() const = 0;
    /*!
     * @param[in] fieldName  The name of the scalar field to return.
     * @result A pointer to the field of scalar values.  This is an
     *         array whose dimensions is [\c getScalarFieldSize() ].
     * @throws std::invalid_argument if the scalar field field cannot be found.
     */
    //virtual const T* getScalarFieldPointer(const std::string &fieldName) const = 0;
    /*!
     * @result The geometry type.
     */
    virtual MAssociate::Geometry getGeometry() const noexcept = 0;
};
}
#endif
