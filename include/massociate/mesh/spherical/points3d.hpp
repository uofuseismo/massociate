#ifndef MASSOCIATE_MESH_SPHERICAL_POINTS3D_HPP
#define MASSOCIATE_MESH_SPHERICAL_POINTS3D_HPP
#include <memory>
#include "massociate/mesh/mesh.hpp"
namespace MAssociate::Mesh::Spherical
{
template<class T = float>
class Points3D : public MAssociate::Mesh::IMesh<T>
{
public:
    /*!
     * @brief Constructor.
     */
    Points3D();
    /*!
     * @brief Copy constructor.
     * @param[in] points   The points class from which to initialize this class.
     */
    Points3D(const Points3D &points);

    /*!
     * @brief Destructor.
     */
    virtual ~Points3D();
    /*!
     * @brief Clears the memory/resets the class.
     */
    void clear() noexcept override;
    /*!
     * @brief Copy assignment.
     * @param points   The points class to copy to this.
     * @return A deep copy of the input points.
     */
    Points3D& operator=(const Points3D &points);
    /*!
     * @brief Clone self.
     * @return A deep copy of this.
     */
    [[nodiscard]] std::unique_ptr<Points3D<T>> cloneSphericalPoints3D() const override;

    /*!
     * @brief Sets the number of points in the field.
     * @param nPoints   The number of points.
     * @throws std::invalid_argument if nPoints is not positive.
     */
    void setNumberOfPoints(int nPoints);
    /*!
     * @return The number of points in the field.
     */
    [[nodiscard]] int getNumberOfPoints() const noexcept;

    /*!
     * @brief Sets the latitudes for each point.
     * @param[in] nPoints    The number of points.
     * @param[in] latitudes  The latitudes in degrees for each point.  This is
     *                       an array of dimension [nPoints].
     * @throws std::invalid_argument if any latitude is not in the range
     *         [-90,90], nPoints does not match \c getNumberOfPoints(),
     *         or latitudes is NULL.
     */
    template<typename U>
    void setLatitudes(int nPoints,  const U latitudes[]);
    /*!
     * @return True indicates that the latitudes are set.
     */
    [[nodiscard]] bool haveLatitudes() const noexcept;

    /*!
     * @brief Sets the longitudes for each point.
     * @param[in] nPoints     The number of points.
     * @param[in] longitudes  The longitudes in degrees for each point.  This is
     *                        an array of dimension [nPoints].
     * @throws std::invalid_argument if any longitude is not in the range
     *         [-540,540), nPoints does not match \c getNumberOfPoints(),
     *         or longitudes is NULL.
     */
    template<typename U>
    void setLongitudes(int nPoints, const U longitudes[]);
    /*!
     * @return True indicates that the longitudes are set.
     */
    [[nodiscard]] bool haveLongitudes() const noexcept;

    /*!
     * @brief Sets the depths for each point.
     * @param[in] nPoints     The number of points.
     * @param[in] depths      The depths in kilometers for each point.  This is
     *                        an array of dimension [nPoints].
     * @throws std::invalid_argument if nPoints does not match
     *         \c getNumberOfPoints() or depths is NULL.
     */
    template<typename U>
    void setDepths(int nPoints, const U depths[]);
    /*!
     * @return True indicates that the depths are set.
     */
    [[nodiscard]] bool haveDepths() const noexcept;

    /*!
     * @param index  The index of the underlying latitudes vector to access.
     * @return The latitude corresponding to the given index.
     * @throws std::out_of_range if the index exceeds the number of
     *         points.
     * @throws std::runtime_error if \c haveLatitudes() is false.
     * @sa \c setLatitudes()
     */
    [[nodiscard]] T getLatitude(size_t index) const;
    /*!
     * @param index  The index of the underlying longitudes vector to access.
     * @return The longitude corresponding to the given index.
     * @throws std::out_of_range if the index exceeds the number of
     *         points.
     * @throws std::runtime_error if \c haveLongitudes() is false.
     * @sa \c setLongitudes()
     */
    [[nodiscard]] T getLongitude(size_t index) const;
    /*!
     * @param index  The index of the underlying depths vector to access.
     * @return The depth corresponding to the given index.
     * @throws std::out_of_range if the index exceeds the number of
     *         points.
     * @throws std::runtime_error if \c haveDepths() is false.
     * @sa \c setDepths()
     */
    [[nodiscard]] T getDepth(size_t index) const;
    [[nodiscard]] int getScalarFieldSize() const noexcept override;

    /*!
     * @return The geometry - i.e., SPHERICAL_POINTS_3D.
     */
    [[nodiscard]] MAssociate::Geometry getGeometry() const noexcept override;

private:
    class Points3DImpl;
    std::unique_ptr<Points3DImpl> pImpl;
};

}
#endif
