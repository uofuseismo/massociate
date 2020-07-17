#ifndef MASSOCIATE_MESH_CARTESIAN_POINTS3D_HPP
#define MASSOCIATE_MESH_CARTESIAN_POINTS3D_HPP
#include <memory>
#include "massociate/mesh/mesh.hpp"
namespace MAssociate::Mesh::Cartesian
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
    [[nodiscard]] std::unique_ptr<Points3D<T>> cloneCartesianPoints3D() const override;

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
     * @brief Sets the x position for each point.
     * @param[in] nPoints   The number of points.
     * @param[in] x         The x position of each point.  This is an array of
     *                      dimension [nPoints].
     * @throws std::invalid_argument if nPoints does not match
     *         \c getNumberOfPoints() or x is NULL.
     */
    template<typename U>
    void setXPositions(int nPoints, const U x[]);
    /*!
     * @result True indicates that the x positions are set.
     */
    [[nodiscard]] bool haveXPositions() const noexcept;

    /*!
     * @brief Sets the y position for each point.
     * @param[in] nPoints   The number of points.
     * @param[in] y         The y position of each point.  This is an array of
     *                      dimension [nPoints].
     * @throws std::invalid_argument if nPoints does not match
     *         \c getNumberOfPoints() or y is NULL.
     */
    template<typename U>
    void setYPositions(int nPoints, const U y[]);
    /*!
     * @result True indicates that the y positions are set.
     */
    [[nodiscard]] bool haveYPositions() const noexcept;

    /*!
     * @brief Sets the z position for each point.
     * @param[in] nPoints   The number of points.
     * @param[in] z         The z position of each point.  This is an array of
     *                      dimension [nPoints].
     * @throws std::invalid_argument if nPoints does not match
     *         \c getNumberOfPoints() or z is NULL.
     */
    template<typename U>
    void setZPositions(int nPoints, const U z[]);
    /*!
     * @result True indicates that the z positions are set.
     */
    [[nodiscard]] bool haveZPositions() const noexcept;

    /*!
     * @param[in] index   The index of the underlying x vector to access.
     * @result The x position at the given index.
     * @throws std::out_of_range if the index exceeds the number of points.
     * @throws std::runtime_error if \c haveXPositions() is false.
     */
    [[nodiscard]] T getXPosition(size_t index) const;
    /*!
     * @param[in] index   The index of the underlying y vector to access.
     * @result The y position at the given index.
     * @throws std::out_of_range if the index exceeds the number of points.
     * @throws std::runtime_error if \c haveXPositions() is false.
     */
    [[nodiscard]] T getYPosition(size_t index) const;
    /*!
     * @param[in] index   The index of the underlying z vector to access.
     * @result The z position at the given index.
     * @throws std::out_of_range if the index exceeds the number of points.
     * @throws std::runtime_error if \c haveZPositions() is false.
     */
    [[nodiscard]] T getZPosition(size_t index) const;

    [[nodiscard]] int getScalarFieldSize() const noexcept override;
    /*!
     * @return The geometry - i.e., CARTESIAN_POINTS_3D.
     */
    [[nodiscard]] MAssociate::Geometry getGeometry() const noexcept override;
private:
    class Points3DImpl;
    std::unique_ptr<Points3DImpl> pImpl;
};
}
#endif
