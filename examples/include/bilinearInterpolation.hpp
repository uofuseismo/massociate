#ifndef MASSOCIATE_EXAMPLES_BILINEARINTERPOLATION_HPP
#define MASSOCIATE_EXAMPLES_BILINEARINTERPOLATION_HPP
#include <memory>
template<class T=double>
class BilinearInterpolation
{
public:
    /*!
     * @brief Constructor.
     */
    BilinearInterpolation();

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~BilinearInterpolation();
    /*!
     * @brief Releases memory and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Step 1: Initialization
     * @{
     */
    /*!
     * @brief Initializes the bilinear interpolation engine for use with a
     *        non-uniform grid.
     * @param[in] nx   The number of grid points in x.  This must be at least 2.
     * @param[in] x    The x points of the grid.  This is an array whose
     *                 dimension is [nx] and must be in strictly increasing
     *                 order.
     * @param[in] ny   The number of grid points in y.  This must be at least 2.
     * @param[in] y    The y points of the grid.  This is an array whose
     *                 dimension is [ny] and must be in strictly increasing
     *                 order.
     * @param[in] nxy  The number of function values.  This must equal nx*ny.
     * @param[in] f    The function values at each grid point.  This is a row
     *                 major matrix of dimension [ny x nx].   For example,
     *                 x could be offfset and y could be depth.
     * @throws std::invalid_argument if nx or ny is too small, any of the
     *         pointers are NULL, or x and y are not sorted in increasing order.
     */
    template<typename U>
    void initialize(int nx, const U x[], 
                    int ny, const U y[],
                    int nxy, const U f[]);
    /*!
     * @brief Initializes the bilinear interpolation for use with a uniform
     *        grid.
     */
    template<typename U>
    void initialize(const std::pair<U, U> &xLimits,
                    const std::pair<U, U> &yLimits,
                    int nx, int ny, const U f[]);
    /*!
     * @result True indicates that the class is initialized.
     */
    bool isInitialized() const noexcept;
    /*! @} */

    /*! @name Step 2: Interpolation
     * @{
     */
    /*!
     * @brief Bilinearly interpolates f at points \f$ (x_q, y_q) \f$.
     * @param[in] nq   The number of points at which to interpolate.
     * @param[in] xq   The x abscissas at which to interpolate.  This is an
     *                 array whose dimension is [nq].
     * @param[in] yq   The y abscissas at which to interpolate.  This is an
     *                 array whose dimension is [nq].
     * @throws std::runtime_error if \c isInitialized() is false.
     */
    void interpolate(int nq, const T xq[], const T yq[]);
    /*!
     * @result True indicates that the interpolated function was computed.
     */
    bool haveInterpolatedFunction() const noexcept;
    /*! @} */

    /*! @name Step 3: Get Result
     * @{
     */
    std::vector<T> getInterpolatedFunction() const;
    const T *getInterpolatedFunctionPointer() const;
  
    int getInterpolatedFunctionSize() const noexcept;
    /*! @} */
private:
    class BilinearInterpolationImpl;
    std::unique_ptr<BilinearInterpolationImpl> pImpl;
};
#endif
