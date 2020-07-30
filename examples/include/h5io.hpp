#ifndef H5IO_HPP
#define H5IO_HPP
#include <memory>
#include <massociate/mesh/spherical/points3d.hpp>
/*!
 * @brief A class for reading/writing the 3D points fields defining the 
 *        travel time fields.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class H5IO
{
public:
    /*!
     * @brief Constructor.
     */
    H5IO();
    /*!
     * @brief Destructor.
     */
    ~H5IO();

    /*! @name Archive Creation
     * @{
     */
    /*!
     * @brief Open file for writing.
     * @param[in] fileName  The file name to open for writing.
     * @note This will completely overwrite an existin file.
     */
    void openFileForWriting(const std::string &fileName); 
    /*!
     * @brief Sets the 3D spherical points geometry. 
     * @param[in] n           The number of points in the travel time field.
     * @param[in] latitudes   The latitudes of the candidate sources in 
     *                        degrees.  This is an array whose dimension is [n].
     * @param[in] longitudes  The longitudes of the candidate sources in
     *                        degrees.  This is an array whose dimension is [n].
     * @param[in] depths      The depths of the candidate sources in kilometers.
     *                        This is an array whose dimension is [n]. 
     * @throw std::invalid_argument if n is not positive or any of the arrays
     *        are NULL.
     * @throws std::runtime_error if the file is not open for writing.
     */
    void setGeometry(const int n,
                     const double latitudes[],
                     const double longitudes[],
                     const double depths[]);
    /*!
     * @brief Adds the travel time from the receiver to all points.
     * @param[in] network  The network name - likely UU.
     * @param[in] station  The station name.
     * @param[in] phase    The phase of the travel time field - likely P or S.
     * @param[in] n        The number of points in the travel time field.
     * @param[in] travelTimes  The travel times in seconds from the receiver
     *                         to all sources.  This is an array whose dimension
     *                         is [n].
     * @throws std::invalid_argument if n is not positive or travelTimes
     *         is NULL.
     * @throws std::runtime_error if the file is not open for writing.
     */
    void addTravelTimeTable(const std::string &network,
                            const std::string &station,
                            const std::string &phase,
                            const int n, const double travelTimes[]);
    /*! @} */

    /*! @name Archive Reading
     * @{
     */
    /*!
     * @brief Opens a travel time table archive file for reading.
     * @param[in] fileName   The HDF5 file to open for reading.
     */
    void openFileForReading(const std::string &fileName); 
    /*!
     * @param[in] network  The network name.
     * @param[in] station  The station name.
     * @param[in] phase    The phase name.
     * @result True indicates that the corresponding travel time table exists.
     */
    bool haveTravelTimeTable(const std::string &network,
                             const std::string &station,
                             const std::string &phase) const noexcept;
    /*!
     * @result The name of all the travel time tables in the archive.
     */
    std::vector<std::string> getTravelTimeTableNames() const noexcept;
    /*!
     * @brief Gets the geometry.
     * @param[out] geometry  The spherical 3D points geometry defining the
     *                       each travel time table's geometry.
     * @throws std::runtime_error if the file isn't open for reading.
     */
    template<typename U>
    void getGeometry(MAssociate::Mesh::Spherical::Points3D<U> *geometry) const;

    /*!
     * @brief Gets the travel time table corresponding to the
     *        network/station/phase.
     * @param[in] network  The network name.
     * @param[in] station  The station name.
     * @param[in] phase    The phase name.
     * @param[out] ttimes  The travel times in seconds from the station to
     *                     each candidate source point in the geometry.
     * @sa \c haveTravelTimeTable()
     */
    template<typename U>
    void getTravelTimeTable(const std::string &network,
                            const std::string &station,
                            const std::string &phase,
                            std::vector<U> *ttimes) const;
    /*! @} */ 

    /*!
     * @brief Closes a file.
     */
    void close() noexcept;
private:
    H5IO(const H5IO &h5io) = delete;
    H5IO(H5IO &&h5io) noexcept = delete;
    H5IO& operator=(const H5IO &h5io) = delete;
    H5IO& operator=(H5IO &&h5io) noexcept = delete;
    class H5IOImpl;
    std::unique_ptr<H5IOImpl> pImpl; 
};
#endif
