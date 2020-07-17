#ifndef MASSOCIATE_ENUMS_HPP
#define MASSOCIATE_ENUMS_HPP
namespace MAssociate
{
/*
 * @brief Defines the polarity of an arrival.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
enum class Polarity
{
    COMPRESSIONAL = 1,  /*!< Upwards polarity. */
    UNKNOWN = 0,        /*!< Unknown polarity. */
    DILATATIONAL  =-1   /*!< Downwards polarity. */
};
/*!
 * @brief Defines the analytic cross-correlation signal to migrate.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
enum class AnalyticCorrelationFunction
{
    GAUSSIAN, /*! Assumes an underlying Gaussian distribution describing the
                  arrival time.  The downside of this model is that it does
                  not have compact support. */
    BOXCAR    /*! Assumes an underlying uniform distribution describing the
                  arrival time. */
};
/*!
 * @brief Defines the travel time field's geometry.
 */
enum class Geometry
{
    SPHERICAL_POINTS_3D,   /*!< 3D points in a spherical geometry,
                                i.e., (latitude, longitude, depth). */
    CARTESIAN_POINTS_3D    /*!< 3D points in a Cartesian geometry,
                                i.e., (x,y,z) */
};

}
#endif
