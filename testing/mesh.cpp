#include <vector>
#include "massociate/mesh/spherical/points3d.hpp"
#include <gtest/gtest.h>

namespace
{
using namespace MAssociate::Mesh;

TEST(MAssociate, MeshSpherical3D)
{
    Spherical::Points3D<double> points;
    int npts = 5;
    std::vector<double> lats({1, 2, 3, 4, 5});
    std::vector<double> lons({-1, -2, -3, -4, -5});
    std::vector<double> depths({9, 10, 11, 12, 13}); 
    EXPECT_EQ(points.getGeometry(), MAssociate::Geometry::SPHERICAL_POINTS_3D);
    EXPECT_NO_THROW(points.setNumberOfPoints(npts));
    EXPECT_NO_THROW(points.setLatitudes(lats.size(), lats.data()));
    EXPECT_NO_THROW(points.setLongitudes(lons.size(), lons.data()));
    EXPECT_NO_THROW(points.setDepths(depths.size(), depths.data()));
    EXPECT_EQ(points.getNumberOfPoints(), npts);
    EXPECT_EQ(points.getScalarFieldSize(), npts);
    for (int i=0; i<npts; ++i)
    {
        EXPECT_NEAR(points.getLatitude(i), lats[i], 1.e-14);
        EXPECT_NEAR(points.getLongitude(i), lons[i], 1.e-14);
        EXPECT_NEAR(points.getDepth(i), depths[i], 1.e-14);
    }

    auto pointsCopy = points.cloneSphericalPoints3D();
    for (int i=0; i<npts; ++i)
    {
        EXPECT_NEAR(pointsCopy->getLatitude(i), lats[i], 1.e-14);
        EXPECT_NEAR(pointsCopy->getLongitude(i), lons[i], 1.e-14);
        EXPECT_NEAR(pointsCopy->getDepth(i), depths[i], 1.e-14);
    } 
}

TEST(MAssociate, MeshCartesian3D)
{

}

}
