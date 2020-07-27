#include <vector>
#include "massociate/event.hpp"
#include "massociate/waveformIdentifier.hpp"
#include "massociate/arrival.hpp"
#include <gtest/gtest.h>
namespace
{
TEST(MAssociator, Event)
{
    uint64_t evid = 5923;
    double latitude = 40.7608;
    double longitude =-111.8910;
    double depth = 5000;
    double originTime = 50;
    double x = 2;
    double y = 3;
    double z = 4;
    MAssociate::Arrival arrival;
    MAssociate::WaveformIdentifier waveid;
    MAssociate::Event event;
 
    EXPECT_NO_THROW(event.setIdentifier(evid));
    EXPECT_EQ(event.getIdentifier(), evid);

    EXPECT_NO_THROW(event.setLatitude(latitude));
    EXPECT_NEAR(event.getLatitude(), latitude, 1.e-10);

    EXPECT_NO_THROW(event.setLongitude(longitude));
    EXPECT_NEAR(event.getLongitude(), longitude + 360, 1.e-10);
    event.setLongitude(event.getLongitude());
    EXPECT_NEAR(event.getLongitude(), longitude + 360, 1.e-10);

    EXPECT_NO_THROW(event.setDepth(depth));
    EXPECT_NEAR(event.getDepth(), depth, 1.e-10);

    EXPECT_NO_THROW(event.setOriginTime(originTime));
    EXPECT_NEAR(event.getOriginTime(), originTime, 1.e-10);

    event.setXPosition(x);
    EXPECT_NEAR(event.getXPosition(), x, 1.e-14);
    event.setYPosition(y);
    EXPECT_NEAR(event.getYPosition(), y, 1.e-14);
    event.setZPosition(z);
    EXPECT_NEAR(event.getZPosition(), z, 1.e-14);

    // Create some arrivals - check this is copy c'tor
    std::vector<MAssociate::Arrival> arrivals;
    waveid.setNetwork("UU");
    waveid.setStation("LKWY");
    waveid.setChannel("ENZ");
    waveid.setLocationCode("00");
    arrival.setIdentifier(0);
    arrival.setTime(4);
    arrival.setWaveformIdentifier(waveid);
    arrival.setPhaseName("P");
    arrivals.push_back(arrival);
 
    waveid.setStation("VEC");
    arrival.setTime(8);
    arrival.setWaveformIdentifier(waveid);
    arrival.setPhaseName("S");
    arrival.setIdentifier(1);
    arrivals.push_back(arrival);

    waveid.setStation("COY");
    arrival.setTime(6);
    arrival.setWaveformIdentifier(waveid);
    arrival.setPhaseName("P");
    arrival.setIdentifier(2);
    arrivals.push_back(arrival);
    for (const auto &arrival : arrivals)
    {
        EXPECT_NO_THROW(event.addArrival(arrival));
    }
    EXPECT_EQ(event.getNumberOfArrivals(), 3);
    EXPECT_EQ(event.getNumberOfPArrivals(), 2);

    // Copy c'tor
    MAssociate::Event eCopy(event);
    EXPECT_EQ(eCopy.getIdentifier(), evid);
    EXPECT_NEAR(eCopy.getLatitude(), latitude, 1.e-10);
    EXPECT_NEAR(eCopy.getLongitude(), longitude + 360, 1.e-10);
    EXPECT_NEAR(eCopy.getDepth(), depth, 1.e-10);
    EXPECT_NEAR(eCopy.getOriginTime(), originTime, 1.e-10); 
    EXPECT_NEAR(eCopy.getXPosition(), x, 1.e-14);
    EXPECT_NEAR(eCopy.getYPosition(), y, 1.e-14);
    EXPECT_NEAR(eCopy.getZPosition(), z, 1.e-14);
    auto arrivalsBack = eCopy.getArrivals(); 
    for (int ia=0; ia<static_cast<int> (arrivals.size()); ++ia)
    {
        EXPECT_EQ(arrivalsBack[ia].getWaveformIdentifier(),
                  arrivals[ia].getWaveformIdentifier());
        EXPECT_EQ(arrivalsBack[ia].getPhaseName(), arrivals[ia].getPhaseName());
        EXPECT_NEAR(arrivalsBack[ia].getTime(), arrivals[ia].getTime(), 1.e-14);
    }
    event.clearArrivals();
    EXPECT_EQ(event.getNumberOfArrivals(), 0);
}

}
