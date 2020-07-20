#include "massociate/arrival.hpp"
#include "massociate/waveformIdentifier.hpp"
#include <gtest/gtest.h>

namespace
{

TEST(MAssociate, Arrival)
{
    MAssociate::Arrival arrival;
    double time = 10;
    double std = 0.2;
    double weight = 1/std;
    double staticCorrection = 0.2;
    uint64_t id = 66823;
    std::string phase = "P";
    auto polarity = MAssociate::Polarity::DILATATIONAL;
    MAssociate::WaveformIdentifier sncl;
    std::string network = "UU";
    std::string station = "CAPU";
    std::string channel = "EHZ";
    std::string location = "01";
    sncl.setNetwork(network);
    sncl.setStation(station);
    sncl.setChannel(channel);
    sncl.setLocationCode(location);

    EXPECT_NO_THROW(arrival.setWaveformIdentifier(sncl));
    EXPECT_EQ(arrival.getWaveformIdentifier(), sncl);

    EXPECT_NO_THROW(arrival.setPhaseName(phase));
    EXPECT_EQ(arrival.getPhaseName(), phase);

    EXPECT_NO_THROW(arrival.setTime(time));
    EXPECT_NEAR(arrival.getTime(), time, 1.e-14);

    EXPECT_NO_THROW(arrival.setPolarity(polarity));
    EXPECT_EQ(arrival.getPolarity(), polarity);

    EXPECT_NO_THROW(arrival.setStandardDeviation(std));
    EXPECT_NEAR(arrival.getStandardDeviation(), std, 1.e-14);

    EXPECT_NEAR(arrival.getWeight(), weight, 1.e-14);

    EXPECT_NO_THROW(arrival.setIdentifier(id));
    EXPECT_EQ(arrival.getIdentifier(), id);

    EXPECT_NO_THROW(arrival.setStaticCorrection(staticCorrection));
    EXPECT_NEAR(arrival.getStaticCorrection(), staticCorrection, 1.e-14);

    // Copy
    MAssociate::Arrival copyArrival(arrival);
    EXPECT_EQ(copyArrival.getWaveformIdentifier(), sncl);
    EXPECT_EQ(copyArrival.getPhaseName(), phase);
    EXPECT_NEAR(copyArrival.getTime(), time, 1.e-14);
    EXPECT_EQ(copyArrival.getPolarity(), polarity);
    EXPECT_NEAR(copyArrival.getWeight(), weight, 1.e-14);
    EXPECT_NEAR(copyArrival.getStandardDeviation(), std, 1.e-14);
    EXPECT_EQ(copyArrival.getIdentifier(), id);
    EXPECT_NEAR(copyArrival.getStaticCorrection(), staticCorrection, 1.e-14);
}

}
