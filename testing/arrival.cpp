#include "massociate/arrival.hpp"
#include "massociate/pick.hpp"
#include "massociate/waveformIdentifier.hpp"
#include <gtest/gtest.h>

namespace
{

TEST(MAssociate, Arrival)
{
    MAssociate::Pick pick;
    double time = 10;
    double std = 0.2;
    double ttime = 4;
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

    EXPECT_NO_THROW(pick.setWaveformIdentifier(sncl));
    EXPECT_EQ(pick.getWaveformIdentifier(), sncl);

    EXPECT_NO_THROW(pick.setPhaseName(phase));
    EXPECT_EQ(pick.getPhaseName(), phase);

    EXPECT_NO_THROW(pick.setTime(time));
    EXPECT_NEAR(pick.getTime(), time, 1.e-14);

    EXPECT_NO_THROW(pick.setPolarity(polarity));
    EXPECT_EQ(pick.getPolarity(), polarity);

    EXPECT_NO_THROW(pick.setStandardDeviation(std));
    EXPECT_NEAR(pick.getStandardDeviation(), std, 1.e-14);

    EXPECT_NEAR(pick.getWeight(), weight, 1.e-14);

    EXPECT_NO_THROW(pick.setIdentifier(id));
    EXPECT_EQ(pick.getIdentifier(), id);

    EXPECT_NO_THROW(pick.setStaticCorrection(staticCorrection));
    EXPECT_NEAR(pick.getStaticCorrection(), staticCorrection, 1.e-14);

    // Create an arrival from a pick
    MAssociate::Arrival copyArrival(pick);
    EXPECT_NO_THROW(copyArrival.setTravelTime(ttime));

    EXPECT_EQ(copyArrival.getWaveformIdentifier(), sncl);
    EXPECT_EQ(copyArrival.getPhaseName(), phase);
    EXPECT_NEAR(copyArrival.getTime(), time, 1.e-14);
    EXPECT_EQ(copyArrival.getPolarity(), polarity);
    EXPECT_NEAR(copyArrival.getWeight(), weight, 1.e-14);
    EXPECT_NEAR(copyArrival.getStandardDeviation(), std, 1.e-14);
    EXPECT_EQ(copyArrival.getIdentifier(), id);
    EXPECT_NEAR(copyArrival.getStaticCorrection(), staticCorrection, 1.e-14);
    EXPECT_NEAR(copyArrival.getTravelTime(), ttime, 1.e-14);
}

}
