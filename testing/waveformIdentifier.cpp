#include "massociate/waveformIdentifier.hpp"
#include <gtest/gtest.h>

namespace
{
TEST(MAssociate, WaveformIdentifier)
{
    MAssociate::WaveformIdentifier waveid;
    std::string network = "UU";
    std::string station = "CAPU";
    std::string channel = "EHZ";
    std::string location = "01";

    EXPECT_TRUE(waveid.isEmpty());

    waveid.setNetwork(network);
    EXPECT_FALSE(waveid.isEmpty());
    EXPECT_EQ(waveid.getNetwork(), network);
    waveid.setStation(station);
    EXPECT_EQ(waveid.getStation(), station);
    waveid.setChannel(channel);
    EXPECT_EQ(waveid.getChannel(), channel);
    waveid.setLocationCode(location);
    EXPECT_EQ(waveid.getLocationCode(), location);

    MAssociate::WaveformIdentifier waveidCopy(waveid);
    EXPECT_EQ(waveidCopy.getNetwork(), network);
    EXPECT_EQ(waveidCopy.getStation(), station);
    EXPECT_EQ(waveidCopy.getChannel(), channel);
    EXPECT_EQ(waveidCopy.getLocationCode(), location);

    EXPECT_TRUE(waveidCopy == waveid);

    waveid.clear();
    EXPECT_TRUE(waveidCopy != waveid);
    EXPECT_TRUE(waveid.isEmpty());

    MAssociate::WaveformIdentifier waveid2(network, station, channel, location);
    EXPECT_TRUE(waveid2, waveidCopy);
}
}

