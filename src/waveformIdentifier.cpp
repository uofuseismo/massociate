#include <string>
#include "massociate/waveformIdentifier.hpp"

using namespace MAssociate;

class WaveformIdentifier::WaveformIdentifierImpl
{
public:
    WaveformIdentifierImpl()
    {
        mNetwork.reserve(16);
        mStation.reserve(16);
        mChannel.reserve(16);
        mLocation.reserve(16);
    }
    void clear() noexcept
    {
        mNetwork.clear();
        mStation.clear();
        mChannel.clear();
        mLocation.clear();
    } 
    std::string mNetwork;
    std::string mStation;
    std::string mChannel;
    std::string mLocation;
};

/// Constructors
WaveformIdentifier::WaveformIdentifier() :
    pImpl(std::make_unique<WaveformIdentifierImpl> ())
{
}

WaveformIdentifier::WaveformIdentifier(const WaveformIdentifier &waveid)
{
    *this = waveid;
}

[[maybe_unused]]
WaveformIdentifier::WaveformIdentifier(WaveformIdentifier &&waveid) noexcept
{
    *this = std::move(waveid);
}

/// Operators
WaveformIdentifier&
WaveformIdentifier::operator=(const WaveformIdentifier &waveid)
{
    if (&waveid == this){return *this;}
    pImpl = std::make_unique<WaveformIdentifierImpl> (*waveid.pImpl);
    return *this;
}

WaveformIdentifier&
WaveformIdentifier::operator=(WaveformIdentifier &&waveid) noexcept
{
    if (&waveid == this){return *this;}
    pImpl = std::move(waveid.pImpl);
    return *this;
}

bool MAssociate::operator==(const WaveformIdentifier &lhs,
                            const WaveformIdentifier &rhs) noexcept
{
    if (lhs.getNetwork() != rhs.getNetwork()){return false;}
    if (lhs.getStation() != rhs.getStation()){return false;}
    if (lhs.getChannel() != rhs.getChannel()){return false;}
    return (lhs.getLocationCode() == rhs.getLocationCode());
}

bool MAssociate::operator!=(const WaveformIdentifier &lhs,
                            const WaveformIdentifier &rhs) noexcept
{
    return !(lhs == rhs);
}

/// Destructors
WaveformIdentifier::~WaveformIdentifier() = default;

void WaveformIdentifier::clear() noexcept
{
    pImpl->clear();
}

/// Network
std::string WaveformIdentifier::getNetwork() const noexcept
{
    return pImpl->mNetwork;
}

void WaveformIdentifier::setNetwork(const std::string &str) noexcept
{
    pImpl->mNetwork = str;
}

/// Station
std::string WaveformIdentifier::getStation() const noexcept
{
    return pImpl->mStation;
}

void WaveformIdentifier::setStation(const std::string &str) noexcept
{
    pImpl->mStation = str;
}

/// Channel
std::string WaveformIdentifier::getChannel() const noexcept
{
    return pImpl->mChannel;
}

void WaveformIdentifier::setChannel(const std::string &str) noexcept
{
    pImpl->mChannel = str;
}

/// Location code
std::string WaveformIdentifier::getLocationCode() const noexcept
{
    return pImpl->mLocation;
}

void WaveformIdentifier::setLocationCode(const std::string &str) noexcept
{
    pImpl->mLocation = str;
}

bool WaveformIdentifier::isEmpty() const noexcept
{
    auto lenos = pImpl->mNetwork.size() + pImpl->mStation.size()
               + pImpl->mChannel.size() + pImpl->mLocation.size();
    return (lenos == 0);
}
