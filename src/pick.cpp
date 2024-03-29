#include <string>
#include <sstream>
#include "massociate/pick.hpp"
#include "massociate/arrival.hpp"
#include "massociate/waveformIdentifier.hpp"

using namespace MAssociate;

class Pick::PickImpl
{
public:
    MAssociate::WaveformIdentifier mWaveID;
    std::string mPhaseName;
    uint64_t mIdentifier = 0;
    double mPickTime = 0;
    double mStd = 1;
    double mStaticCorrection = 0;
    double mPolarityWeight = 1;
    MAssociate::Polarity mPolarity = MAssociate::Polarity::UNKNOWN;
    bool mHaveIdentifier = false;
    bool mHavePickTime = false;
};

/// C'tor
Pick::Pick() :
    pImpl(std::make_unique<PickImpl> ())
{
}

/// Copy c'tor
Pick::Pick(const Pick &pick)
{
    *this = pick;
}

/// Move c'tor
[[maybe_unused]]
Pick::Pick(Pick &&pick) noexcept
{
    *this = std::move(pick);
}

/// Copy assignment
Pick& Pick::operator=(const Pick &pick)
{
    if (&pick == this){return *this;}
    pImpl = std::make_unique<PickImpl> (*pick.pImpl);
    return *this;
}

/// Move assignment
Pick& Pick::operator=(Pick &&pick) noexcept
{
    if (&pick == this){return *this;}
    pImpl = std::move(pick.pImpl);
    return *this;
}

/// Destructor
Pick::~Pick() = default;

/// Resets the class
[[maybe_unused]]
void Pick::clear() noexcept
{
    pImpl->mWaveID.clear();
    pImpl->mPhaseName.clear();
    pImpl->mIdentifier = 0;
    pImpl->mPickTime = 0;
    pImpl->mStd = 1;
    pImpl->mStaticCorrection = 0;
    pImpl->mPolarityWeight = 1;
    pImpl->mPolarity = MAssociate::Polarity::UNKNOWN;
    pImpl->mHaveIdentifier = false;
    pImpl->mHavePickTime = false;
}

/// Set/get phase name
void Pick::setPhaseName(const std::string &phaseName)
{
    if (phaseName.empty()){throw std::invalid_argument("Phase name is blank");}
    pImpl->mPhaseName = phaseName;
}

std::string Pick::getPhaseName() const
{
    if (!havePhaseName()){throw std::runtime_error("Phase name not set.");}
    return pImpl->mPhaseName;
} 

bool Pick::havePhaseName() const noexcept
{
    return !pImpl->mPhaseName.empty();
}

/// Get/set pick time
void Pick::setTime(const double pickTime) noexcept
{
    pImpl->mPickTime = pickTime;
    pImpl->mHavePickTime = true;
}

double Pick::getTime() const
{
    if (!haveTime()){throw std::invalid_argument("Pick time not set");}
    return pImpl->mPickTime;
}

bool Pick::haveTime() const noexcept
{
    return pImpl->mHavePickTime;
}

/// Polarity
void Pick::setPolarity(const MAssociate::Polarity polarity) noexcept
{
    pImpl->mPolarity = polarity;
}

MAssociate::Polarity Pick::getPolarity() const noexcept
{
    return pImpl->mPolarity;
}

/// Polarity weight
void Pick::setPolarityWeight(const double weight)
{
    if (weight < 0)
    {
        throw std::invalid_argument("weight must be positive");
    }
    pImpl->mPolarityWeight = weight;
}

double Pick::getPolarityWeight() const noexcept
{
    return pImpl->mPolarityWeight;
}

/// Standard deviation
void Pick::setStandardDeviation(const double std)
{
    if (std <= 0)
    {
        throw std::invalid_argument("Standard deviation must be positive");
    }
    pImpl->mStd = std;
}

double Pick::getStandardDeviation() const noexcept
{
    return pImpl->mStd;
}

double Pick::getWeight() const noexcept
{
    return 1/pImpl->mStd;
}

/// WaveformIdentifier
void Pick::setWaveformIdentifier(const WaveformIdentifier &waveid)
{
    if (waveid.isEmpty())
    {
        throw std::invalid_argument("WaveformIdentifier is empty");
    }
    pImpl->mWaveID = waveid;
}

MAssociate::WaveformIdentifier Pick::getWaveformIdentifier() const
{
    if (!haveWaveformIdentifier())
    {
        throw std::runtime_error("WaveformIdentifier not set");
    }
    return pImpl->mWaveID;
}

bool Pick::haveWaveformIdentifier() const noexcept
{
    return !pImpl->mWaveID.isEmpty();
}

/// Identifier
void Pick::setIdentifier(uint64_t identifier) noexcept
{
    pImpl->mIdentifier = identifier;
    pImpl->mHaveIdentifier = true;
}

uint64_t Pick::getIdentifier() const
{
    if (!haveIdentifier()){throw std::runtime_error("identifier not yet set");}
    return pImpl->mIdentifier;
}

bool Pick::haveIdentifier() const noexcept
{
    return pImpl->mHaveIdentifier;
}

/// Static correction
void Pick::setStaticCorrection(const double correction) noexcept
{
    pImpl->mStaticCorrection = correction;
}

double Pick::getStaticCorrection() const noexcept
{
    return pImpl->mStaticCorrection;
}

/// Print the pick
std::ostream&
MAssociate::operator<<(std::ostream &os,
                       const MAssociate::Pick &pick)
{
    std::string result = "Pick\n";
    if (pick.haveWaveformIdentifier())
    {
        std::stringstream ss;
        ss << pick.getWaveformIdentifier();
        result = result + "   Waveform: " + ss.str() + "\n";
    } 
    if (pick.haveTime())
    {
        result = result + "   Time: " + std::to_string(pick.getTime()) + "\n";
    }
    result = result + "   Standard Deviation: "
                    + std::to_string(pick.getStandardDeviation())
                    + " (s)\n";
    result = result + "   Weight: " + std::to_string(pick.getWeight())
                    + " (1/s)\n";
    if (pick.havePhaseName())
    {
        result = result + "   Phase: " + pick.getPhaseName() + "\n";
    }
    auto polarity = pick.getPolarity();
    if (polarity == MAssociate::Polarity::COMPRESSIONAL)
    {
        result = result + "   Polarity: Compression\n";
    }
    else if (polarity == MAssociate::Polarity::DILATATIONAL)
    {
        result = result + "   Polarity: Dilatational\n";
    }
    else
    {
        result = result + "   Polarity: Unknown\n";
    }
    if (pick.haveIdentifier())
    {
        result = result + "   Identifier: "
               + std::to_string(pick.getIdentifier()) + "\n";
    }
    return os << result;
}
