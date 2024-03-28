#include <cmath>
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
    std::string mPhaseHint;
    uint64_t mIdentifier{0};
    std::chrono::microseconds mTime{0};
    double mStandardError{1};
    double mFirstMotionWeight{1};
    Arrival::FirstMotion mFirstMotion{Arrival::FirstMotion::Unknown};
    bool mHaveIdentifier{false};
    bool mHavePickTime{false};
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
    pImpl = std::make_unique<PickImpl> ();
}

/// Set/get phase name
void Pick::setPhaseHint(const Pick::PhaseHint phaseHint) noexcept
{
    if (phaseHint == Pick::PhaseHint::P)
    {
        setPhaseHint("P");
    }
    else if (phaseHint == Pick::PhaseHint::S)
    {
        setPhaseHint("S");
    }
    else if (phaseHint == Pick::PhaseHint::Unknown)
    {
        setPhaseHint("Unknown");
    }
    else if (phaseHint == Pick::PhaseHint::Noise)
    {
        setPhaseHint("Noise");
    }
#ifndef NDEBUG
    else
    {
        assert(false);
    }
#endif
}

void Pick::setPhaseHint(const std::string &phaseHint)
{
    if (phaseHint.empty()){throw std::invalid_argument("Phase hint is blank");}
    pImpl->mPhaseHint = phaseHint;
}

std::string Pick::getPhaseHint() const
{
    if (!havePhaseHint()){throw std::runtime_error("Phase hint not set.");}
    return pImpl->mPhaseHint;
} 

bool Pick::havePhaseHint() const noexcept
{
    return !pImpl->mPhaseHint.empty();
}

/// Get/set pick time
void Pick::setTime(const double time) noexcept
{
    auto iTime = static_cast<int64_t> (std::round(time*1.e6));
    setTime(std::chrono::microseconds {iTime});
}

void Pick::setTime(const std::chrono::microseconds &time) noexcept
{
    pImpl->mTime = time;
    pImpl->mHavePickTime = true;
}

std::chrono::microseconds Pick::getTime() const
{
    if (!haveTime()){throw std::invalid_argument("Pick time not set");}
    return pImpl->mTime;
}

bool Pick::haveTime() const noexcept
{
    return pImpl->mHavePickTime;
}

/// First motion
void Pick::setFirstMotion(const Arrival::FirstMotion firstMotion) noexcept
{
    pImpl->mFirstMotion = firstMotion;
}

MAssociate::Arrival::FirstMotion Pick::getFirstMotion() const noexcept
{
    return pImpl->mFirstMotion;
}

/// First motion weight
void Pick::setFirstMotionWeight(const double weight)
{
    if (weight < 0)
    {
        throw std::invalid_argument("Weight must be positive");
    }
    pImpl->mFirstMotionWeight = weight;
}

double Pick::getFirstMotionWeight() const noexcept
{
    return pImpl->mFirstMotionWeight;
}

/// Standard error
void Pick::setStandardError(const double standardError)
{
    if (standardError <= 0)
    {
        throw std::invalid_argument("Standard error must be positive");
    }
    pImpl->mStandardError = standardError;
}

double Pick::getStandardError() const noexcept
{
    return pImpl->mStandardError;
}

double Pick::getWeight() const noexcept
{
    return 1/pImpl->mStandardError;
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
        result = result + "   Time: "
               + std::to_string(pick.getTime().count()*1.e-6) + "\n";
    }
    result = result + "   Standard Error: "
                    + std::to_string(pick.getStandardError())
                    + " (s)\n";
    result = result + "   Weight: " + std::to_string(pick.getWeight())
                    + " (1/s)\n";
    if (pick.havePhaseHint())
    {
        result = result + "   Phase Hint: " + pick.getPhaseHint() + "\n";
    }
    auto firstMotion = pick.getFirstMotion();
    if (firstMotion == MAssociate::Arrival::FirstMotion::Up)
    {
        result = result + "   First Motion: Up\n";
    }
    else if (firstMotion == MAssociate::Arrival::FirstMotion::Down)
    {
        result = result + "   First Motion: Down\n";
    }
    else
    {
        result = result + "   First Motion: Unknown\n";
    }
    if (pick.haveIdentifier())
    {
        result = result + "   Identifier: "
               + std::to_string(pick.getIdentifier()) + "\n";
    }
    return os << result;
}
