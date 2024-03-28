#include <cmath>
#include <string>
#include <chrono>
#include <stdexcept>
#include "massociate/arrival.hpp"
#include "massociate/pick.hpp"
#include "massociate/waveformIdentifier.hpp"

using namespace MAssociate;

class Arrival::ArrivalImpl
{
public:
    MAssociate::WaveformIdentifier mWaveIdentifier;
    std::string mPhase;
    std::chrono::microseconds mArrivalTime{0};
    uint64_t mIdentifier{0};
    double mStandardError{1};
    double mTravelTime{-1};
    double mFirstMotionWeight{1};
    Arrival::FirstMotion mFirstMotion{Arrival::FirstMotion::Unknown};
    bool mHaveIdentifier{false};
    bool mHaveArrivalTime{false};
};

/// Constructor
Arrival::Arrival() :
    pImpl(std::make_unique<ArrivalImpl> ())
{
}

/// Copy constructor
Arrival::Arrival(const Arrival &arrival)
{
    *this = arrival;
}

Arrival::Arrival(const Pick &pick)
{
    *this = pick;
}

/// Move constructor
[[maybe_unused]]
Arrival::Arrival(Arrival &&arrival) noexcept
{
    *this = std::move(arrival);
}

/// Copy assignment
Arrival& Arrival::operator=(const Arrival &arrival)
{
    if (&arrival == this){return *this;}
    pImpl = std::make_unique<ArrivalImpl> (*arrival.pImpl);
    return *this;
}


/// Copy
Arrival& Arrival::operator=(const Pick &pick)
{
    Arrival arrival;
    if (pick.haveWaveformIdentifier())
    {
        arrival.setWaveformIdentifier(pick.getWaveformIdentifier());
    }
    if (pick.havePhaseHint())
    {
        auto phaseHint = pick.getPhaseHint();
        if (phaseHint == "P" || phaseHint == "S")
        {
            arrival.setPhase(phaseHint);
        }
    }
    if (pick.haveIdentifier()){arrival.setIdentifier(pick.getIdentifier());}
    if (pick.haveTime()){arrival.setTime(pick.getTime());}
    arrival.setStandardError(pick.getStandardError());
    arrival.setFirstMotion(pick.getFirstMotion());
    arrival.setFirstMotionWeight(pick.getFirstMotionWeight());
    *this = arrival;
    return *this; 
} 

/// Move assignment
Arrival& Arrival::operator=(Arrival &&arrival) noexcept
{
    if (&arrival == this){return *this;}
    pImpl = std::move(arrival.pImpl);
    return *this;
}

/// Destructor
Arrival::~Arrival() = default;

/// Resets the class
[[maybe_unused]]
void Arrival::clear() noexcept
{
    pImpl = std::make_unique<ArrivalImpl> ();
}

/// Set/get phase name
void Arrival::setPhase(const Arrival::Phase phase) noexcept
{
    if (phase == Arrival::Phase::P)
    {
        setPhase("P");
    }
    else if (phase == Arrival::Phase::S)
    {
        setPhase("S");
    }
#ifndef NDEBUG
    else
    {
        assert(false);
    }
#endif
}

void Arrival::setPhase(const std::string &phase)
{
    if (phase.empty()){throw std::invalid_argument("Phase name is blank");}
    pImpl->mPhase = phase;
}

std::string Arrival::getPhase() const
{
    if (!havePhase()){throw std::runtime_error("Phase not set.");}
    return pImpl->mPhase;
} 

bool Arrival::havePhase() const noexcept
{
    return !pImpl->mPhase.empty();
}

/// Get/set arrival time
void Arrival::setTime(const double time) noexcept
{
    auto iTime = static_cast<int64_t> (std::round(time*1.e6));
    setTime(std::chrono::microseconds {iTime});
}

void Arrival::setTime(const std::chrono::microseconds &arrivalTime) noexcept
{
    pImpl->mArrivalTime = arrivalTime;
    pImpl->mHaveArrivalTime = true;
}

std::chrono::microseconds Arrival::getTime() const
{
    if (!haveTime()){throw std::invalid_argument("Arrival time not set");}
    return pImpl->mArrivalTime;
}

bool Arrival::haveTime() const noexcept
{
    return pImpl->mHaveArrivalTime;
}

/// First motion
void Arrival::setFirstMotion(const Arrival::FirstMotion firstMotion) noexcept
{
    pImpl->mFirstMotion = firstMotion;
}

Arrival::FirstMotion Arrival::getFirstMotion() const noexcept
{
    return pImpl->mFirstMotion;
}

/// First motion weight
void Arrival::setFirstMotionWeight(const double weight)
{
    if (weight < 0)
    {
        throw std::invalid_argument("weight must be positive");
    }
    pImpl->mFirstMotionWeight = weight;
}

double Arrival::getFirstMotionWeight() const noexcept
{
    return pImpl->mFirstMotionWeight;
}

/// Standard error
void Arrival::setStandardError(const double standardError)
{
    if (standardError <= 0)
    {
        throw std::invalid_argument("Standard error must be positive");
    }
    pImpl->mStandardError = standardError;
}

double Arrival::getStandardError() const noexcept
{
    return pImpl->mStandardError;
}

double Arrival::getWeight() const noexcept
{
    return 1.0/pImpl->mStandardError;
}

/// Waveform identifier 
void Arrival::setWaveformIdentifier(const WaveformIdentifier &waveid)
{
    if (waveid.isEmpty())
    {
        throw std::invalid_argument("waveform identifier is empty");
    }
    pImpl->mWaveIdentifier= waveid;
}

MAssociate::WaveformIdentifier Arrival::getWaveformIdentifier() const
{
    if (!haveWaveformIdentifier())
    {
        throw std::runtime_error("waveform identifier not set");
    }
    return pImpl->mWaveIdentifier;
}

bool Arrival::haveWaveformIdentifier() const noexcept
{
    return !pImpl->mWaveIdentifier.isEmpty();
}

/// Identifier
void Arrival::setIdentifier(uint64_t identifier) noexcept
{
    pImpl->mIdentifier = identifier;
    pImpl->mHaveIdentifier = true;
}

uint64_t Arrival::getIdentifier() const
{
    if (!haveIdentifier()){throw std::runtime_error("identifier not yet set");}
    return pImpl->mIdentifier;
}

bool Arrival::haveIdentifier() const noexcept
{
    return pImpl->mHaveIdentifier;
}

/// Travel time
void Arrival::setTravelTime(const double travelTime)
{
    if (travelTime < 0)
    {
        throw std::invalid_argument("Travel time must be positive");
    }
    pImpl->mTravelTime = travelTime;
}

double Arrival::getTravelTime() const
{
    return pImpl->mTravelTime;
}

bool Arrival::haveTravelTime() const noexcept
{
    return (pImpl->mTravelTime >=0);
}
