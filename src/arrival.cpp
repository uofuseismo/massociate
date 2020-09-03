#include <string>
#include "massociate/arrival.hpp"
#include "massociate/pick.hpp"
#include "massociate/waveformIdentifier.hpp"

using namespace MAssociate;

class Arrival::ArrivalImpl
{
public:
    MAssociate::WaveformIdentifier mWaveID;
    std::string mPhaseName;
    uint64_t mIdentifier = 0;
    double mArrivalTime = 0;
    double mStd = 1;
    double mTravelTime =-1;
    double mStaticCorrection = 0;
    double mPolarityWeight = 1;
    MAssociate::Polarity mPolarity = MAssociate::Polarity::UNKNOWN;
    bool mHaveIdentifier = false;
    bool mHaveArrivalTime = false;
};

/// C'tor
Arrival::Arrival() :
    pImpl(std::make_unique<ArrivalImpl> ())
{
}

/// Copy c'tor
Arrival::Arrival(const Arrival &arrival)
{
    *this = arrival;
}

Arrival::Arrival(const Pick &pick)
{
    *this = pick;
}

/// Move c'tor
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
    if (pick.havePhaseName()){arrival.setPhaseName(pick.getPhaseName());}
    if (pick.haveIdentifier()){arrival.setIdentifier(pick.getIdentifier());}
    if (pick.haveTime()){arrival.setTime(pick.getTime());}
    arrival.setStandardDeviation(pick.getStandardDeviation());
    arrival.setStaticCorrection(pick.getStaticCorrection());
    arrival.setPolarity(pick.getPolarity());
    arrival.setPolarityWeight(pick.getPolarityWeight());
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
    pImpl->mWaveID.clear();
    pImpl->mPhaseName.clear();
    pImpl->mIdentifier = 0;
    pImpl->mArrivalTime = 0;
    pImpl->mStd = 1;
    pImpl->mStaticCorrection = 0;
    pImpl->mTravelTime =-1;
    pImpl->mPolarityWeight = 1;
    pImpl->mPolarity = MAssociate::Polarity::UNKNOWN;
    pImpl->mHaveIdentifier = false;
    pImpl->mHaveArrivalTime = false;
}

/// Set/get phase name
void Arrival::setPhaseName(const std::string &phaseName)
{
    if (phaseName.empty()){throw std::invalid_argument("Phase name is blank");}
    pImpl->mPhaseName = phaseName;
}

std::string Arrival::getPhaseName() const
{
    if (!havePhaseName()){throw std::runtime_error("Phase name not set.");}
    return pImpl->mPhaseName;
} 

bool Arrival::havePhaseName() const noexcept
{
    return !pImpl->mPhaseName.empty();
}

/// Get/set arrival time
void Arrival::setTime(const double arrivalTime) noexcept
{
    pImpl->mArrivalTime = arrivalTime;
    pImpl->mHaveArrivalTime = true;
}

double Arrival::getTime() const
{
    if (!haveTime()){throw std::invalid_argument("Arrival time not set");}
    return pImpl->mArrivalTime;
}

bool Arrival::haveTime() const noexcept
{
    return pImpl->mHaveArrivalTime;
}

/// Polarity
void Arrival::setPolarity(const MAssociate::Polarity polarity) noexcept
{
    pImpl->mPolarity = polarity;
}

MAssociate::Polarity Arrival::getPolarity() const noexcept
{
    return pImpl->mPolarity;
}

/// Polarity weight
void Arrival::setPolarityWeight(const double weight)
{
    if (weight < 0)
    {
        throw std::invalid_argument("weight must be positive");
    }
    pImpl->mPolarityWeight = weight;
}

double Arrival::getPolarityWeight() const noexcept
{
    return pImpl->mPolarityWeight;
}

/// Standard deviation
void Arrival::setStandardDeviation(const double std)
{
    if (std <= 0)
    {
        throw std::invalid_argument("Standard deviation must be positive");
    }
    pImpl->mStd = std;
}

double Arrival::getStandardDeviation() const noexcept
{
    return pImpl->mStd;
}

double Arrival::getWeight() const noexcept
{
    return 1/pImpl->mStd;
}

/// Waveform identifier 
void Arrival::setWaveformIdentifier(const WaveformIdentifier &waveid)
{
    if (waveid.isEmpty())
    {
        throw std::invalid_argument("waveform identifier is empty");
    }
    pImpl->mWaveID = waveid;
}

MAssociate::WaveformIdentifier Arrival::getWaveformIdentifier() const
{
    if (!haveWaveformIdentifier())
    {
        throw std::runtime_error("waveform identifier not set");
    }
    return pImpl->mWaveID;
}

bool Arrival::haveWaveformIdentifier() const noexcept
{
    return !pImpl->mWaveID.isEmpty();
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

/// Static correction
void Arrival::setStaticCorrection(const double correction) noexcept
{
    pImpl->mStaticCorrection = correction;
}

double Arrival::getStaticCorrection() const noexcept
{
    return pImpl->mStaticCorrection;
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
