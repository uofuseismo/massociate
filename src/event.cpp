#include <iostream>
#include <string>
#include <vector>
#include "massociate/event.hpp"
#include "massociate/waveformIdentifier.hpp"
#include "massociate/arrival.hpp"

using namespace MAssociate;

class Event::EventImpl
{
public:
    std::vector<Arrival> mArrivals;
    uint64_t mIdentifier = 0;
    double mLatitude = 0;
    double mLongitude = 0;
    double mDepth = 0;
    double mOriginTime = 0;
    bool mHaveIdentifier = false;
    bool mHaveLatitude = false;
    bool mHaveLongitude = false;
    bool mHaveDepth = false;
    bool mHaveOriginTime = false;
};

/// C'tor
Event::Event() :
    pImpl(std::make_unique<EventImpl>())
{
}

/// Copy c'tor
Event::Event(const Event &event)
{
    *this = event;
}

/// Move c'tor
Event::Event(Event &&event) noexcept
{
    *this = std::move(event);
}

/// Copy assignment
Event& Event::operator=(const Event &event)
{
    if (&event == this){return *this;}
    pImpl = std::make_unique<EventImpl> (*event.pImpl);
    return *this;
}

/// Move assignment
Event& Event::operator=(Event &&event) noexcept
{
    if (&event == this){return *this;}
    pImpl = std::move(event.pImpl);
    return *this;
}

/// Destructor
Event::~Event() = default;

/// Reset the class
void Event::clear() noexcept
{
    pImpl->mArrivals.clear();
    pImpl->mIdentifier = 0;
    pImpl->mLatitude = 0;
    pImpl->mLongitude = 0;
    pImpl->mDepth = 0;
    pImpl->mOriginTime = 0;
    pImpl->mHaveIdentifier = false;
    pImpl->mHaveLatitude = false;
    pImpl->mHaveLongitude = false;
    pImpl->mHaveDepth = false;
    pImpl->mHaveOriginTime = false;
}

/// Latitude
void Event::setLatitude(const double latitude)
{
    if (latitude <-90 || latitude > 90)
    {
        throw std::invalid_argument("latitude = " + std::to_string(latitude)
                                  + " must be in range [-90,90]");
    }
    pImpl->mLatitude = latitude;
    pImpl->mHaveLatitude = true;
}

double Event::getLatitude() const
{
    if (!haveLatitude()){throw std::runtime_error("latitude not set");}
    return pImpl->mLatitude;
}

bool Event::haveLatitude() const noexcept
{
    return pImpl->mHaveLatitude;
}

/// Longitude
void Event::setLongitude(const double longitude)
{
    if (longitude <-540 || longitude >= 540)
    {
        throw std::invalid_argument("longitude = " + std::to_string(longitude)
                                 +  " must be in range [-540,540)");
    }
    double lon = longitude;
    while (lon < 0){lon = lon + 360;}
    while (lon >= 360){lon = lon - 360;}
    pImpl->mLongitude = lon;
    pImpl->mHaveLongitude = true;
}

double Event::getLongitude() const
{
    if (!haveLongitude()){throw std::runtime_error("longitude not set");}
    return pImpl->mLongitude;
}

bool Event::haveLongitude() const noexcept
{
    return pImpl->mHaveLongitude;
}

/// Depth
void Event::setDepth(const double depth) noexcept
{
    pImpl->mDepth = depth;
    pImpl->mHaveDepth = true;
}

double Event::getDepth() const
{
    if (!haveDepth()){throw std::runtime_error("depth not set");}
    return pImpl->mDepth;
}

bool Event::haveDepth() const noexcept
{
    return pImpl->mHaveDepth;
}

/// Origin time
void Event::setOriginTime(const double originTime) noexcept
{
    pImpl->mOriginTime = originTime;
    pImpl->mHaveOriginTime = true;
}

double Event::getOriginTime() const
{
    if (!haveOriginTime()){throw std::runtime_error("origin time not set");}
    return pImpl->mOriginTime;
}

bool Event::haveOriginTime() const noexcept
{
    return pImpl->mHaveOriginTime;
}

/// Identifier
void Event::setIdentifier(uint64_t evid) noexcept
{
    pImpl->mIdentifier = evid;
    pImpl->mHaveIdentifier = true;
}

uint64_t Event::getIdentifier() const
{
    if (!haveIdentifier()){throw std::runtime_error("identifier not set");}
    return pImpl->mIdentifier;
}

bool Event::haveIdentifier() const noexcept
{
    return pImpl->mHaveIdentifier;
}

/// Arrivals
void Event::addArrival(const Arrival &arrival)
{
    if (!arrival.haveTime())
    {
        throw std::invalid_argument("Arrival time not set");
    }
    if (!arrival.haveWaveformIdentifier())
    {
        throw std::invalid_argument("Waveform identifier not set");
    }
    if (!arrival.havePhaseName())
    {
        throw std::invalid_argument("Arrival phase name not set");
    }
    if (!arrival.haveIdentifier())
    {
        throw std::invalid_argument("Arrival identifier not set");
    }
    // Make sure adding this arrival makes sense 
    bool isNew = true;
    auto nArrival = getNumberOfArrivals();
    auto newWaveID = arrival.getWaveformIdentifier();
    auto newPhase = arrival.getPhaseName();
    auto newTime = arrival.getTime();
    for (int i=0; i<getNumberOfArrivals(); ++i)
    {
        auto waveid = pImpl->mArrivals[i].getWaveformIdentifier();
        auto phase = pImpl->mArrivals[i].getPhaseName();
        auto time = pImpl->mArrivals[i].getTime();
        if (waveid == newWaveID)
        {
            if (newPhase == phase)
            {
                std::cerr << "Overwriting arrival" << std::endl;
                pImpl->mArrivals[i] = arrival;
                isNew = false;
                break;
            }
            if (newPhase == "P" && phase == "S")
            {
                if (newTime >= time)
                {
                    throw std::invalid_argument("P can't arrive after S");
                }
            }
            if (newPhase == "S" && phase == "P")
            {
                if (newTime < time)
                {
                    throw std::invalid_argument("S can't arrive before P");
                }
           }
        }
    }
    if (isNew){pImpl->mArrivals.push_back(arrival);}
}

std::vector<Arrival> Event::getArrivals() const
{
    return pImpl->mArrivals;
}

void Event::clearArrivals() noexcept
{
    pImpl->mArrivals.clear();
}

int Event::getNumberOfArrivals() const noexcept
{
    return static_cast<int> (pImpl->mArrivals.size());
}
