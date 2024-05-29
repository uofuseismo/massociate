#include <cmath>
#include <chrono>
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
    std::chrono::microseconds mOriginTime{0};
    uint64_t mIdentifier{0};
    double mLatitude{0};
    double mLongitude{0};
    double mDepth{0};
    //double mX{0};
    //double mY{0};
    //double mZ{0};
    int mNumberOfPArrivals{-1};
    Type mType{Type::Event};
    bool mHaveIdentifier{false};
    bool mHaveLatitude{false};
    bool mHaveLongitude{false};
    bool mHaveDepth{false};
    bool mHaveOriginTime{false};
    //bool mHaveX{false};
    //bool mHaveY{false};
    //bool mHaveZ{false};
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
    pImpl = std::make_unique<EventImpl> ();
}

/*
/// x position
void Event::setXPosition(const double x) noexcept
{
    pImpl->mX = x;
    pImpl->mHaveX = true;
}

double Event::getXPosition() const
{
    if (!haveXPosition()){throw std::runtime_error("x position not set");}
    return pImpl->mX;
}

bool Event::haveXPosition() const noexcept
{
    return pImpl->mHaveX;
}

// y position
void Event::setYPosition(const double y) noexcept
{
    pImpl->mY = y;
    pImpl->mHaveY = true;
}

double Event::getYPosition() const
{
    if (!haveYPosition()){throw std::runtime_error("y position not set");}
    return pImpl->mY;
}

bool Event::haveYPosition() const noexcept
{
    return pImpl->mHaveY;
}

// z position
void Event::setZPosition(const double z) noexcept
{
    pImpl->mZ = z;
    pImpl->mHaveZ = true;
}

double Event::getZPosition() const
{
    if (!haveZPosition()){throw std::runtime_error("z position not set");}
    return pImpl->mZ;
}

bool Event::haveZPosition() const noexcept
{
    return pImpl->mHaveZ;
}
*/

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
    auto iOriginTime = static_cast<int64_t> (std::round(originTime*1.e6));
    setOriginTime(std::chrono::microseconds {iOriginTime});
}

void Event::setOriginTime(const std::chrono::microseconds &originTime) noexcept
{
    pImpl->mOriginTime = originTime;
    pImpl->mHaveOriginTime = true;
}

std::chrono::microseconds Event::getOriginTime() const
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
    pImpl->mNumberOfPArrivals =-1;
    if (!arrival.haveTime())
    {
        throw std::invalid_argument("Arrival time not set");
    }
    if (!arrival.haveWaveformIdentifier())
    {
        throw std::invalid_argument("Waveform identifier not set");
    }
    if (!arrival.havePhase())
    {
        throw std::invalid_argument("Arrival phase name not set");
    }
    if (!arrival.haveIdentifier())
    {
        throw std::invalid_argument("Arrival identifier not set");
    }
    // Make sure adding this arrival makes sense 
    auto nArrivals = getNumberOfArrivals();
    auto index = canAddArrival(arrival, true);
    if (index < 0)
    {
        if (index ==-1){return;}
        if (index ==-2)
        {
            throw std::invalid_argument("P can't arrive after S");
        }
        throw std::invalid_argument("S can't arrive before P");
    }
    if (index < nArrivals)
    {
        std::cerr << "Overwriting arrival" << std::endl;
        pImpl->mArrivals[index] = arrival;
        return;
    }
    pImpl->mArrivals.push_back(arrival);
    std::sort(pImpl->mArrivals.begin(), 
              pImpl->mArrivals.end(),
              [](const Arrival &lhs, const Arrival &rhs)
              {
                  return lhs.getTime() < rhs.getTime();
              });
}

int Event::canAddArrival(const Arrival &arrival,
                         const bool overwriteIfExists) const noexcept
{
    if (!arrival.haveTime()){return -1;}
    if (!arrival.haveWaveformIdentifier()){return -1;}
    if (!arrival.havePhase()){return -1;}
    if (!arrival.haveIdentifier()){return -1;}
    // Arrival can't come in before origin time
    auto newTime = arrival.getTime();
    if (haveOriginTime())
    {
        if (newTime < getOriginTime()){return -1;} 
    }
    // Get the other arrival information 
    auto nArrivals = getNumberOfArrivals();
    if (nArrivals == 0){return nArrivals;}
    auto newWaveID = arrival.getWaveformIdentifier();
    auto newPhase = arrival.getPhase();
    for (int i=0; i<nArrivals; ++i)
    {
        auto waveid = pImpl->mArrivals[i].getWaveformIdentifier();
        auto phase = pImpl->mArrivals[i].getPhase();
        auto time = pImpl->mArrivals[i].getTime();
        if (waveid == newWaveID)
        {
            if (newPhase == phase)
            {
                if (overwriteIfExists){return i;}
                return -1;
            }
            if (newPhase == "P" && phase == "S")
            {
                if (newTime >= time){return -2;}
            }
            if (newPhase == "S" && phase == "P")
            {
                if (newTime <= time){return -3;}
           }
        }
    }
    return nArrivals;
}

const std::vector<Arrival> &Event::getArrivalsReference() const
{
    return *&pImpl->mArrivals;
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

int Event::getNumberOfPArrivals() const noexcept
{
    // Easy case - nothing to do
    if (getNumberOfArrivals() == 0){return 0;}
    // If not already computed get the number of P arrivals
    if (pImpl->mNumberOfPArrivals < 0)
    {
        pImpl->mNumberOfPArrivals = 0;
        for (const auto &arrival : pImpl->mArrivals)
        {
            auto phase = arrival.getPhase();
            if (phase == "P")
            {
                pImpl->mNumberOfPArrivals = pImpl->mNumberOfPArrivals + 1;
            }
        }
    }
    return pImpl->mNumberOfPArrivals;
}

bool Event::haveArrival(const uint64_t arrivalIdentifier) const noexcept
{
    auto result = false;
    for (const auto &arrival : pImpl->mArrivals)
    {
        if (arrival.getIdentifier() == arrivalIdentifier){return true;}
    }
    return result;
}

/// Event type
void Event::setType(Type type) noexcept
{
    pImpl->mType = type;
}

Event::Type Event::getType() const noexcept
{
    return pImpl->mType;
}
