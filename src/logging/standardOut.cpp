#include <iostream>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_sinks.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include "standardOut.hpp"

using namespace UMPS::Logging;

class StandardOut::StandardOutImpl
{
public:
    StandardOutImpl() :
        mStandardOutSink(std::make_shared<spdlog::sinks::stdout_color_sink_mt> ()),
        mLogger(std::make_shared<spdlog::logger> ("stdout", mStandardOutSink))
    {
//        mLogger = spdlog::create("stdout", mStandardOutSink);
    }
    //StandardOutImpl() :
    //    mSink(std::make_shared<spdlog::sinks::stdout_color_sink_mt>() )
    //{   
    //}
    //std::shared_ptr<spdlog::sinks::stdout_color_sink_mt> mSink;
    Level mLevel = Level::Info;
    std::shared_ptr<spdlog::sinks::stdout_color_sink_mt> mStandardOutSink{nullptr};
    std::shared_ptr<spdlog::logger> mLogger = nullptr;
};

/// C'tor
StandardOut::StandardOut(const Level level) :
    pImpl(std::make_unique<StandardOutImpl> ())
{
    setLevel(level);
}

/// Copy c'tor
StandardOut::StandardOut(const StandardOut &logger)
{
    *this = logger;
}

/// Move c'tor
StandardOut::StandardOut(StandardOut &&logger) noexcept
{
    *this = std::move(logger);
}

/// Copy assignment
StandardOut& StandardOut::operator=(const StandardOut &logger)
{
    if (&logger == this){return *this;}
    pImpl = std::make_unique<StandardOutImpl> (*logger.pImpl);
    return *this;
}

/// Move assignment
StandardOut& StandardOut::operator=(StandardOut &&logger) noexcept
{
    if (&logger == this){return *this;}
    pImpl = std::move(logger.pImpl);
    return *this;
}

/// Destructor
StandardOut::~StandardOut() = default;

/// Level
void StandardOut::setLevel(const Level level) noexcept
{
    pImpl->mLevel = level;
    if (level == Level::Error)
    {
        pImpl->mLogger->set_level(spdlog::level::err);
    }
    else if (level == Level::Warn)
    {
        pImpl->mLogger->set_level(spdlog::level::warn);
    }
    else if (level == Level::Info)
    {
        pImpl->mLogger->set_level(spdlog::level::info);
    }
    else if (level == Level::Debug)
    {
        pImpl->mLogger->set_level(spdlog::level::debug);
    }
    else
    {
        pImpl->mLogger->set_level(spdlog::level::off);
    }
}

Level StandardOut::getLevel() const noexcept
{   
    return pImpl->mLevel;
}

/// Info
void StandardOut::info(const std::string &message)
{
    if (pImpl->mLevel >= Level::Info)
    {
        //std::cout << message << std::endl;
        //spdlog::info(message);
        pImpl->mLogger->info(message);
    }
}

/// Warn
void StandardOut::warn(const std::string &message)
{
    if (pImpl->mLevel >= Level::Warn)
    {
        //std::cout << message << std::endl;
        //spdlog::warn(message);
        pImpl->mLogger->warn(message);
    }
}

/// Error
void StandardOut::error(const std::string &message)
{
    if (pImpl->mLevel >= Level::Error)
    {
        //std::cerr << message << std::endl;
        //spdlog::error(message);
        pImpl->mLogger->error(message);
    }
}

/// Debug
void StandardOut::debug(const std::string &message)
{
    if (pImpl->mLevel >= Level::Debug)
    {
        //std::cerr << message << std::endl;
        //spdlog::debug(message);
        pImpl->mLogger->debug(message);
    }
}
