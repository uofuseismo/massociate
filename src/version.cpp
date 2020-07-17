#include <string>
#include "massociate/version.hpp"

using namespace MAssociate;

int Version::getMajor() noexcept
{
    return MASSOCIATE_MAJOR;
}

int Version::getMinor() noexcept
{
    return MASSOCIATE_MINOR;
}

int Version::getPatch() noexcept
{
    return MASSOCIATE_PATCH;
}

[[maybe_unused]] bool Version::isAtLeast(const int major, const int minor,
                        const int patch) noexcept
{
    if (MASSOCIATE_MAJOR < major){return false;}
    if (MASSOCIATE_MAJOR > major){return true;}
    if (MASSOCIATE_MINOR < minor){return false;}
    if (MASSOCIATE_MINOR > minor){return true;}
    return MASSOCIATE_PATCH >= patch;
}

[[maybe_unused]] std::string Version::getVersion() noexcept
{
    std::string version(std::to_string(getMajor()) + "."
                      + std::to_string(getMinor()) + "."
                      + std::to_string(getPatch()));
    return version;
}
