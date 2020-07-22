#ifndef BLACKLIST_HPP
#define BLACKLIST_HPP
#include <string>
/// @brief These are the nodal stations that are known to have timing issues.
/// @param[in] network  The network name.
/// @param[in] station  The station name.
/// @result True indicates that the given station is in the black list.
bool isBlackListed(const std::string &network, const std::string &station)
{
    if (network == "UU")
    {
        if (station == "136"){return true;}
        if (station == "229"){return true;}
    }
    return false;
}
#endif
