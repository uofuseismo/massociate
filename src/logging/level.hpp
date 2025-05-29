#ifndef ULOCATOR_UMPS_LOGGING_LEVEL_HPP
#define ULOCATOR_UMPS_LOGGING_LEVEL_HPP
#include <cstdint>
namespace UMPS::Logging
{
/// @class Level "level.hpp" "umps/logging/level.hpp"
/// @brief Defines the logging level.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
/// @ingroup Logging_BaseClass
enum class Level : uint32_t
{
    None  = 0,
    Error = 1,
    Warn  = 2,
    Info  = 3,
    Debug = 4
};
}
#endif
