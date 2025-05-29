#ifndef ULOCATOR_UMPS_LOGGING_STANDARD_OUT_HPP
#define ULOCATOR_UMPS_LOGGING_STANDARD_OUT_HPP
#include <memory>
#include "log.hpp"
namespace UMPS::Logging
{
/// StandardOut "stdout.hpp" "umps/logging/stdout.hpp"
/// @brief A naive logger that prints all messages to standard out and,
///        in the case of errors, standard error.
/// @note This class is not recommended for typical use.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
/// @ingroup Logging_Loggers
class StandardOut : public ILog
{
public:
    /// @brief Constructor.
    explicit StandardOut(Logging::Level level = Logging::Level::Info);
    /// @brief Copy constructor.
    StandardOut(const StandardOut &logger);
    /// @brief Move constructor.
    StandardOut(StandardOut &&logger) noexcept;

    /// @brief Copy assignment operator.
    /// @param[in] logger   The logging class to copy to this.
    /// @result A deep copy of the logger class.
    StandardOut &operator=(const StandardOut &logger); 
    /// @brief Move assignment operator.
    /// @param[in,out] logger  The logging class whose memory will be moved to
    ///                        this.  On exit, logger's behavior is undefined.
    /// @result The memory from logger moved to this.
    StandardOut &operator=(StandardOut &&logger) noexcept;

    /// @brief Destructor.
    ~StandardOut() override;

    /// @brief Sets the logging level.
    /// @param[in] level  The logging level.
    void setLevel(Level level) noexcept;
    /// @result Gets the logging level.
    [[nodiscard]] Level getLevel() const noexcept override;

    /// @brief Writes an error message.
    /// @note This requires \c getLevel() be >= Level::ERROR.
    void error(const std::string &message) override;
    /// @brief Writes a warning message.
    /// @note This requires \c getLevel() be >= Level::WARN.
    void warn(const std::string &message) override;
    /// @brief Writes an info message.
    /// @note This requires \c getLevel() be >= least Level::INFO.
    void info(const std::string &message) override;
    /// @brief Writes a debug message.
    /// @note This requires \c getLevel() be >= Level::DEBUG.
    void debug(const std::string &message) override;
private:
    class StandardOutImpl;
    std::unique_ptr<StandardOutImpl> pImpl; 
};
}
#endif
