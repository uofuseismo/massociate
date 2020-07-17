#ifndef MASSOCIATOR_ASSOCIATOR_HPP
#define MASSOCIATOR_ASSOCIATOR_HPP
#include <memory>
namespace MAssociate
{
class Pick;
class AssociatorParameters;
namespace Mesh
{
template<class T> class IMesh;
}
};
namespace MAssociate
{
/*!
 * @brief This class performs the association.
 */
template<class T>
class Associator
{
public:
    /*!
     * @brief Constructor.
     */
    Associator();
    /*!
     * @brief Destructor.
     */
    ~Associator();
    /*!
     * @brief Releases all memory and resets the class.
     */
    void clear() noexcept;

    /*!
     * @brief Initializes the associator class.
     * @param parameters   The parameters from which to initialize this class.
     * @throws std::invalid_argument if required parameters are not set.
     */
    void initialize(const AssociatorParameters &parameters,
                    const Mesh::IMesh<T> &geometry);
    /*!
     * @return True indicates that the class is initialized.
     */
    [[nodiscard]] bool isInitialized() const noexcept;

    /*! @name Step 3: Set Picks
     * @{
     */
    /*!
     * @brief Adds an pick whose differential travel time will be migrated.
     * @param[in] pick   The pick whose time will be migrated.
     * @throws std::invalid_argument if the pick's network/station/phase
     *         do not correspond to an existing travel time table or the
     *         pick time is not set.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c haveTravelTimeTable().
     */
    void addPick(const Pick &pick);
    /*!
     * @return The number of picks.
     */
    [[nodiscard]] int getNumberOfPicks() const noexcept;
    /*!
     * @brief Clears the picks.
     * @note This should be used when seeking to migrate a new batch of
     *       picks using the existing travel time tables.
     */
    void clearPicks() noexcept;
    /*! @} */


    void associate();
private:
    class AssociatorImpl;
    std::unique_ptr<AssociatorImpl> pImpl;
};
}
#endif
