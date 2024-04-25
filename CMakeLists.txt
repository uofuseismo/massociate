cmake_minimum_required(VERSION 3.16)
project(MAssociate VERSION 1.0.0 LANGUAGES CXX)
enable_testing()

option(WRAP_PYTHON "Compile the Python bindings" OFF)
#cmake -DBUILD_SHARED_LIBS=YES /path/to/source
set(BUILD_SHARED_LIBS YES)

include(CheckCXXCompilerFlag)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(THREADS_PREFER_PTHREAD_FLAG ON)
set(FindDAAL_DIR ${CMAKE_SOURCE_DIR}/cmake)
set(FindTBB_DIR ${CMAKE_SOURCE_DIR}/cmake)
set(FindMKL_DIR ${CMAKE_SOURCE_DIR}/cmake)
find_package(NLopt REQUIRED)
find_package(pagmo REQUIRED)
find_package(Catch2 3 REQUIRED)
find_package(Boost REQUIRED)
find_package(FindDAAL REQUIRED)
#find_package(FindMKL REQUIRED)
find_package(FindTBB)
find_package(Threads REQUIRED)
if (TBB_FOUND)
   message("TBB found")
   add_compile_definitions(USE_TBB)
endif()

#include_directories(
#    ${CMAKE_CURRENT_SOURCE_DIR}/include
#    ${Boost_INCLUDE_DIRS}
#    ${DAAL_INCLUDE_DIR}
#    ${TBB_INCLUDE_DIR}
#    ${MKL_INCLUDE_DIR}
#)
configure_file(${CMAKE_SOURCE_DIR}/include/private/version.hpp.in
               ${CMAKE_SOURCE_DIR}/include/massociate/version.hpp)
set(PUBLIC_HEADER_DIRECTORIES
    include/massociate)
set(PRIVATE_HEADER_DIRECTORIES
    include/private)

set(SRC
    src/arrival.cpp
    src/associator.cpp
    #src/associatorParameters.cpp
    src/dbscan.cpp
    #src/direct.cpp
    src/event.cpp
    src/optimizer.cpp
    #src/migrate.cpp
    #src/migrationParameters.cpp
    #src/pageRank.cpp
    src/migrator.cpp
    src/pick.cpp
    src/version.cpp
    src/waveformIdentifier.cpp
    #src/mesh/spherical/points3d.cpp
    #src/mesh/cartesian/points3d.cpp
)

if (${pagmo_FOUND})
   add_compile_definitions(USE_PAGMO)
   set(SRC ${SRC} src/particleSwarm.cpp)
endif()
if (${NLopt_FOUND})
   add_compile_definitions(USE_NLOPT)
   set(SRC ${SRC} src/dividedRectangles.cpp)
endif()

add_library(massociate ${SRC})
target_link_libraries(massociate
                      PRIVATE ${DAAL_LIBRARY} ${MKL_LIBRARY} ${TBB_LIBRARY} ${UMPS_LIBRARY}
                              ${uLocator_LIBRARY} Threads::Threads)
target_include_directories(massociate PUBLIC  ${CMAKE_SOURCE_DIR}/include
                                      PRIVATE ${CMAKE_SOURCE_DIR}/include
                                              ${NLOPT_INCLUDE_DIRS}
                                              ${DAAL_INCLUDE_DIR}
                                              ${UMPS_INCLUDE_DIR}
                                              ${uLocator_INCLUDE_DIR})
if (${pagmo_FOUND})
   target_include_directories(massociate PRIVATE Pagmo::pagmo)
   target_link_libraries(massociate PRIVATE Pagmo::pagmo)
endif()
if (${NLopt_FOUND})
   target_include_directories(massociate PRIVATE NLopt::nlopt)
   target_link_libraries(massociate PRIVATE NLopt::nlopt)
endif() 
set_target_properties(massociate PROPERTIES
                      CXX_STANDARD 20
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)

##########################################################################################
#                                       Unit Tests                                       #
##########################################################################################
add_executable(tests
               #testing/main.cpp
               testing/waveformIdentifier.cpp
               testing/pick.cpp
               testing/arrival.cpp
#               testing/mesh.cpp
               testing/event.cpp
               testing/migrate.cpp
               testing/dbscan.cpp
               testing/associator.cpp
               )
target_link_libraries(tests PRIVATE massociate ${uLocator_LIBRARY} ${UMPS_LIBRARY} Catch2::Catch2 Catch2::Catch2WithMain)
target_include_directories(tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include
                                         ${CMAKE_CURRENT_SOURCE_DIR}/src
                                         ${UMPS_INCLUDE_DIR}
                                         ${uLocator_INCLUDE_DIR}
                                         Catch2::Catch2)
set_target_properties(tests PROPERTIES
                      CXX_STANDARD 20
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)
add_test(NAME tests
         COMMAND tests --rng-seed=2323432)

###############################################################################
#                                 Installation                                #
###############################################################################
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
configure_package_config_file(
    cmake/MAssociateConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION "${version}"
    COMPATIBILITY AnyNewerVersion
)

#if (MASSOCIATE_WRAP_PYTHON)
#   install(TARGETS massociate pymassociate
#           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
#           LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#           ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#           PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
#else()
   install(TARGETS massociate
           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
           LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
           ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
           PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
#endif()

##########################################################################################
#                                     CPACK Packaging                                    #
##########################################################################################
set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
set(CPACK_PACKAGE_VENDOR "UUSS")
set(CPACK_PACKAGE_CONTACT "ben.baker@utah.edu")
set(CPACK_PACKAGE_LICENSE "MIT")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "An earthquake assocation utility used at UUSS.")
set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})
set(CPACK_VERBATIM_VARIABLES YES)
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_COMPONENTS_ALL libraries headers)
if (WIN32)
    set(CPACK_GENERATOR ZIP WIX)
elseif (APPLE)
    set(CPACK_GENERATOR TGZ productbuild)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CPACK_GENERATOR TGZ RPM)
else()
    set(CPACK_GENERATOR TGZ)
endif()
set(CPACK_SOURCE_IGNORE_FILES
  /\\.git/
  \\.swp
  \\.orig
  /CMakeLists\\.txt\\.user
  /private/
)
include(CPack) # Put this last!

