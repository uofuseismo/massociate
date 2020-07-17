  
cmake_minimum_required(VERSION 3.10)
project(MAssociate VERSION 0.1.0 LANGUAGES CXX)
enable_testing()

include(CheckCXXCompilerFlag)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(FindOpenMP QUIET)
if (OPENMP_FOUND)
   message("OpenMP found")
   string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS}")
endif()

set(FindDAAL_DIR ${CMAKE_SOURCE_DIR}/CMakeModules)
set(FindTBB_DIR ${CMAKE_SOURCE_DIR}/CMakeModules)
set(FindMKL_DIR ${CMAKE_SOURCE_DIR}/CMakeModules)
find_package(GTest REQUIRED)
find_package(Boost REQUIRED)
find_package(FindDAAL REQUIRED)
find_package(FindMKL REQUIRED)
find_package(FindTBB)
if (TBB_FOUND)
   message("TBB found")
   add_compile_definitions(USE_TBB)
endif()

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${Boost_INCLUDE_DIRS}
    ${DAAL_INCLUDE_DIR}
    ${TBB_INCLUDE_DIR}
    ${MKL_INCLUDE_DIR}
)
configure_file(${CMAKE_SOURCE_DIR}/include/private/version.hpp.in
               ${CMAKE_SOURCE_DIR}/include/massociate/version.hpp)
set(PUBLIC_HEADER_DIRECTORIES
    include/massociate)
set(PRIVATE_HEADER_DIRECTORIES
    include/private)

set(SRC
    src/arrival.cpp
    src/associator.cpp
    src/associatorParameters.cpp
    src/dbscan.cpp
    src/event.cpp
    src/migrate.cpp
    src/migrationParameters.cpp
    src/pageRank.cpp
    src/pick.cpp
    src/version.cpp
    src/waveformIdentifier.cpp
    src/mesh/spherical/points3d.cpp
    src/mesh/cartesian/points3d.cpp)

#cmake -DBUILD_SHARED_LIBS=YES /path/to/source
set(BUILD_SHARED_LIBS YES)
add_library(massociate ${SRC})
target_link_libraries(massociate ${DAAL_LIBRARY} ${TBB_LIBRARY})
target_include_directories(massociate PUBLIC  ${PUBLIC_HEADER_FILES})
target_include_directories(massociate PRIVATE ${PRIVATE_HEADER_FILES})
set_target_properties(massociate PROPERTIES
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)

##########################################################################################
#                                       Unit Tests                                       #
##########################################################################################
add_executable(tests
               testing/main.cpp
               testing/waveformIdentifier.cpp
               testing/arrival.cpp
               testing/mesh.cpp
               testing/event.cpp
               testing/migrate.cpp
               testing/dbscan.cpp
               testing/associator.cpp)
target_link_libraries(tests PRIVATE massociate ${GTEST_BOTH_LIBRARIES})
target_include_directories(tests PRIVATE ${GTEST_INCLUDE_DIRS})
set_target_properties(tests PROPERTIES
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)
add_test(NAME tests
         COMMAND tests)

##########################################################################################
#                                     CPACK Packaging                                    #
##########################################################################################
set(CPACK_PACKAGE_NAME "massociate")
set(CPACK_PACKAGE_VENDOR "UUSS")
set(CPACK_PACKAGE_CONTACT "bbaker@seis.utah.edu")
set(CPACK_PACKAGE_LICENSE "MIT")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A migration-based toolkit for associating seismic picks.")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
