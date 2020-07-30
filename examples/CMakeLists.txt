cmake_minimum_required(VERSION 3.10)
project(MAssociateExamples VERSION 0.1.0 LANGUAGES CXX)
enable_testing()

include(CheckCXXCompilerFlag)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(FindOpenMP QUIET)
if (OPENMP_FOUND)
   message("OpenMP found")
   string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS}")
endif()

set(FindDAAL_DIR ${CMAKE_SOURCE_DIR}/../CMakeModules)
set(FindTBB_DIR ${CMAKE_SOURCE_DIR}/../CMakeModules)
set(FindMKL_DIR ${CMAKE_SOURCE_DIR}/../CMakeModules)
set(FindMAssociate_DIR ${CMAKE_SOURCE_DIR}/CMakeModules)
set(FindGeographicLib_DIR ${CMAKE_SOURCE_DIR}/CMakeModules)
set(FindHDF5_DIR ${CMAKE_SOURCE_DIR}/CMakeModules)
find_package(GTest REQUIRED)
find_package(Boost REQUIRED)
#find_package(FindDAAL REQUIRED)
find_package(FindMKL REQUIRED)
find_package(FindMAssociate REQUIRED)
find_package(FindGeographicLib REQUIRED)
find_package(FindHDF5 REQUIRED)
find_package(FindTBB)
if (TBB_FOUND)
   message("TBB found")
   add_compile_definitions(USE_TBB)
endif()

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${MASSOCIATE_INCLUDE_DIR}
    ${GEOGRAPHICLIB_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIR}
    ${MKL_INCLUDE_DIR}
    ${H5_INCLUDE_DIR}
    ${SFF_INCLUDE_DIR}
)
#configure_file(${CMAKE_SOURCE_DIR}/include/private/version.hpp.in
#               ${CMAKE_SOURCE_DIR}/include/massociate/version.hpp)
#set(PUBLIC_HEADER_DIRECTORIES
#    include/massociate)
set(PRIVATE_HEADER_DIRECTORIES
    include)

set(UTILITIES_SRC utilities/bilinearInterpolation.cpp utilities/h5io.cpp)
add_executable(ttables magna/ttables.cpp ${UTILITIES_SRC})
target_link_libraries(ttables ${MASSOCIATE_LIBRARY} ${HDF5_LIBRARY} ${MKL_LIBRARY} ${GEOGRAPHICLIB_LIBRARIES})
target_include_directories(ttables PRIVATE ${PRIVATE_HEADER_FILES})
set_target_properties(ttables PROPERTIES
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)

add_executable(associate magna/associate.cpp ${UTILITIES_SRC})
target_link_libraries(associate ${SFF_LIBRARY} ${MASSOCIATE_LIBRARY} ${HDF5_LIBRARY} ${MKL_LIBRARY})
target_include_directories(associate PRIVATE ${PRIVATE_HEADER_FILES})
set_target_properties(associate PROPERTIES
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)

##########################################################################################
#                                       Unit Tests                                       #
##########################################################################################
add_executable(tests
               testing/main.cpp
               testing/utilities.cpp
               ${UTILITIES_SRC})
target_link_libraries(tests PRIVATE ${MASSOCIATE_LIBRARY} ${HDF5_LIBRARY} ${MKL_LIBRARY} ${GTEST_BOTH_LIBRARIES})
target_include_directories(tests PRIVATE ${GTEST_INCLUDE_DIRS})
set_target_properties(tests PROPERTIES
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED YES
                      CXX_EXTENSIONS NO)
add_test(NAME tests
         COMMAND tests)
