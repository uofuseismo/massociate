# Already in cache, be silent
if (GEOGRAPHICLIB_INCLUDE_DIR AND GEOGRAPHICLIB_LIBRARY)
   set (GEOGRAPHICLIB_FIND_QUIETLY TRUE)
endif()

find_path(GEOGRAPHICLIB_INCLUDE_DIR
          NAMES GeographicLib
          HINTS $ENV{GEOGRAPHICLIB_ROOT}/include
                /usr/local/include
                /usr/include)
find_library(GEOGRAPHICLIB_LIBRARIES
             NAMES Geographic
             PATHS $ENV{GEOGRAPHICLIB_ROOT}/lib
                   $ENV{GEOGRAPHICLIB_ROOT}/lib64
                   /usr/local/lib
                   /usr/local/lib64)


# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GEOGRAPHICLIB DEFAULT_MSG GEOGRAPHICLIB_LIBRARIES GEOGRAPHICLIB_INCLUDE_DIR)
mark_as_advanced(GEOGRAPHICLIB_INCLUDE_DIR GeographicLib_LIBRARIES)
