# Already in cache, be silent
if (MASSOCIATE_INCLUDE_DIR AND MASSOCIATE_LIBRARY)
   set (MASSOCIATE_FIND_QUIETLY TRUE)
endif()

set(MASSOCIATE "massociate")

find_path(MASSOCIATE_INCLUDE_DIR
          NAMES massociate
          HINTS $ENV{MASSOCIATE_ROOT}/include
                $ENV{MASSOCIATE_INC_DIR}
                /usr/local/include)
find_library(MASSOCIATE_LIBRARY
             NAMES ${MASSOCIATE}
             PATHS $ENV{MASSOCIATE_ROOT}/lib/
                   $ENV{MASSOCIATE_LIB_DIR}
                   /usr/local/lib
                   /usr/local/lib64)

# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MASSOCIATE DEFAULT_MSG MASSOCIATE_LIBRARY MASSOCIATE_INCLUDE_DIR)
mark_as_advanced(MASSOCIATE_INCLUDE_DIR MASSOCIATE_LIBRARY)
