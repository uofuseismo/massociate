# Already in cache, be silent
if (HDF5_INCLUDE_DIR AND HDF5_LIBRARY)
   set (HDF5_FIND_QUIETLY TRUE)
endif()

find_path(HDF5_INCLUDE_DIR
          NAMES hdf5.h
          HINTS $ENV{HDF5_ROOT}/include
                /usr/local/include
                /usr/include)
find_library(HDF5_LIBRARY
             NAMES hdf5
             PATHS $ENV{HDF5_ROOT}/lib
                   $ENV{HDF5_ROOT}/lib64
                   /usr/local/lib
                   /usr/local/lib64)


# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDF5 DEFAULT_MSG HDF5_LIBRARY HDF5_INCLUDE_DIR)
mark_as_advanced(HDF5_INCLUDE_DIR HDF5_LIBRARY)
