# Already in cache, be silent
if (MKL_INCLUDE_DIR AND MKL_LIBRARY)
    set(MKL_FIND_QUIETLY TRUE)
endif()

#if (NOT BUILD_SHARED_LIBS)
    #set(INTEL "libmkl_intel_lp64.a")
    #set(SEQUENTIAL "libmkl_sequential.a")
    #set(CORE "libmkl_core.a")
#else()
    set(INTEL mkl_intel_lp64)
    set(SEQUENTIAL mkl_sequential)
    set(CORE mkl_core)
#endif

# Find the include directory
find_path(MKL_INCLUDE_DIR
          NAMES mkl.h
          HINTS /opt/intel/mkl/include
                /opt/intel/oneapi/mkl/latest/include/
                $ENV{MKL_ROOT}/include)
# Find the library components
find_library(MKL_INTEL_LIBRARY
             NAMES ${INTEL}
             PATHS /opt/intel/oneapi/mkl/latest/lib/intel64
                   /opt/intel/mkl/lib/intel64
                   /opt/intel/mkl/lib 
                   $ENV{MKL_ROOT}/lib/intel64
                   $ENV{MKL_ROOT}/lib
                   $ENV{MKL_LIB_DIR}
            )
find_library(MKL_SEQUENTIAL_LIBRARY
             NAMES ${SEQUENTIAL}
             PATHS /opt/intel/oneapi/mkl/latest/lib/intel64
                   /opt/intel/mkl/lib/intel64
                   /opt/intel/mkl/lib
                   $ENV{MKL_ROOT}/lib/intel64
                   $ENV{MKL_ROOT}/lib
                   $ENV{MKL_LIB_DIR}
            )
find_library(MKL_CORE_LIBRARY
             NAMES ${CORE}
             PATHS /opt/intel/oneapi/mkl/latest/lib/intel64
                   /opt/intel/mkl/lib/intel64
                   /opt/intel/mkl/lib
                   $ENV{MKL_ROOT}/lib/intel64
                   $ENV{MKL_ROOT}/lib
                   $ENV{MKL_LIB_DIR}
            )

set(MKL_LIBRARY ${MKL_INTEL_LIBRARY} ${MKL_SEQUENTIAL_LIBRARY} ${MKL_CORE_LIBRARY})

# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FindMKL DEFAULT_MSG MKL_INCLUDE_DIR MKL_LIBRARY MKL_INTEL_LIBRARY MKL_SEQUENTIAL_LIBRARY MKL_CORE_LIBRARY)
mark_as_advanced(MKL_INCLUDE_DIR MKL_LIBRARY)
