# Already in cache, be silent
if (DAAL_INCLUDE_DIR AND DAAL_CORE_LIBRARY AND DAAL_SEQUENTIAL_LIBRARY)
   set (DAAL_FIND_QUIETLY TRUE)
endif()

#if (NOT BUILD_SHARED_LIBS)
#   set(IPPS "libipps.a")
#   set(VM "libippvm.a")
#   set(CORE "libippcore.a")
#else()
   set(CORE "onedal_core")
   #set(SEQUENTIAL "onedal_sequential")
   set(THREAD "onedal_thread")
#endif()

find_path(DAAL_INCLUDE_DIR
          NAMES daal.h
          HINTS $ENV{DAAL_ROOT}/include
                $ENV{DAAL_INCLUDE_DIR}
                /opt/intel/daal/include)
find_library(DAAL_CORE_LIBRARY
             NAMES ${CORE}
             PATHS $ENV{DAAL_ROOT}/lib/intel64
                   $ENV{DAAL_ROOT}/lib/
                   $ENV{DAAL_LIB_DIR}
                   /opt/intel/daal/lib/intel64
                   /opt/intel/daal/lib)
#find_library(DAAL_SEQUENTIAL_LIBRARY
#             NAMES ${SEQUENTIAL}
#             PATHS $ENV{DAAL_ROOT}/lib/intel64
#                   $ENV{DAAL_ROOT}/lib/
#                   $ENV{DAAL_LIB_DIR}
#                   /opt/intel/daal/lib/intel64
#                   /opt/intel/daal/lib)
find_library(DAAL_THREAD_LIBRARY
             NAMES ${THREAD}
             PATHS $ENV{DAAL_ROOT}/lib/intel64
                   $ENV{DAAL_ROOT}/lib/
                   $ENV{DAAL_LIB_DIR}
                   /opt/intel/daal/lib/intel64
                   /opt/intel/daal/lib)

#set(DAAL_LIBRARY ${DAAL_SEQUENTIAL_LIBRARY} ${DAAL_CORE_LIBRARY})
set(DAAL_LIBRARY ${DAAL_CORE_LIBRARY} ${DAAL_THREAD_LIBRARY})

# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FindDAAL DEFAULT_MSG DAAL_LIBRARY DAAL_INCLUDE_DIR)
mark_as_advanced(DAAL_INCLUDE_DIR DAAL_LIBRARY)
