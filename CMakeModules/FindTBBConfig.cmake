# Already in cache, be silent
if (TBB_INCLUDE_DIR AND TBB_LIBRARY)
   set (TBB_FIND_QUIETLY TRUE)
endif()

set(TBB "tbb")

find_path(TBB_INCLUDE_DIR
          NAMES tbb
          HINTS /opt/intel/oneapi/tbb/latest/include
                /opt/intel/tbb/include
                $ENV{TBB_ROOT}/include
                $ENV{TBB_INC_DIR})
find_library(TBB_LIBRARY
             NAMES ${TBB}
             PATHS /opt/intel/oneapi/tbb/latest/lib/intel64/gcc4.8/
                   /opt/intel/tbb/lib/intel64
                   /opt/intel/tbb/lib
                   $ENV{TBB_ROOT}/lib/intel64
                   $ENV{TBB_ROOT}/lib/
                   $ENV{TBB_LIB_DIR})

#set(TBB_LIBRARY ${TBB_IPPS_LIBRARY} ${TBB_VM_LIBRARY} ${TBB_CORE_LIBRARY})

# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FindTBB DEFAULT_MSG TBB_LIBRARY TBB_INCLUDE_DIR)
mark_as_advanced(TBB_INCLUDE_DIR TBB_LIBRARY)
