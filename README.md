# Migration-based Assocation

Herein are some tools for associating pick arrivals whose phase is known.  The method effectively migrates differential travel times.

# Build

This section outlines the build instructions though your mileage may vary depending on particulars of your system.
 
## Dependencies

The following dependencies must be satisified prior to building the code.

    1.  A C++17 compliant compiler.
    2.  [CMake](https://cmake.org/) version 3.10 or greater.
    3.  [GTest](https://github.com/google/googletest) for unit testing.
    4.  [Boost](https://www.boost.org/).
    5.  [DAAL](https://software.seek.intel.com/performance-libraries) for the DBSCAN clustering algorithm.
    6.  [MKL](https://software.seek.intel.com/performance-libraries) for a LAPACK/BLAS implementation.

Optional
    
    1.  [TBB](https://software.seek.intel.com/performance-libraries) for a portable threading implementation.

The build instructions are as follows:

   1.  Download the code:

    git clone https://github.com/uofuseismo/massociate.git


   2.  Descend into massociate

    cd massociate

   3.  Configure the software.  Typically I'll use a script such as

    #!/bin/bash
    export CXX=/usr/bin/clang++
    export MKL_ROOT=/opt/intel/mkl
    export DAAL_ROOT=/opt/intel/parallel_studio_xe_2020.0.015/compilers_and_libraries_2020/mac/daal/
    export BUILD_DIR=build
    export BOOST_ROOT=/usr/local/boost-1.70.0
    if [ -d ${BUILD_DIR} ]; then
       rm -rf ${BUILD_DIR}
    fi
    mkdir ${BUILD_DIR}
    cd ${BUILD_DIR}
    cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-g -Wall" \
    -DDAAL_INCLUDE_DIR=${DAAL_ROOT}/include \
    -DDAAL_CORE_LIBRARY=${DAAL_ROOT}/lib/libdaal_core.dylib \
    -DDAAL_SEQUENTIAL_LIBRARY=${DAAL_ROOT}/lib/libdaal_sequential.dylib

  4.  Following a successful configuration just type 

    make

  

