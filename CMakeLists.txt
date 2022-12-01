cmake_minimum_required(VERSION 3.0.0)
project(pTFCE_C VERSION 0.1.0)

#include(CTest)
#enable_testing()
set(CMAKE_VERBOSE_MAKEFILE TRUE)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)


LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")
FIND_PACKAGE(FSL REQUIRED)
FIND_PACKAGE(GSL REQUIRED)

# For CMake < 3.9, we need to make the target ourselves
if(NOT TARGET OpenMP::OpenMP_CXX)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)

endif()
#target_link_libraries(MyTarget PUBLIC OpenMP::OpenMP_CXX)

MESSAGE(STATUS "FSL include dirs: ${FSL_INCLUDE_DIRS} ${FSLEXTRAS_INCLUDE_DIRS}")
MESSAGE(STATUS "FSL library dirs: ${FSL_LIBRARY_DIR} ${FSLEXTRAS_LIBRARY_DIR}")
MESSAGE(STATUS "GSL include dirs: ${GSL_INCLUDE_DIRS}")
MESSAGE(STATUS "GSL libraries: ${GSL_LIBRARIES}")
MESSAGE(STATUS "OpenMP flags: ${OpenMP_CXX_FLAGS}")

LINK_DIRECTORIES(${FSL_LIBRARY_DIR} ${FSLEXTRAS_LIBRARY_DIR})


# Create static library from sources
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -no-pie -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0")

SET(ptfce_SOURCES
    mathutil.cc
    pvalutil.cc
    grfClust.cc
    smoothest_ext.cc
    ptfce.cc
)

FILE(GLOB ptfce_HDR
    RELATIVE ${PROJECT_SOURCE_DIR}
    "*.h"
)

#ADD_LIBRARY( ptfce SHARED ${ptfce_SOURCES} ${ptfce_HDR})
#TARGET_INCLUDE_DIRECTORIES( ptfce PRIVATE ${PROJECT_SOURCE_DIR} ${FSL_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} )


# Create executables
ADD_EXECUTABLE( ptfce-bin ptfce_cmd.cc ${ptfce_SOURCES} ${ptfce_HDR})
TARGET_INCLUDE_DIRECTORIES( ptfce-bin PRIVATE ${PROJECT_SOURCE_DIR} ${FSL_INCLUDE_DIRS} )
TARGET_LINK_LIBRARIES( ptfce-bin utils warpfns basisfield newimage miscmaths prob NewNifti znz newmat z blas lapack gsl gslcblas )
SET_TARGET_PROPERTIES( ptfce-bin PROPERTIES
    OUTPUT_NAME ptfce
)

ADD_EXECUTABLE( ptfce-test pTFCE-test.cc ${ptfce_SOURCES} ${ptfce_HDR} )
TARGET_INCLUDE_DIRECTORIES( ptfce-test PRIVATE ${PROJECT_SOURCE_DIR} ${FSL_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} )
TARGET_LINK_LIBRARIES( ptfce-test utils warpfns basisfield newimage miscmaths prob NewNifti znz newmat z blas lapack gsl gslcblas )

ADD_EXECUTABLE( smoothest smoothest_app.cc smoothest_ext.cc mathutil.cc )
TARGET_INCLUDE_DIRECTORIES( smoothest PRIVATE ${PROJECT_SOURCE_DIR} ${FSL_INCLUDE_DIRS} )
TARGET_LINK_LIBRARIES( smoothest utils warpfns basisfield newimage miscmaths prob NewNifti znz newmat z blas lapack )
SET_TARGET_PROPERTIES( smoothest PROPERTIES
    OUTPUT_NAME smoothest-app
)

ADD_EXECUTABLE( connectedcompsize connectedcompsize.cc )
TARGET_INCLUDE_DIRECTORIES( connectedcompsize PRIVATE ${PROJECT_SOURCE_DIR} ${FSL_INCLUDE_DIRS} )
TARGET_LINK_LIBRARIES( connectedcompsize utils warpfns basisfield newimage miscmaths prob NewNifti znz newmat z blas lapack )
