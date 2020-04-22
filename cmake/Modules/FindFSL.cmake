#[=======================================================================[.rst:
FindFSL
--------

Find the native FMRIB's Software Library (FSL) includes and libraries.

The FMRIB's Software Library (FSL) is a comprehensive library of tools
for FMRI, MRI and DTI brain imaging data.

Imported Targets
^^^^^^^^^^^^^^^^

If FSL is found, this module defines the following :prop_tgt:`IMPORTED`
targets::

 FSL::fsl       - The main FSL library.
 FSL::fslextras - The FSL wrapper for third party libraries.

Result Variables
^^^^^^^^^^^^^^^^

This module will set the following variables in your project::

 FSL_FOUND                - True if FSL found on the local system
 FSL_INCLUDE_DIRS         - Location of FSL header files.
 FSL_LIBRARY_DIR_DIRS     - Location of FSL library files.
 FSL_VERSION              - The version of the discovered FSL install.

Hints
^^^^^

Set ``FSLDIR`` to a directory that contains a FSL installation.

This script expects to find libraries at ``$FSLDIR/lib`` and the FSL
headers at ``$FSLDIR/include``, and for fslextras at ``$FSLDIR/extras/lib`` and the
headers at ``$FSLDIR/extras/include``.

Cache Variables
^^^^^^^^^^^^^^^

This module may set the following variables depending on platform and type
of FSL installation discovered.  These variables may optionally be set to
help this module find the correct files::

 FSL_LIBRARY_DIR             - Location of the FSL library.

#]=======================================================================]

#include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
include(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

#=============================================================================
# If the user has provided ``FSLDIR``, use it!  Choose items found
# at this location over system locations.
if( EXISTS "$ENV{FSLDIR}" )
  file( TO_CMAKE_PATH "$ENV{FSLDIR}" FSLDIR )
  set( FSLDIR "${FSLDIR}" CACHE PATH "Prefix for FSL installation." )
endif()
if( NOT EXISTS "${FSLDIR}" )
  set( FSL_USE_PKGCONFIG ON )
endif()

#=============================================================================
# As a first try, use the PkgConfig module.  This will work on many
# *NIX systems.  See :module:`findpkgconfig`
# This will return ``FSL_INCLUDEDIR`` and ``FSL_LIBDIR`` used below.
if( FSL_USE_PKGCONFIG )
  find_package(PkgConfig)
  pkg_check_modules( FSL QUIET fsl )

  if( EXISTS "${FSL_INCLUDEDIR}" )
    get_filename_component( FSLDIR "${FSL_INCLUDEDIR}" DIRECTORY CACHE)
  endif()
endif()

#=============================================================================
# Set FSL_INCLUDE_DIRS and FSL_LIBRARIES. If we skipped the PkgConfig step, try
# to find the libraries at $FSLDIR (if provided) or in standard system
# locations.  These find_library and find_path calls will prefer custom
# locations over standard locations (HINTS).  If the requested file is not found
# at the HINTS location, standard system locations will be still be searched
# (/usr/lib64 (Redhat), lib/i386-linux-gnu (Debian)).

set( FSL_INCLUDE_DIR ${FSLDIR}/include )
set( FSL_LIBRARY_DIR ${FSLDIR}/lib )
set( FSLEXTRAS_INCLUDE_DIR ${FSLDIR}/extras/include )
set( FSLEXTRAS_LIBRARY_DIR ${FSLDIR}/extras/lib )
set( ARMAWRAP_INCLUDE_DIR ${FSLDIR}/extras/include/armawrap/armawrap )

set( FSL_INCLUDE_DIRS ${FSL_INCLUDE_DIR} ${FSLEXTRAS_INCLUDE_DIR} ${ARMAWRAP_INCLUDE_DIR} )
#set( FSL_LIBRARIES ${FSL_LIBRARY_DIR} ${FSLEXTRAS_LIBRARY_DIR} )

# If we didn't use PkgConfig, try to find the version via fsl-config or by
# reading fsl_version.h.
if( NOT FSL_VERSION )
  if( EXISTS "${FSLDIR}/etc/fslversion" )
    file( STRINGS "${FSLDIR}/etc/fslversion" FSL_VERSION )
  endif()
endif()

#=============================================================================
# handle the QUIETLY and REQUIRED arguments and set FSL_FOUND to TRUE if all
# listed variables are TRUE
find_package_handle_standard_args( FSL
  FOUND_VAR
    FSL_FOUND
  REQUIRED_VARS
    FSL_INCLUDE_DIR
    FSL_LIBRARY_DIR
    FSLEXTRAS_INCLUDE_DIR
    FSLEXTRAS_LIBRARY_DIR
  VERSION_VAR
    FSL_VERSION
    )

mark_as_advanced( FSLDIR FSL_VERSION FSL_LIBRARY_DIR FSL_INCLUDE_DIR
  FSLEXTRAS_LIBRARY_DIR FSLEXTRAS_INCLUDE_DIR
  FSL_USE_PKGCONFIG FSL_CONFIG )

#=============================================================================
# Register imported libraries:
# 1. If we can find a Windows .dll file, we will set appropriate target
#    properties for these.
# 2. However, for most systems, we will only register the import location and
#    include directory.

# Look for dlls.
if(WIN32)
  string( REPLACE ".lib" ".dll" FSL_LIBRARY_DIR_DLL       "${FSL_LIBRARY_DIR}" )
  string( REPLACE ".lib" ".dll" FSLEXTRAS_LIBRARY_DIR_DLL "${FSLEXTRAS_LIBRARY_DIR}" )
endif()

if( FSL_FOUND AND NOT TARGET FSL::fsl )
  if( EXISTS "${FSL_LIBRARY_DIR_DLL}" AND EXISTS "${FSLEXTRAS_LIBRARY_DIR_DLL}")

    # Windows systems with dll libraries.
    add_library( FSL::fsl       SHARED IMPORTED )
    add_library( FSL::fslextras SHARED IMPORTED )

    # Windows with dlls, but only Release libraries.
    set_target_properties( FSL::fslextras PROPERTIES
      IMPORTED_LOCATION_RELEASE         "${FSLEXTRAS_LIBRARY_DIR_DLL}"
      IMPORTED_IMPLIB                   "${FSLEXTRAS_LIBRARY_DIR}"
      INTERFACE_INCLUDE_DIRECTORIES     "${FSL_INCLUDE_DIRS}"
      IMPORTED_CONFIGURATIONS           Release
      IMPORTED_LINK_INTERFACE_LANGUAGES "CXX" )
    set_target_properties( FSL::fsl PROPERTIES
      IMPORTED_LOCATION_RELEASE         "${FSL_LIBRARY_DIR_DLL}"
      IMPORTED_IMPLIB                   "${FSL_LIBRARY_DIR}"
      INTERFACE_INCLUDE_DIRECTORIES     "${FSL_INCLUDE_DIRS}"
      IMPORTED_CONFIGURATIONS           Release
      IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
      INTERFACE_LINK_LIBRARIES          FSL::fslextras )

  else()

    # For all other environments (ones without dll libraries), create
    # the imported library targets.
    add_library( FSL::fsl       UNKNOWN IMPORTED )
    add_library( FSL::fslextras UNKNOWN IMPORTED )
    set_target_properties( FSL::fslextras PROPERTIES
      IMPORTED_LOCATION                 "${FSLEXTRAS_LIBRARY_DIR}"
      INTERFACE_INCLUDE_DIRECTORIES     "${FSL_INCLUDE_DIRS}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "CXX" )
    set_target_properties( FSL::fsl PROPERTIES
      IMPORTED_LOCATION                 "${FSL_LIBRARY_DIR}"
      INTERFACE_INCLUDE_DIRECTORIES     "${FSL_INCLUDE_DIRS}"
      IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
      INTERFACE_LINK_LIBRARIES          FSL::fslextras )
  endif()
endif()
