###############################################################################
# Find SOCP
#
# This sets the following variables:
# SOCP_FOUND - True if SOCP was found.
# SOCP_INCLUDE_DIRS - Directories containing the SOCP include files.
# SOCP_LIBRARIES - Libraries needed to use SOCP.
# SOCP_LIBRARIES_DEBUG - Libraries needed to use SOCP, debug version.
# SOCP_DEFINITIONS - Compiler flags for SOCP.

find_package(PkgConfig)
pkg_check_modules(PC_SOCP socp)
set(SOCP_DEFINITIONS ${PC_SOCP_CFLAGS_OTHER})

find_path(SOCP_INCLUDE_DIR socp/shooting.hpp
          HINTS ${PC_SOCP_INCLUDEDIR} ${PC_SOCP_INCLUDE_DIRS}
		  PATHS "$ENV{PROGRAMFILES}/SOCP" "$ENV{PROGRAMW6432}/SOCP"
          PATH_SUFFIXES include)


			 
find_library(SOCP_LIBRARY_DEBUG NAMES socp_d
             HINTS ${PC_SOCP_LIBDIR} ${PC_SOCP_LIBRARY_DIRS} 
			 PATHS "$ENV{PROGRAMFILES}/SOCP" "$ENV{PROGRAMW6432}/SOCP"
			 PATH_SUFFIXES lib)

# Prefer static libraries in Windows over shared ones
if(WIN32)
	find_library(SOCP_LIBRARY NAMES socp
             HINTS ${PC_SOCP_LIBDIR} ${PC_SOCP_LIBRARY_DIRS} 
			 PATHS "$ENV{PROGRAMFILES}/SOCP" "$ENV{PROGRAMW6432}/SOCP"
			 PATH_SUFFIXES lib)
			 
	find_library(SOCP_LIBRARY_DEBUG NAMES socp_d
             HINTS ${PC_SOCP_LIBDIR} ${PC_SOCP_LIBRARY_DIRS} 
			 PATHS "$ENV{PROGRAMFILES}/SOCP" "$ENV{PROGRAMW6432}/SOCP"
			 PATH_SUFFIXES lib)
			 
	find_library(CMINPACK_LIBRARY 
			 NAMES cminpack_s cminpack
             HINTS ${PC_SOCP_LIBDIR} ${PC_SOCP_LIBRARY_DIRS} 
			 PATHS "$ENV{PROGRAMFILES}/SOCP" "$ENV{PROGRAMW6432}/SOCP"
			 PATH_SUFFIXES lib)
			 
	find_library(CMINPACK_LIBRARY_DEBUG 
               NAMES cminpack_s_d cminpack_d
               HINTS ${PC_SOCP_LIBDIR} ${PC_SOCP_LIBRARY_DIRS} 
               PATHS "$ENV{PROGRAMFILES}/SOCP" "$ENV{PROGRAMW6432}/SOCP"
               PATH_SUFFIXES lib)
else(WIN32)
	find_library(SOCP_LIBRARY NAMES socp
             HINTS ${PC_SOCP_LIBDIR} ${PC_SOCP_LIBRARY_DIRS} 
			 PATHS "$ENV{PROGRAMFILES}/SOCP" "$ENV{PROGRAMW6432}/SOCP"
			 PATH_SUFFIXES lib lib64)
			 
	find_library(SOCP_LIBRARY_DEBUG NAMES socp_d
             HINTS ${PC_SOCP_LIBDIR} ${PC_SOCP_LIBRARY_DIRS} 
			 PATHS "$ENV{PROGRAMFILES}/SOCP" "$ENV{PROGRAMW6432}/SOCP"
			 PATH_SUFFIXES lib lib64)

	find_library(CMINPACK_LIBRARY 
               NAMES cminpack_s cminpack
               HINTS ${PC_SOCP_LIBDIR} ${PC_SOCP_LIBRARY_DIRS}
               PATH_SUFFIXES lib lib64)

	find_library(CMINPACK_LIBRARY_DEBUG 
               NAMES cminpack_s_d cminpack_d
               HINTS ${PC_SOCP_LIBDIR} ${PC_SOCP_LIBRARY_DIRS}
               PATH_SUFFIXES lib lib64)
endif(WIN32)			 


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(socp  DEFAULT_MSG
                                  SOCP_LIBRARY SOCP_INCLUDE_DIR)

mark_as_advanced(SOCP_INCLUDE_DIR SOCP_LIBRARY)

set(SOCP_LIBRARIES optimized ${SOCP_LIBRARY} ${CMINPACK_LIBRARY} debug ${SOCP_LIBRARY_DEBUG} ${CMINPACK_LIBRARY_DEBUG})
set(SOCP_INCLUDE_DIRS ${SOCP_INCLUDE_DIR})