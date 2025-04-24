cmake_minimum_required(VERSION 3.10)

# SOCP project
project(SOCP)
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)
set(SOCP_VERSION 0.1.0)
set(SOCP_INCLUDE_INSTALL_DIR "include/${PROJECT_NAME_LOWER}")
set(SOCP_LIB_INSTALL_DIR "lib")

if(POLICY CMP0078) 	# SWIG
  cmake_policy(SET CMP0078 NEW)
endif()
if(POLICY CMP0086)	# SWIG
  cmake_policy(SET CMP0086 NEW)
endif()
if(POLICY CMP0135)	# DOWNLOAD_EXTRACT_TIMESTAMP dans externalproject_add()
  cmake_policy(SET CMP0135  NEW)
endif()
if(POLICY CMP0148)	# Python
  cmake_policy(SET CMP0148 OLD)
endif()
if(POLICY CMP0167)	# Boost
  cmake_policy(SET CMP0167 OLD)
endif()

# Set the CMake module path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# type for real numbers (default is double)
#add_definitions(-DdoubleType)
#add_definitions(-DfloatType)

# Set the various project paths
set(BINARIES_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/binaries)
set(EXECUTABLE_OUTPUT_PATH ${BINARIES_OUTPUT_PATH}/${CMAKE_BUILD_TYPE})
set(LIBRARY_OUTPUT_PATH ${BINARIES_OUTPUT_PATH}/${CMAKE_BUILD_TYPE})

# To disable PkgConfig
set(CMAKE_DISABLE_FIND_PACKAGE_PkgConfig "ON")
macro( pkg_check_modules ) # a macro to replace pkg_check_modules (error otherwise since PkgConfig is disabled)
endmacro()

if (WIN32)
	add_definitions(-D_USE_MATH_DEFINES)
endif (WIN32)

if(UNIX)	
	if( NOT CMAKE_BUILD_TYPE )
		set( CMAKE_BUILD_TYPE Release CACHE STRING
			"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
			FORCE )
	endif()
	SET(GCC_11_COMPILE_FLAGS "-std=gnu++11")
	SET(ENABLE_MULTITHREADING_COMPILE_FLAGS "-pthread")
	#SET(GCC_FPERMISSIVE_COMPILE_FLAGS "-fpermissive")
	SET(GCC_POSITION_INDEPENDENT_CODE "-fPIC")
	SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_FPERMISSIVE_COMPILE_FLAGS} ${GCC_POSITION_INDEPENDENT_CODE} ${GCC_11_COMPILE_FLAGS} ${ENABLE_MULTITHREADING_COMPILE_FLAGS}" )
endif(UNIX)

# By default all wrappers are disabled
option(BUILD_PYTHON "Build Python Package" OFF)
if(BUILD_PYTHON)
	find_package(PythonInterp REQUIRED)
	find_package(PythonLibs REQUIRED)
	find_package(SWIG REQUIRED)
	if(${SWIG_FOUND})
		add_subdirectory(python)
	endif(${SWIG_FOUND})
endif()

# Add boost
option(Boost_USE "Use Boost" OFF)
if (${Boost_USE})
	add_definitions(-D_USE_BOOST)
	find_package(Boost REQUIRED) 
	if(Boost_FOUND)
		include_directories(${Boost_INCLUDE_DIRS}) 
	endif(Boost_FOUND)
endif()

# Add sources
add_subdirectory(src)

# Add tests
option(BUILD_TESTING "Build Tests" OFF)
if (BUILD_TESTING)
	add_subdirectory(tests)
endif()

# Add documentation
add_subdirectory(documentation)

# configuration file
add_subdirectory(cmake)