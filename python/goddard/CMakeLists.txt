# Setup SWIG
include(UseSWIG)

# include path for swig foobar.i 'target'
include_directories(../src)
set_property(SOURCE Goddard.i PROPERTY CPLUSPLUS ON)
# outdir for .py file
set(CMAKE_SWIG_OUTDIR ${BINARIES_OUTPUT_PATH}/python/goddard)
# outdir for .so file
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${BINARIES_OUTPUT_PATH}/python/goddard)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
# Second, for multi-config builds (e.g. msvc)
foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${BINARIES_OUTPUT_PATH}/python/goddard)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${BINARIES_OUTPUT_PATH}/python/goddard)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

# this will generate pyGoddard.py AND _%module%.so
swig_add_library(pyGoddard
	TYPE SHARED
	LANGUAGE python
	OUTPUT_DIR ${CMAKE_SWIG_OUTDIR}
	OUTFILE_DIR ${CMAKE_CURRENT_BINARY_DIR}
	SOURCES	Goddard.i)

# ${SWIG_MODULE_py${PROJECT_NAME}_REAL_NAME} will resolve to _pyGoddard.so TARGET
target_include_directories(${SWIG_MODULE_pyGoddard_REAL_NAME} PRIVATE
	${PYTHON_INCLUDE_DIRS})
target_link_libraries(${SWIG_MODULE_pyGoddard_REAL_NAME} PRIVATE
	goddard socp ${PYTHON_LIBRARIES})
if(UNIX)
	if(NOT APPLE)
		set_target_properties(${SWIG_MODULE_pyGoddard_REAL_NAME} PROPERTIES
			INSTALL_RPATH	"$ORIGIN")
	else()
		set_target_properties(${SWIG_MODULE_pyGoddard_REAL_NAME} PROPERTIES
			SUFFIX ".so"
			INSTALL_RPATH	"@loader_path")
		set_property(TARGET ${SWIG_MODULE_pyGoddard_REAL_NAME} APPEND PROPERTY
			LINK_FLAGS "-flat_namespace -undefined suppress")
	endif()
endif()

# Install directories
install(TARGETS ${SWIG_MODULE_pyGoddard_REAL_NAME}
		RUNTIME DESTINATION bin
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib)
install(FILES ${CMAKE_SWIG_OUTDIR}/models/pyGoddard.py DESTINATION bin)
