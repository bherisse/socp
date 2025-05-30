file(GLOB Source_Files *.cpp)
source_group("Source Files" FILES ${Source_Files})

file(GLOB Header_Files *.hpp)
source_group("Header Files" FILES ${Header_Files})

add_library(socp ${Source_Files} ${Header_Files})

# Add the external dependencies
include(ExternalProject)
externalproject_add(CMinPack
					GIT_REPOSITORY https://github.com/devernay/cminpack.git
					#URL https://github.com/bherisse/cminpack/archive/patch-1.zip
					#DOWNLOAD_EXTRACT_TIMESTAMP ON
					CMAKE_ARGS	-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_SOURCE_DIR}/extDep/CMinPack
							#-DCMAKE_DEBUG_POSTFIX:STRING=_d
							#-DCMAKE_RELEASE_POSTFIX:STRING=
							-DCMAKE_POLICY_VERSION_MINIMUM=3.10
							-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
							-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
							#-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
							-DBUILD_EXAMPLES:BOOL=OFF
							-DUSE_BLAS:BOOL=OFF
					)

add_dependencies(socp CMinPack)

set(CMinpack_DIR ${CMAKE_SOURCE_DIR}/extDep/CMinPack)
set(CMAKE_PREFIX_PATH ${CMAKE_SOURCE_DIR}/extDep/CMinPack)
if(EXISTS ${CMinpack_DIR})
	find_package(CMinpack REQUIRED)
	get_target_property(CMINPACK_INCLUDE_DIRS cminpack::cminpack INTERFACE_INCLUDE_DIRECTORIES)
	include_directories(${CMINPACK_INCLUDE_DIRS})
	target_link_libraries(socp cminpack::cminpack)
endif(EXISTS ${CMinpack_DIR})

set_target_properties(socp PROPERTIES DEBUG_OUTPUT_NAME socp_d)
set_target_properties(socp PROPERTIES RELEASE_OUTPUT_NAME socp)

install (TARGETS socp
   		LIBRARY DESTINATION lib COMPONENT library
   		ARCHIVE DESTINATION lib COMPONENT library
   		RUNTIME DESTINATION bin COMPONENT library)
install (FILES ${CMINPACK_LIBRARY} DESTINATION "lib")
install (FILES ${CMINPACK_LIBRARY_DEBUG} DESTINATION "lib")
install (DIRECTORY ${CMAKE_SOURCE_DIR}/src/socp DESTINATION "include" FILES_MATCHING PATTERN "*.hpp")
