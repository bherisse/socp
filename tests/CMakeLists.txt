file(GLOB Source_Files *.cpp)
source_group("Source Files" FILES ${Source_Files})

if (BUILD_GODDARD)
	add_executable(testGoddard testGoddard.cpp)
	target_link_libraries(testGoddard socp goddard)
	install (TARGETS testGoddard DESTINATION ${CMAKE_BINARY_DIR}/testsBIN)
endif ()