file(GLOB Source_Files *.cpp)
source_group("Source Files" FILES ${Source_Files})

if (BUILD_GODDARD)
	add_executable(testGoddard testGoddard.cpp)
	target_link_libraries(testGoddard socp goddard)
	install (TARGETS testGoddard DESTINATION bin)
endif ()

if (BUILD_VTOLUAV)
	add_executable(testVtolUAV testVtolUAV.cpp)
	target_link_libraries(testVtolUAV socp vtolUAV obstacle)
	install (TARGETS testVtolUAV DESTINATION bin)
endif ()
