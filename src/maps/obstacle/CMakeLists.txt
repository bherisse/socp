file(GLOB Source_Files *.cpp)
source_group("Source Files" FILES ${Source_Files})

file(GLOB Header_Files *.hpp)
source_group("Header Files" FILES ${Header_Files})

add_library(obstacle ${Source_Files} ${Header_Files})

set_target_properties(obstacle PROPERTIES DEBUG_OUTPUT_NAME obstacle_d)
set_target_properties(obstacle PROPERTIES RELEASE_OUTPUT_NAME obstacle)

install (TARGETS obstacle
   		LIBRARY DESTINATION lib COMPONENT library
   		ARCHIVE DESTINATION lib COMPONENT library
   		RUNTIME DESTINATION bin COMPONENT library)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/maps/obstacle
				DESTINATION "include/maps" 
				FILES_MATCHING
				PATTERN "*.hpp"
				)
