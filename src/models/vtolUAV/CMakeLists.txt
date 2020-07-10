file(GLOB Source_Files *.cpp)
source_group("Source Files" FILES ${Source_Files})

file(GLOB Header_Files *.hpp)
source_group("Header Files" FILES ${Header_Files})

add_library(vtolUAV ${Source_Files} ${Header_Files})
target_link_libraries(vtolUAV socp)

set_target_properties(vtolUAV PROPERTIES DEBUG_OUTPUT_NAME vtolUAV_d)
set_target_properties(vtolUAV PROPERTIES RELEASE_OUTPUT_NAME vtolUAV)

install (TARGETS vtolUAV
   		LIBRARY DESTINATION lib COMPONENT library
   		ARCHIVE DESTINATION lib COMPONENT library
   		RUNTIME DESTINATION bin COMPONENT library)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/models/vtolUAV
				DESTINATION "include/models" 
				FILES_MATCHING
				PATTERN "*.hpp"
				PATTERN "*.pdf"
				PATTERN "*.bib"
				)
