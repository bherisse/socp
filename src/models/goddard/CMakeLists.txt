file(GLOB Source_Files *.cpp)
source_group("Source Files" FILES ${Source_Files})

file(GLOB Header_Files *.hpp)
source_group("Header Files" FILES ${Header_Files})

add_library(goddard ${Source_Files} ${Header_Files})
target_link_libraries(goddard socp)

set_target_properties(goddard PROPERTIES DEBUG_OUTPUT_NAME goddard_d)
set_target_properties(goddard PROPERTIES RELEASE_OUTPUT_NAME goddard)

install (TARGETS goddard
   		LIBRARY DESTINATION lib COMPONENT library
   		ARCHIVE DESTINATION lib COMPONENT library
   		RUNTIME DESTINATION bin COMPONENT library)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/models/goddard
				DESTINATION "include/models" 
				FILES_MATCHING
				PATTERN "*.hpp"
				PATTERN "*.pdf"
				PATTERN "*.bib"
				)
