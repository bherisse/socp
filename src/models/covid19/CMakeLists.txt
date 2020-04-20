file(GLOB Source_Files *.cpp)
source_group("Source Files" FILES ${Source_Files})

file(GLOB Header_Files *.hpp)
source_group("Header Files" FILES ${Header_Files})

add_library(covid19 ${Source_Files} ${Header_Files})
target_link_libraries(covid19 socp)

set_target_properties(covid19 PROPERTIES DEBUG_OUTPUT_NAME covid19_d)
set_target_properties(covid19 PROPERTIES RELEASE_OUTPUT_NAME covid19)

install (TARGETS covid19
   		LIBRARY DESTINATION lib COMPONENT library
   		ARCHIVE DESTINATION lib COMPONENT library
   		RUNTIME DESTINATION bin COMPONENT library)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/models/covid19
				DESTINATION "include/models" 
				FILES_MATCHING
				PATTERN "*.hpp"
				PATTERN "*.pdf"
				PATTERN "*.bib"
				)
