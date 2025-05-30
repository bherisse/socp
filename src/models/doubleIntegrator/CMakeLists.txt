# Crée des variables avec les fichiers à compiler
file(GLOB Source_Files *.cpp)
source_group("Source Files" FILES ${Source_Files})

file(GLOB Header_Files *.hpp)
source_group("Header Files" FILES ${Header_Files})

# On indique que l'on veut un exécutable "socp" compilé à partir des fichiers décrits par les variables SRCS et HEADERS
add_library(doubleIntegrator ${Source_Files} ${Header_Files})
target_link_libraries(doubleIntegrator socp)

set_target_properties(doubleIntegrator PROPERTIES DEBUG_OUTPUT_NAME doubleIntegrator_d)
set_target_properties(doubleIntegrator PROPERTIES RELEASE_OUTPUT_NAME doubleIntegrator)

install (TARGETS doubleIntegrator
   		LIBRARY DESTINATION lib COMPONENT library
   		ARCHIVE DESTINATION lib COMPONENT library
   		RUNTIME DESTINATION bin COMPONENT library)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/models/doubleIntegrator
				DESTINATION "include/models" 
				FILES_MATCHING
				PATTERN "*.hpp"
				PATTERN "*.pdf"
				PATTERN "*.bib"
				)

