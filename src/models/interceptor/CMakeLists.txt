find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
	include_directories(${EIGEN3_INCLUDE_DIRS}) 
endif(Eigen3_FOUND)

# Crée des variables avec les fichiers à compiler
file(GLOB Source_Files *.cpp)
source_group("Source Files" FILES ${Source_Files})

file(GLOB Header_Files *.hpp)
source_group("Header Files" FILES ${Header_Files})

# On indique que l'on veut un exécutable "socp" compilé à partir des fichiers décrits par les variables SRCS et HEADERS
add_library(interceptor ${Source_Files} ${Header_Files})
target_link_libraries(interceptor socp)

set_target_properties(interceptor PROPERTIES DEBUG_OUTPUT_NAME interceptor_d)
set_target_properties(interceptor PROPERTIES RELEASE_OUTPUT_NAME interceptor)

install (TARGETS interceptor
   		LIBRARY DESTINATION lib COMPONENT library
   		ARCHIVE DESTINATION lib COMPONENT library
   		RUNTIME DESTINATION bin COMPONENT library)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/models/interceptor
				DESTINATION "include/models" 
				FILES_MATCHING
				PATTERN "*.hpp"
				PATTERN "*.pdf"
				PATTERN "*.bib"
				)
