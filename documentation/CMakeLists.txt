
set(AUTO_DOXYGEN TRUE CACHE BOOL "If FALSE, don't generate Doxygen documentation unless you type make doxygen.")

find_package(Doxygen)

#if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/html)
#	file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/html)
#endif(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/html)

if(DOXYGEN_FOUND)
	if(AUTO_DOXYGEN)
		add_custom_target(	DOCUMENTATION
							ALL
							${DOXYGEN_EXECUTABLE} Doxyfile
							WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
	else(AUTO_DOXYGEN)
		add_custom_target(	DOCUMENTATION
							${DOXYGEN_EXECUTABLE} Doxyfile
							WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
	endif(AUTO_DOXYGEN)

	install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/SOCP DESTINATION documentation)
	install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/CMinpack DESTINATION documentation)

endif(DOXYGEN_FOUND)