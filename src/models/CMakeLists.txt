##############################################
######             MODELS          ###########
##############################################

if (BUILD_GODDARD)
	add_subdirectory (goddard)
endif ()

if (BUILD_VTOLUAV)
	add_subdirectory (vtolUAV)
endif ()