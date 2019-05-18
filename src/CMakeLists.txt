# add directories
add_subdirectory (socp)

option(BUILD_GODDARD "Build Goddard" OFF)
option(BUILD_VTOLUAV "Build VTOL UAV" OFF)
option(BUILD_INTERCEPTOR "Build Interceptor" OFF)
add_subdirectory (maps)
add_subdirectory (models)



