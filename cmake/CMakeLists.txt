set(PKG_DESC "SOCP")
set(PKG_EXTERNAL_DEPS "")

set(pkg_conf_file ${CMAKE_CURRENT_BINARY_DIR}/socp.pc)
configure_file(socp.pc.in ${pkg_conf_file} @ONLY)

install(FILES ${pkg_conf_file}
    DESTINATION ${SOCP_LIB_INSTALL_DIR}/pkgconfig/ COMPONENT pkgconfig)

install(FILES SOCPConfig.cmake DESTINATION share/cmake/SOCP)
