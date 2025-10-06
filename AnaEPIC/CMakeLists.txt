set(INCLUDES
  ${ROOT_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/hipo/hipo4
  ${CLAS12AnaTools_INCLUDE_DIRS})
include_directories(${INCLUDES})
link_directories(${CMAKE_BINARY_DIR}/lib)
set(CMAKE_VERBOSE_MAKEFILE ON)

################################  Compile Executables  ##########################################
add_executable(AnaEPIC AnaEPIC.cc)
set_target_properties(AnaEPIC PROPERTIES OUTPUT_NAME "AnaEPIC.exe")
target_link_libraries(AnaEPIC ${ROOT_LIBRARIES})

################################ Installing Executable targets #############################
install( TARGETS AnaEPIC
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaEPIC/bin )

################################ Installing Non-Target files ##############################
install( FILES CompareGrrapeEpic.cc NotWorkingCode.cc
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaEPIC/ )


################################ Installing needed directores ##############################
install( DIRECTORY Figs
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaEPIC/ )
