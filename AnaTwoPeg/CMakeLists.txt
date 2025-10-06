set(INCLUDES
  ${ROOT_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/hipo/hipo4
  ${CLAS12AnaTools_INCLUDE_DIRS})
include_directories(${INCLUDES})
link_directories(${CMAKE_BINARY_DIR}/lib)
set(CMAKE_VERBOSE_MAKEFILE ON)

################################  Compile Executables  ##########################################
add_executable(AnaTwoPeg AnaTwoPeg.cc)
set_target_properties(AnaTwoPeg PROPERTIES OUTPUT_NAME "AnaTwoPeg.exe")
target_link_libraries(AnaTwoPeg ${ROOT_LIBRARIES} )
#################################################################################################

################################ Installing Executable targets #############################
install( TARGETS AnaTwoPeg
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaTwoPeg/bin )

################################ Installing needed directores ##############################
install( DIRECTORY Data
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaTwoPeg/ )
