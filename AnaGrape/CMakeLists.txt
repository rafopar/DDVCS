set(INCLUDES
  ${ROOT_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/hipo/hipo4
  ${CLAS12AnaTools_INCLUDE_DIRS})
include_directories(${INCLUDES})
link_directories(${CMAKE_BINARY_DIR}/lib)
set(CMAKE_VERBOSE_MAKEFILE ON)

################################  Compile Executables  ##########################################
add_executable(AnaGrape AnaGrape.cc)
set_target_properties(AnaGrape PROPERTIES OUTPUT_NAME "AnaGrape.exe")
target_link_libraries(AnaGrape ${ROOT_LIBRARIES} )

add_executable(AnaGrape_22GeV AnaGrape_22GeV.cc)
set_target_properties(AnaGrape_22GeV PROPERTIES OUTPUT_NAME "AnaGrape_22GeV.exe")
target_link_libraries(AnaGrape_22GeV ${ROOT_LIBRARIES} )

#################################################################################################

################################ Installing Executable targets #############################
install( TARGETS AnaGrape AnaGrape_22GeV
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaGrape/bin )

################################ Installing Non target files ###############################
install( FILES DrawPlots.cc DrawPlots_22GeV.cc DrawAccEffects.cc
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaGrape/bin)

################################ Installing needed directores ##############################
install( DIRECTORY Figs Dumps
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaGrape/ )
