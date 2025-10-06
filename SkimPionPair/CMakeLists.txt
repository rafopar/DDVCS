set(INCLUDES
  ${ROOT_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/hipo/hipo4
  ${CLAS12AnaTools_INCLUDE_DIRS})
include_directories(${INCLUDES})
link_directories(${CMAKE_BINARY_DIR}/lib)
set(CMAKE_VERBOSE_MAKEFILE ON)

################################  SkimPionPair.cc  ##########################################
add_executable(SkimPionPair.exe SkimPionPair.cc)
target_link_libraries(SkimPionPair.exe Hipo4 ${ROOT_LIBRARIES} ${CLAS12AnaTools_LIBRARIES})
#################################################################################################

################################  SkimPionPair.cc  ##########################################
add_executable(SkimSinglePion.exe SkimSinglePion.cc)
target_link_libraries(SkimSinglePion.exe Hipo4 ${ROOT_LIBRARIES} ${CLAS12AnaTools_LIBRARIES})
#################################################################################################


################################ Installing Executable targets #############################
install( TARGETS SkimPionPair.exe SkimSinglePion.exe
  DESTINATION ${CMAKE_INSTALL_PREFIX}/SkimPionPair/bin )

################################ Installing Non target files ###############################
install( FILES miscCodes.cc RunManySkims.py
  DESTINATION ${CMAKE_INSTALL_PREFIX}/SkimPionPair/bin)

################################ Installing needed directores ##############################
install( DIRECTORY Data Out_LUND Out_LUND_SinglePion
  DESTINATION ${CMAKE_INSTALL_PREFIX}/SkimPionPair/ )
