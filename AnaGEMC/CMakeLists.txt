set(INCLUDES
  ${ROOT_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/hipo/hipo4
  ${CLAS12AnaTools_INCLUDE_DIRS})
include_directories(${INCLUDES})
link_directories(${CMAKE_BINARY_DIR}/lib)
set(CMAKE_VERBOSE_MAKEFILE ON)

################################  Compile Executables  ##########################################
add_executable(AnaGEMC.exe AnaGEMC.cc)
target_link_libraries(AnaGEMC.exe Hipo4 ${ROOT_LIBRARIES} ${CLAS12AnaTools_LIBRARIES} )

add_executable(UpdateConstants.exe UpdateConstants.cc)
target_link_libraries(UpdateConstants.exe Hipo4 ${ROOT_LIBRARIES} ROOT::Spectrum )

add_executable(DrawPlots.exe DrawPlots.cc)
target_link_libraries(DrawPlots.exe ${ROOT_LIBRARIES} )

add_executable(AnaPionRejection.exe  AnaPionRejection.cc)
target_link_libraries(AnaPionRejection.exe Hipo4 ${ROOT_LIBRARIES} ${CLAS12AnaTools_LIBRARIES})

add_executable(DrawPionRejection.exe  DrawPionRejection.cc)
target_link_libraries(DrawPionRejection.exe ${ROOT_LIBRARIES})

add_executable(DrawComparison.exe  DrawComparison.cc)
target_link_libraries(DrawComparison.exe ${ROOT_LIBRARIES})

#################################################################################################


################################ Installing Executable targets #############################
install( TARGETS AnaGEMC.exe UpdateConstants.exe DrawPlots.exe AnaPionRejection.exe DrawPionRejection.exe DrawComparison.exe
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaGEMC/bin )

################################ Installing Non target files ###############################
install( FILES DrawBgr_Estimation.cc DrawAngularResolutions.cc Draw_ISR_Effects.cc
  DrawInclusiveCoincidence.cc RunPionDataMany.py linkPars.py RunSettings.dat
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaGEMC/bin)

install( FILES RunSettings.dat
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaGEMC)


################################ Installing needed directores ##############################
install( DIRECTORY Data Latex Figs Pars Dumps
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaGEMC/ )
