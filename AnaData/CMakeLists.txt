set(INCLUDES
  ${ROOT_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/hipo/hipo4
  ${CLAS12AnaTools_INCLUDE_DIRS})
include_directories(${INCLUDES})
link_directories(${CMAKE_BINARY_DIR}/lib)
set(CMAKE_VERBOSE_MAKEFILE ON)

################################  AnaData.cc  ##########################################
add_executable(AnaData.exe AnaData.cc)
target_link_libraries(AnaData.exe Hipo4 ${ROOT_LIBRARIES} ${CLAS12AnaTools_LIBRARIES} )
#################################################################################################

################################  AnaDDVCS.cc ##########################################
add_executable(AnaDDVCS.exe AnaDDVCS.cc)
target_link_libraries(AnaDDVCS.exe Hipo4 ${ROOT_LIBRARIES} ${CLAS12AnaTools_LIBRARIES} )
#################################################################################################

################################  AnaDDVCSOSG.cc ##########################################
add_executable(AnaDDVCSOSG.exe AnaDDVCSOSG.cc)
target_link_libraries(AnaDDVCSOSG.exe Hipo4 ${ROOT_LIBRARIES} ${CLAS12AnaTools_LIBRARIES} )
#################################################################################################

################################ Installing Executable targets #############################
install( TARGETS AnaData.exe AnaDDVCS.exe AnaDDVCSOSG.exe
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaData/bin )

################################ Installing Non target files ###############################
install( FILES OSG_Validate.cc DrawPlots1.cc SkimOverAllRuns.py RunAllFilesOfRun.py  SkimOverAllRuns_v2.py
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaData/bin)

################################ Installing needed directores ##############################
install( DIRECTORY Skim_Data Skims_F18_Outbending Skims_F18_InbendingEarly Skims_F18_Inbending Figs OSG_ValidateFigs inpHipoFiles
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaData/ )
