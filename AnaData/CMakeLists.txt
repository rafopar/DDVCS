set(INCLUDES
  ${ROOT_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/hipo/hipo4
  ${CLAS12AnaTools_INCLUDE_DIRS})
include_directories(${INCLUDES})
link_directories(${CMAKE_BINARY_DIR}/lib)
set(CMAKE_VERBOSE_MAKEFILE ON)

################################  Executables  ##########################################
add_executable(AnaData.exe AnaData.cc)
target_link_libraries(AnaData.exe Hipo4 ${ROOT_LIBRARIES} ${CLAS12AnaTools_LIBRARIES} )

add_executable(AnaDDVCS.exe AnaDDVCS.cc)
target_link_libraries(AnaDDVCS.exe Hipo4 ${ROOT_LIBRARIES} ${CLAS12AnaTools_LIBRARIES} DDVCSToolslib )

add_executable(AnaDDVCSOSG.exe AnaDDVCSOSG.cc)
target_link_libraries(AnaDDVCSOSG.exe Hipo4 ${ROOT_LIBRARIES} ${CLAS12AnaTools_LIBRARIES} )

add_executable(TriggerStudy TriggerStudy.cc)
set_target_properties(TriggerStudy PROPERTIES OUTPUT_NAME "TriggerStudy.exe")
target_link_libraries(TriggerStudy Hipo4 ${ROOT_LIBRARIES} ${CLAS12AnaTools_LIBRARIES})

add_executable(AnaElectronDDVCS AnaElectronDDVCS.cc)
set_target_properties(AnaElectronDDVCS PROPERTIES OUTPUT_NAME "AnaElectronDDVCS")
target_link_libraries(AnaElectronDDVCS Hipo4 ${ROOT_LIBRARIES} ${CLAS12AnaTools_LIBRARIES} )
#################################################################################################

################################ Installing Executable targets #############################
install( TARGETS AnaData.exe AnaDDVCS.exe AnaDDVCSOSG.exe TriggerStudy AnaElectronDDVCS
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaData/bin )

################################ Installing Non target files ###############################
install( FILES OSG_Validate.cc DrawPlots1.cc SkimOverAllRuns.py RunAllFilesOfRun.py  SkimOverAllRuns_v2.py
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaData/bin)

install( FILES RunSettings.dat
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaData)

################################ Installing needed directores ##############################
install( DIRECTORY Skim_Data Skims_F18_Outbending Skims_F18_InbendingEarly Skims_F18_Inbending Figs OSG_ValidateFigs inpHipoFiles HistData
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaData/ )
