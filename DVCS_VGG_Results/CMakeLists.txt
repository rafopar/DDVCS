################################ Installing Non target files ###############################
install( FILES DrawPlots1.cc Draw22GeVPlots.cc DrawTCSPlots.cc
  DESTINATION ${CMAKE_INSTALL_PREFIX}/DVCS_VGG_Results/bin)

################################ Installing needed directores ##############################
install( DIRECTORY Figs
  DESTINATION ${CMAKE_INSTALL_PREFIX}/DVCS_VGG_Results/bin/ )
