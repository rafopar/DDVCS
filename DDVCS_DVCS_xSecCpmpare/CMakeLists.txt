################################ Installing Non target files ###############################
install( FILES DrawDVCS_DDVCS_Compare.cc
  DESTINATION ${CMAKE_INSTALL_PREFIX}/DDVCS_DVCS_xSecCpmpare/bin)

################################ Installing needed directores ##############################
install( DIRECTORY Figs
  DESTINATION ${CMAKE_INSTALL_PREFIX}/DDVCS_DVCS_xSecCpmpare )
