set(INCLUDES
  ${ROOT_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/hipo/hipo4
  ${CLAS12AnaTools_INCLUDE_DIRS}
  ${TORCH_INCLUDE_DIRS}
  ~/work/include)
include_directories(${INCLUDES})
link_directories(${CMAKE_BINARY_DIR}/lib)
set(CMAKE_VERBOSE_MAKEFILE ON)

################################  Executables  ##########################################
add_executable(SkimDDVCS Skim_DDVCS.cc)
set_target_properties(SkimDDVCS PROPERTIES OUTPUT_NAME "SkimDDVCS.exe")
target_link_libraries(SkimDDVCS Hipo4 ${ROOT_LIBRARIES} ${CLAS12AnaTools_LIBRARIES} DDVCSToolslib AnaDB::anadb)
#################################################################################################

################################ Installing Executable targets #############################
install( TARGETS SkimDDVCS
  DESTINATION ${CMAKE_INSTALL_PREFIX}/elDDVCS/bin )

################################ Installing Non target files ###############################
install( FILES
  DESTINATION ${CMAKE_INSTALL_PREFIX}/elDDVCS/bin)

################################ Installing needed directores ##############################
install( DIRECTORY Train SkimDDVCS OutHists
  DESTINATION ${CMAKE_INSTALL_PREFIX}/elDDVCS/ )
