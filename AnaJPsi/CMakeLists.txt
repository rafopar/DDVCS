set(INCLUDES
  ${ROOT_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/hipo/hipo4
  ${CLAS12AnaTools_INCLUDE_DIRS})
include_directories(${INCLUDES})
link_directories(${CMAKE_BINARY_DIR}/lib)
set(CMAKE_VERBOSE_MAKEFILE ON)

################################  Compile Executables  ##########################################
add_executable(AnaJPsi.exe AnaJPsi.cc)
target_link_libraries(AnaJPsi.exe Hipo4 ${ROOT_LIBRARIES} ${CLAS12AnaTools_LIBRARIES} )

add_executable(UpdateConstants_JPsi UpdateConstants.cc)
set_target_properties(UpdateConstants_JPsi PROPERTIES OUTPUT_NAME "UpdateConstants.exe")
target_link_libraries(UpdateConstants_JPsi ${ROOT_LIBRARIES} ROOT::Spectrum)

add_executable(DrawPlots_JPsi DrawPlots.cc)
set_target_properties(DrawPlots_JPsi PROPERTIES OUTPUT_NAME "DrawPlots.exe")
target_link_libraries(DrawPlots_JPsi ${ROOT_LIBRARIES} )

#################################################################################################

################################ Installing Executable targets #############################
install( TARGETS AnaJPsi.exe UpdateConstants_JPsi DrawPlots_JPsi
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaJPsi/bin )

################################ Installing needed directores ##############################
install( DIRECTORY Figs Pars Data
  DESTINATION ${CMAKE_INSTALL_PREFIX}/AnaJPsi/ )
