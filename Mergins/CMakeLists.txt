set(INCLUDES
  ${ROOT_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/hipo/hipo4
  ${CLAS12AnaTools_INCLUDE_DIRS})
include_directories(${INCLUDES})
link_directories(${CMAKE_BINARY_DIR}/lib)
set(CMAKE_VERBOSE_MAKEFILE ON)

################################  Compile Executables  ##########################################
add_executable(MergeRootFiles.exe MergeRootFiles.cc)
target_link_libraries(MergeRootFiles.exe ${ROOT_LIBRARIES} )
#################################################################################################

################################ Installing Executable targets #############################
install( TARGETS MergeRootFiles.exe 
  DESTINATION ${CMAKE_INSTALL_PREFIX}/Mergins/bin )

################################ Installing needed directores ##############################
install( DIRECTORY InpData OutData
  DESTINATION ${CMAKE_INSTALL_PREFIX}/Mergins/ )
