set(INCLUDES
  ${ROOT_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/hipo/hipo4
  ${CLAS12AnaTools_INCLUDE_DIRS})
include_directories(${INCLUDES})
link_directories(${CMAKE_BINARY_DIR}/lib)
set(CMAKE_VERBOSE_MAKEFILE ON)

################################  Compile Executables  ##########################################
add_executable( CreateTrainingSample CreateTrainingSample.cc )
set_target_properties(CreateTrainingSample PROPERTIES OUTPUT_NAME "CreateTrainingSample.exe")
target_link_libraries(CreateTrainingSample Hipo4 ${ROOT_LIBRARIES} ${CLAS12AnaTools_LIBRARIES} DDVCSToolslib)
#target_link_libraries(CreateTrainingSample Hipo4 ${ROOT_LIBRARIES})

#################################################################################################

################################ Installing Executable targets #############################
install(TARGETS CreateTrainingSample
DESTINATION ${CMAKE_INSTALL_PREFIX}/AIML_Train/el_DDVCS/bin )

install(FILES train_GEMC_elDDCVC.py Plot_Training_summary.py DrawNN.py
DESTINATION ${CMAKE_INSTALL_PREFIX}/AIML_Train/el_DDVCS/bin )

################################ Installing Directories #############################
install(DIRECTORY Figs Data TrainingSample pars
DESTINATION ${CMAKE_INSTALL_PREFIX}/AIML_Train/el_DDVCS )