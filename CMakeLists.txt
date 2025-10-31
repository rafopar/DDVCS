cmake_minimum_required (VERSION VERSION 3.16)

# CMakeLists for tools
project(DDVCS_Studies VERSION 1.0.1 DESCRIPTION "DDVCS related Analysis codes and tools")

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_FLAGS "-std=c++17")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
set(CMAKE_PREFIX_PATH $ENV{ROOTSYS}/cmake)

find_package(ROOT REQUIRED)
find_package(LZ4 REQUIRED)
find_package(CLAS12AnaTools REQUIRED)
add_definitions(-D__LZ4__) # I might need to remove this completely (Rafo)
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)   # shared libs
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)   # static libs

include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})
set(DDVCS_SOURCE_FILES
  "${PROJECT_SOURCE_DIR}/src/DDVCSTools.cc" "${PROJECT_SOURCE_DIR}/src/DDVCSKine.cc" "${PROJECT_SOURCE_DIR}/src/ElectronDDVCSKine.cc")
set(DDVCS_HEADERS
  "${PROJECT_SOURCE_DIR}/include/DDVCSTools.h"  )
add_library(DDVCSToolslib SHARED ${DDVCS_SOURCE_FILES})
set_target_properties(DDVCSToolslib PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION 1)
set (CMAKE_INSTALL_LIBDIR "lib")

target_include_directories(DDVCSToolslib
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

include(GNUInstallDirs)

install(TARGETS DDVCSToolslib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_file(DDVCS_Studies.pc.in DDVCS_StudiesLib.pc @ONLY)

message("ROOT_INCLUDE_DIRS is " ${ROOT_INCLUDE_DIRS})

add_subdirectory(hipo/hipo4)
add_subdirectory(AnaData)
add_subdirectory(SkimPionPair)
add_subdirectory(AnaGEMC)
add_subdirectory(Mergins)
add_subdirectory(AnaJPsi)
add_subdirectory(AnaTwoPeg)
add_subdirectory(KinematicPlots)
add_subdirectory(AnaGrape)
add_subdirectory(RecoilDimensions)
add_subdirectory(DVCS_VGG_Results)
add_subdirectory(DDVCS_DVCS_xSecCpmpare)
add_subdirectory(AnaEPIC)
