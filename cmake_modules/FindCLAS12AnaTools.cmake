# - Try to find CLAS12AnaTools
# Once done, this will define
#
#  CLAS12AnaTools_FOUND - system has CLAS12AnaTools
#  CLAS12AnaTools_INCLUDE_DIRS - the CLAS12AnaTools include directories
#  CLAS12AnaTools_LIBRARIES - link these to use CLAS12AnaTools

# Look for the header file
find_path(CLAS12AnaTools_INCLUDE_DIR
  NAMES clas12AnaTools.h RecParticle.h
  PATHS /home/rafopar/work/git/clas12AnaTools/include
)

# Look for the library file
find_library(CLAS12AnaTools_LIBRARY
  NAMES CLAS12AnaTools Hipo4
  PATHS /local/work/builds/CLAS12AnaTools/lib/
)

# Check if we found everything
if(CLAS12AnaTools_INCLUDE_DIR AND CLAS12AnaTools_LIBRARY)
  set(CLAS12AnaTools_FOUND TRUE)
  message(STATUS "FOUND CLAS12AnaTools")
else()
  set(CLAS12AnaTools_FOUND FALSE)
endif()

# Provide the include directories and libraries
if(CLAS12AnaTools_FOUND)
  message(STATUS "Setting CLAS12AnaTools_INCLUDE_DIRS and  CLAS12AnaTools_LIBRARIES variables...")
  set(CLAS12AnaTools_INCLUDE_DIRS ${CLAS12AnaTools_INCLUDE_DIR})
  set(CLAS12AnaTools_LIBRARIES ${CLAS12AnaTools_LIBRARY})

  message(STATUS "... CLAS12AnaTools_INCLUDE_DIRS = ${CLAS12AnaTools_INCLUDE_DIRS}")
  message(STATUS "... CLAS12AnaTools_LIBRARIES = ${CLAS12AnaTools_LIBRARIES}")
endif()

# Mark the cache entries as advanced
mark_as_advanced(CLAS12AnaTools_INCLUDE_DIR CLAS12AnaTools_LIBRARY)
