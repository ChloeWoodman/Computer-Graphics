# Generated by CMake

if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" LESS 2.5)
   message(FATAL_ERROR "CMake >= 2.6.0 required")
endif()
cmake_policy(PUSH)
cmake_policy(VERSION 2.6...3.19)
#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Protect against multiple inclusion, which would fail when already imported targets are added once more.
set(_targetsDefined)
set(_targetsNotDefined)
set(_expectedTargets)
foreach(_expectedTarget glog)
  list(APPEND _expectedTargets ${_expectedTarget})
  if(NOT TARGET ${_expectedTarget})
    list(APPEND _targetsNotDefined ${_expectedTarget})
  endif()
  if(TARGET ${_expectedTarget})
    list(APPEND _targetsDefined ${_expectedTarget})
  endif()
endforeach()
if("${_targetsDefined}" STREQUAL "${_expectedTargets}")
  unset(_targetsDefined)
  unset(_targetsNotDefined)
  unset(_expectedTargets)
  set(CMAKE_IMPORT_FILE_VERSION)
  cmake_policy(POP)
  return()
endif()
if(NOT "${_targetsDefined}" STREQUAL "")
  message(FATAL_ERROR "Some (but not all) targets in this export set were already defined.\nTargets Defined: ${_targetsDefined}\nTargets not yet defined: ${_targetsNotDefined}\n")
endif()
unset(_targetsDefined)
unset(_targetsNotDefined)
unset(_expectedTargets)


# Create imported target glog
add_library(glog STATIC IMPORTED)

set_target_properties(glog PROPERTIES
  INTERFACE_COMPILE_DEFINITIONS "GLOG_NO_ABBREVIATED_SEVERITIES;GOOGLE_GLOG_DLL_DECL="
  INTERFACE_INCLUDE_DIRECTORIES "V:/GitHub/milestones-ChloeWoodman/Lab1/pbrt-v3/build/src/ext/glog;V:/GitHub/milestones-ChloeWoodman/Lab1/pbrt-v3/src/ext/glog/src;V:/GitHub/milestones-ChloeWoodman/Lab1/pbrt-v3/src/ext/glog/src/windows"
)

# Import target "glog" for configuration "Debug"
set_property(TARGET glog APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(glog PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "V:/GitHub/milestones-ChloeWoodman/Lab1/pbrt-v3/build/src/ext/glog/Debug/glog.lib"
  )

# Import target "glog" for configuration "Release"
set_property(TARGET glog APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(glog PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "V:/GitHub/milestones-ChloeWoodman/Lab1/pbrt-v3/build/src/ext/glog/Release/glog.lib"
  )

# Import target "glog" for configuration "MinSizeRel"
set_property(TARGET glog APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(glog PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "CXX"
  IMPORTED_LOCATION_MINSIZEREL "V:/GitHub/milestones-ChloeWoodman/Lab1/pbrt-v3/build/src/ext/glog/MinSizeRel/glog.lib"
  )

# Import target "glog" for configuration "RelWithDebInfo"
set_property(TARGET glog APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(glog PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELWITHDEBINFO "CXX"
  IMPORTED_LOCATION_RELWITHDEBINFO "V:/GitHub/milestones-ChloeWoodman/Lab1/pbrt-v3/build/src/ext/glog/RelWithDebInfo/glog.lib"
  )

# This file does not depend on other imported targets which have
# been exported from the same project but in a separate export set.

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
cmake_policy(POP)
