#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "NLopt::nlopt" for configuration "Release"
set_property(TARGET NLopt::nlopt APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(NLopt::nlopt PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;CXX"
  IMPORTED_LOCATION_RELEASE "/Users/wheelermw/Downloads/ToxicR/src/nlopt/lib/libnlopt.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS NLopt::nlopt )
list(APPEND _IMPORT_CHECK_FILES_FOR_NLopt::nlopt "/Users/wheelermw/Downloads/ToxicR/src/nlopt/lib/libnlopt.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
