#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "TBB::tbb" for configuration "RelWithDebInfo"
set_property(TARGET TBB::tbb APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(TBB::tbb PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELWITHDEBINFO "CXX;RC"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/tbb12.lib"
  )

list(APPEND _cmake_import_check_targets TBB::tbb )
list(APPEND _cmake_import_check_files_for_TBB::tbb "${_IMPORT_PREFIX}/lib/tbb12.lib" )

# Import target "TBB::tbbmalloc" for configuration "RelWithDebInfo"
set_property(TARGET TBB::tbbmalloc APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(TBB::tbbmalloc PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELWITHDEBINFO "CXX;RC"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/tbbmalloc.lib"
  )

list(APPEND _cmake_import_check_targets TBB::tbbmalloc )
list(APPEND _cmake_import_check_files_for_TBB::tbbmalloc "${_IMPORT_PREFIX}/lib/tbbmalloc.lib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
