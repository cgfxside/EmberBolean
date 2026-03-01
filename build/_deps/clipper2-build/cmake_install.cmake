# Install script for directory: E:/Tools/EMBER/build/_deps/clipper2-src/CPP

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "E:/Tools/EMBER/install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "E:/Tools/EMBER/build/_deps/clipper2-build/Debug/Clipper2.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "E:/Tools/EMBER/build/_deps/clipper2-build/Release/Clipper2.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "E:/Tools/EMBER/build/_deps/clipper2-build/MinSizeRel/Clipper2.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "E:/Tools/EMBER/build/_deps/clipper2-build/RelWithDebInfo/Clipper2.lib")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/clipper2" TYPE FILE FILES
    "E:/Tools/EMBER/build/_deps/clipper2-src/CPP/Clipper2Lib/include/clipper2/clipper.h"
    "E:/Tools/EMBER/build/_deps/clipper2-src/CPP/Clipper2Lib/include/clipper2/clipper.version.h"
    "E:/Tools/EMBER/build/_deps/clipper2-src/CPP/Clipper2Lib/include/clipper2/clipper.core.h"
    "E:/Tools/EMBER/build/_deps/clipper2-src/CPP/Clipper2Lib/include/clipper2/clipper.engine.h"
    "E:/Tools/EMBER/build/_deps/clipper2-src/CPP/Clipper2Lib/include/clipper2/clipper.export.h"
    "E:/Tools/EMBER/build/_deps/clipper2-src/CPP/Clipper2Lib/include/clipper2/clipper.minkowski.h"
    "E:/Tools/EMBER/build/_deps/clipper2-src/CPP/Clipper2Lib/include/clipper2/clipper.offset.h"
    "E:/Tools/EMBER/build/_deps/clipper2-src/CPP/Clipper2Lib/include/clipper2/clipper.rectclip.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "E:/Tools/EMBER/build/_deps/clipper2-build/Clipper2.pc")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "E:/Tools/EMBER/build/_deps/clipper2-build/Debug/Clipper2.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "E:/Tools/EMBER/build/_deps/clipper2-build/Release/Clipper2.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "E:/Tools/EMBER/build/_deps/clipper2-build/MinSizeRel/Clipper2.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "E:/Tools/EMBER/build/_deps/clipper2-build/RelWithDebInfo/Clipper2.lib")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/clipper2" TYPE FILE FILES
    "E:/Tools/EMBER/build/_deps/clipper2-src/CPP/Clipper2Lib/include/clipper2/clipper.h"
    "E:/Tools/EMBER/build/_deps/clipper2-src/CPP/Clipper2Lib/include/clipper2/clipper.version.h"
    "E:/Tools/EMBER/build/_deps/clipper2-src/CPP/Clipper2Lib/include/clipper2/clipper.core.h"
    "E:/Tools/EMBER/build/_deps/clipper2-src/CPP/Clipper2Lib/include/clipper2/clipper.engine.h"
    "E:/Tools/EMBER/build/_deps/clipper2-src/CPP/Clipper2Lib/include/clipper2/clipper.export.h"
    "E:/Tools/EMBER/build/_deps/clipper2-src/CPP/Clipper2Lib/include/clipper2/clipper.minkowski.h"
    "E:/Tools/EMBER/build/_deps/clipper2-src/CPP/Clipper2Lib/include/clipper2/clipper.offset.h"
    "E:/Tools/EMBER/build/_deps/clipper2-src/CPP/Clipper2Lib/include/clipper2/clipper.rectclip.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/clipper2" TYPE FILE FILES
    "E:/Tools/EMBER/build/_deps/clipper2-build/Clipper2Config.cmake"
    "E:/Tools/EMBER/build/_deps/clipper2-build/Clipper2ConfigVersion.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/clipper2/Clipper2Targets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/clipper2/Clipper2Targets.cmake"
         "E:/Tools/EMBER/build/_deps/clipper2-build/CMakeFiles/Export/d8256b042a676486d1b2b19831a8c686/Clipper2Targets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/clipper2/Clipper2Targets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/clipper2/Clipper2Targets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/clipper2" TYPE FILE FILES "E:/Tools/EMBER/build/_deps/clipper2-build/CMakeFiles/Export/d8256b042a676486d1b2b19831a8c686/Clipper2Targets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/clipper2" TYPE FILE FILES "E:/Tools/EMBER/build/_deps/clipper2-build/CMakeFiles/Export/d8256b042a676486d1b2b19831a8c686/Clipper2Targets-debug.cmake")
  endif()
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/clipper2" TYPE FILE FILES "E:/Tools/EMBER/build/_deps/clipper2-build/CMakeFiles/Export/d8256b042a676486d1b2b19831a8c686/Clipper2Targets-minsizerel.cmake")
  endif()
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/clipper2" TYPE FILE FILES "E:/Tools/EMBER/build/_deps/clipper2-build/CMakeFiles/Export/d8256b042a676486d1b2b19831a8c686/Clipper2Targets-relwithdebinfo.cmake")
  endif()
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/clipper2" TYPE FILE FILES "E:/Tools/EMBER/build/_deps/clipper2-build/CMakeFiles/Export/d8256b042a676486d1b2b19831a8c686/Clipper2Targets-release.cmake")
  endif()
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "E:/Tools/EMBER/build/_deps/clipper2-build/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
