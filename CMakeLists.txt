cmake_minimum_required (VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

project (artimagen)
set (CMAKE_BUILD_TYPE release)

set (artimagen_VERSION_MAJOR 0)
set (artimagen_VERSION_MINOR 91)

option(USE_LUA "Compile with Lua support" ON)
option(COMPILE_SHARED "Compile libartimagen as a shared lib." OFF)

if (COMPILE_SHARED)
   set (BUILD_SHARED_LIBS ON)
else (COMPILE_SHARED)
   set (BUILD_SHARED_LIBS OFF)
endif (COMPILE_SHARED)

set (PACKAGE_VERSION ${artimagen_VERSION_MAJOR}.${artimagen_VERSION_MINOR} )

include (CheckFunctionExists.cmake)
CHECK_FUNCTION_EXISTS (rand HAVE_RAND)


find_library(TIFF_LIB NAMES tiff)
if (NOT TIFF_LIB)
	message(FATAL_ERROR, "Cannot find libtiff library.")
endif (NOT TIFF_LIB)
set(LIBS ${LIBS} ${TIFF_LIB})

find_library(FFTW_LIB NAMES fftw3)
if (NOT FFTW_LIB)
	message(FATAL_ERROR, "Cannot find libfftw3.")
endif (NOT FFTW_LIB)
set(LIBS ${LIBS} ${FFTW_LIB})

if (USE_LUA)
   find_library(LUA_LIB NAMES lua5.1)
   if (NOT LUA_LIB)
	   message(FATAL_ERROR, "Cannot find liblua5.1. Switch off Lua support.")
   endif (NOT LUA_LIB)
   set(LIBS ${LIBS} ${LUA_LIB})
   set (HAVE_LUA 1)
endif (USE_LUA)

set (CMAKE_CXX_FLAGS_RELEASE "-O3 -fstrict-aliasing -ffast-math")
set (CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb  -fstrict-aliasing -ffast-math -Wall")

configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )

add_subdirectory (src)
