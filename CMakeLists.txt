cmake_minimum_required (VERSION 2.8)
project (artimagen)

set (artimagen_VERSION_MAJOR 0)
set (artimagen_VERSION_MINOR 91)

configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )

include (CheckFunctionExists.cmake)
check_function_exists (rand HAVE_RAND)

add_subdirectory (src)
