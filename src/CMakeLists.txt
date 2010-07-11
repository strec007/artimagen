cmake_minimum_required (VERSION 2.8)

add_library (artimagen libartimagen/geometry.cxx libartimagen/image.cxx
libartimagen/object.cxx libartimagen/sample.cxx libartimagen/wrapper.cxx)
target_link_libraries (artimagen ${LIBS})

if (USE_LUA)
   add_executable (artimagenl artimagenl.cxx)
   target_link_libraries (artimagenl artimagen ${LIBS})
endif (USE_LUA)

add_executable (aigtest artimagen.cxx)
target_link_libraries (aigtest artimagen ${LIBS})
