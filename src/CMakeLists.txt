cmake_minimum_required (VERSION 2.8)

add_library (artimagen libartimagen/geometry.cxx libartimagen/image.cxx
libartimagen/object.cxx libartimagen/sample.cxx libartimagen/wrapper.cxx)
target_link_libraries (artimagen ${LIBS})
SET_TARGET_PROPERTIES(artimagen PROPERTIES VERSION 0.9.3 SOVERSION 0)


if (USE_LUA)
   add_executable (artimagenl artimagenl.cxx)
   target_link_libraries (artimagenl artimagen ${LIBS})
   install (TARGETS artimagenl DESTINATION bin)
endif (USE_LUA)

add_executable (aigtest artimagen.cxx)
target_link_libraries (aigtest artimagen ${LIBS})

install (TARGETS artimagen DESTINATION lib)
install (FILES "${PROJECT_SOURCE_DIR}/src/libartimagen/artimagen.h" DESTINATION include)
install (FILES "${PROJECT_SOURCE_DIR}/src/libartimagen/artimagen_i.h" DESTINATION include)
