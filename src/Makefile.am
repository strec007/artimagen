SUBDIRS=. cexample
AUTOMAKE_OPTIONS = foreign
ARTIMAGEN_LIBRARY_VERSION = 0:4:0

lib_LTLIBRARIES = libartimagen.la
libartimagen_la_SOURCES = libartimagen/image.cxx libartimagen/geometry.cxx libartimagen/sample.cxx \
   libartimagen/wrapper.cxx libartimagen/object.cxx
include_HEADERS=libartimagen/artimagen.h libartimagen/artimagen_i.h
libartimagen_la_LIBADD=$(TIFFLIB) $(FFTW3LIB) $(LUALIB)
libartimagen_la_LDFLAGS= -version-info $(ARTIMAGEN_LIBRARY_VERSION) -no-undefined


noinst_PROGRAMS=artimagen
artimagen_SOURCES=artimagen.cxx
artimagen_LDADD=libartimagen.la

#CXXFLAGS=-O3 -ggdb -fomit-frame-pointer -fstrict-aliasing -ffast-math -Wall
CXXFLAGS=-O3 -ggdb -fstrict-aliasing -ffast-math -Wall
