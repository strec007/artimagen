AUTOMAKE_OPTIONS = foreign
ARTIMAGEN_LIBRARY_VERSION = 0:3:0

lib_LTLIBRARIES = libartimagen.la
libartimagen_la_SOURCES = libartimagen/image.cxx libartimagen/geometry.cxx libartimagen/sample.cxx
include_HEADERS=libartimagen/artimagen.h
libartimagen_la_LIBADD=$(TIFFLIB) $(FFTW3LIB)
libartimagen_la_LDFLAGS= -version-info $(ARTIMAGEN_LIBRARY_VERSION)


noinst_PROGRAMS=artimagen
artimagen_SOURCES=artimagen.cxx
artimagen_LDADD=libartimagen.la

bin_PROGRAMS=
gaig_SOURCES=gaig.cxx
gaig_CXXFLAGS=$(WX_CXXFLAGS)
gaig_LDADD=$(WX_LIBS) libartimagen.la
if WX_OK
bin_PROGRAMS += gaig
endif

CXXFLAGS=-O3 -ggdb -fomit-frame-pointer -malign-double -fstrict-aliasing -ffast-math -Wall
