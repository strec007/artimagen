SUBDIRS=. cexample
AUTOMAKE_OPTIONS = foreign
ARTIMAGEN_LIBRARY_VERSION = 0:9:0

lib_LTLIBRARIES = libartimagen.la
libartimagen_la_SOURCES = libartimagen/image.cxx libartimagen/geometry.cxx libartimagen/sample.cxx \
   libartimagen/wrapper.cxx libartimagen/object.cxx
include_HEADERS=libartimagen/artimagen.h libartimagen/artimagen_i.h
libartimagen_la_LIBADD=$(TIFFLIB) $(FFTW3LIB) $(FFTW3TLIB) $(LUALIB)
libartimagen_la_LDFLAGS= -version-info $(ARTIMAGEN_LIBRARY_VERSION) -no-undefined

if HAVE_FFTW3_THREADS
libartimagen_la_LIBADD += $(FFTW3TLIB)
endif

if HAVE_LUA
bin_PROGRAMS=artimagenl
artimagenl_SOURCES=artimagenl.cxx
artimagenl_LDADD=libartimagen.la
endif

noinst_PROGRAMS=artimagen
artimagen_SOURCES=artimagen.cxx
artimagen_LDADD=libartimagen.la

CXXFLAGS=-O3 -ggdb -fstrict-aliasing -ffast-math -Wall
#Flags with the GNU profiler analysis enabled
#CXXFLAGS=-O3 -ggdb -pg -fstrict-aliasing -ffast-math -Wall
