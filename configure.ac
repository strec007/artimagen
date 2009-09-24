#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.63)
define(ARTIMAGEN_VERSION,0.4rc1)
AC_SUBST([ARTIMAGEN_VERSION])
AC_INIT(artimagen, ARTIMAGEN_VERSION, petr.cizmar@nist.gov)
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_DISABLE_SHARED
AC_PROG_LIBTOOL

# Checks for libraries.
AC_CHECK_LIB([tiff], [TIFFOpen], [TIFFLIB=-ltiff],
AC_MSG_ERROR([libtiff check failed]))
AC_SUBST([TIFFLIB])

AC_CHECK_LIB([fftw3], [fftw_execute], [FFTW3LIB=-lfftw3],
AC_MSG_ERROR([fftw3 check failed]))
AC_SUBST([FFTW3LIB])

AM_OPTIONS_WXCONFIG
reqwx=2.8.0
AM_PATH_WXCONFIG($reqwx, wxWin=1)
if test "$wxWin" != 1; then
#	AC_MSG_ERROR([
#		wxWidgets must be installed on your system.
# 
#		Please check that wx-config is in path, the directory
#		where wxWidgets libraries are installed (returned by
#		'wx-config --libs' or 'wx-config --static --libs' command)
#		is in LD_LIBRARY_PATH or equivalent variable and
#		wxWidgets version is $reqwx or above.
#		])
	AC_MSG_NOTICE([
		wxWidgets library is not installed properly.

		Not making the GUI wrapper (gaig).
		])
fi
if test "$wxWin" = 1; then
	wx_ok=true
fi
AM_CONDITIONAL([WX_OK], [test x$wx_ok = xtrue])


# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.

AC_HEADER_STDC
AC_CHECK_FUNCS([floor sqrt pow random srandom memset])
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
