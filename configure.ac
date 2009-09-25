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
