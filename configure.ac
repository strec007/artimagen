#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.63)
define(ARTIMAGEN_VERSION,0.5a1)
AC_SUBST([ARTIMAGEN_VERSION])
AC_INIT(artimagen, ARTIMAGEN_VERSION, petr.cizmar@nist.gov)
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_DISABLE_SHARED
AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL

dnl initialize libtool
LT_INIT([win32-dll])
LT_LANG([C++])

# Compilation options
AC_ARG_WITH(lua, 
  AS_HELP_STRING([--without-lua],[disable LUA support (needed for saving)]),
 [],
 [with_lua=yes])

# Checks for libraries.
AC_CHECK_LIB([tiff], [TIFFOpen], [TIFFLIB=-ltiff],
AC_MSG_ERROR([libtiff check failed]))
AC_SUBST([TIFFLIB])

AC_CHECK_LIB([fftw3], [fftw_execute], [FFTW3LIB=-lfftw3],
AC_MSG_ERROR([fftw3 check failed]))
AC_SUBST([FFTW3LIB])

AC_CHECK_LIB([fftw3_threads], [fftw_init_threads], 
      [FFTW3TLIB="-lfftw3_threads -lm -lpthread"]
      AC_DEFINE([HAVE_FFTW3_THREADS], [1], [Define if you have fftw3_threads]),
AC_MSG_ERROR([fftw3 check failed]),[-lfftw3 -lpthread -lm])
AC_SUBST([FFTW3TLIB])

LUALIB=
AS_IF([test "x$with_lua" != xno],
    [AC_CHECK_LIB([lua5.1], [main],
     [AC_SUBST([LUALIB], ["-llua5.1"])
       AC_DEFINE([HAVE_LUA], [1], [Define if you have liblua5.1])
        [have_lua=yes]
       ],
     [AC_MSG_FAILURE([lua test failed (--without-lua to disable)])],
     [])
    ]
)

AM_CONDITIONAL(HAVE_LUA, test x$have_lua = xyes)


# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.

AC_HEADER_STDC
AC_CHECK_FUNCS([floor sqrt pow random srandom memset])
AC_CONFIG_FILES([Makefile src/Makefile src/cexample/Makefile])
AC_OUTPUT
