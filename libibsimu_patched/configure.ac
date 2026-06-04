AC_INIT([libibsimu-1.0.6dev], [1.0.6dev], [taneli.kalvas@jyu.fi])
AC_CANONICAL_HOST
#AM_CONFIG_HEADER(src/config.h)
AC_CONFIG_HEADERS(src/config.h)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE
AC_LANG(C++)

dnl AC_USE_SYSTEM_EXTENSIONS

AC_DEFUN(
	AC_GCC_SUPPORTS_RDYNAMIC,
	[
		AC_PATH_PROG([grep], [grep])
                echo -n "checking if the compiler supports -rdynamic... "
                if test "X${ac_compiler_gnu}X" = "XyesX"; then
                        if test "X`${CC} -dumpspecs | ${grep} rdynamic`X" != "XX"; then
                                echo "yes"
                                AC_SUBST([rdynamic_flag], ["-rdynamic"])
                        else
                                echo "no"
                                AC_SUBST([rdynamic_flag], [])
                        fi
                else
                        echo "no"
                        AC_SUBST([rdynamic_flag], [])
                fi
        ]
)
AC_GCC_SUPPORTS_RDYNAMIC([rdynamic_flag])

lt_enable_auto_import=""
case "$host_os" in
   mingw* | cegcc*)
      lt_enable_auto_import="-Wl,--enable-auto-import"
esac
AC_SUBST(lt_enable_auto_import)

dnl ---------------------------------
dnl Check for programs
dnl ---------------------------------

dnl m4_pattern_allow([AC_PROG_LIBTOOL])
dnl AC_PROG_LIBTOOL
LT_INIT
AC_PROG_CXX
dnl AC_PROG_RANLIB

dnl ---------------------------------
dnl pkgconfig
dnl ---------------------------------

PKG_PROG_PKG_CONFIG()
if test "x$PKG_CONFIG" = "x" ; then
   AC_MSG_ERROR([pkg-config >= $PKGCONFIG_REQUIRED required but not
found (http://pkgconfig.freedesktop.org/)])
fi

dnl -------------
dnl zlib

CHECK_ZLIB()

dnl ---------------------------------
dnl Check for selectable compilation
dnl ---------------------------------

AC_ARG_ENABLE([sigsegv_stack],
  AS_HELP_STRING([--disable-sigsegv_stack], [Disable stack printing on SEGV]))

AS_IF([test "x$enable_sigsegv_stack" != "xno"], [
  CPPFLAGS="$CPPFLAGS -DSIGSEGV_STACK"
])

dnl Gtk
dnl Check by default, can be disabled
dnl ---------------------------------
AC_ARG_WITH([gtk3],
  [AS_HELP_STRING([--with-gtk3],
                  [enable interactive gtk plotter @<:@default=check@:>@])],
                  [],
                  [with_gtk3=check])
AS_IF([test "x$with_gtk3" != "xno"],
      [PKG_CHECK_MODULES([GTK3], [gtk+-3.0], [], [echo "Disabling functionality"])
       if $PKG_CONFIG --exists gtk+-3.0 ; then
         USEGTK3="yes"
         AC_DEFINE([GTK3], [1], [Define to 1 if Gtk3 dependent code is to be compiled.])
         LIBS="$LIBS $GTK3_LIBS"
         CPPFLAGS="$CPPFLAGS $GTK3_CFLAGS"
         AC_SUBST(GTK3_COND,[gtk+-3.0])
       elif test "x$with_gtk3" != "xcheck" ; then
         AC_MSG_ERROR([--with-gtk3 was given, but pkgconfig check for gtk+-3.0 library failed])
       else
         AC_SUBST(GTK3_COND,"")
       fi]
)
AM_CONDITIONAL([GTK3], [test x$USEGTK3 = xyes])

dnl OpenGL/GtkGLExt
dnl Check by default, can be disabled
dnl ---------------------------------
AC_ARG_WITH([opengl],
  [AS_HELP_STRING([--with-opengl],
                  [enable opengl 3d geometry plotter @<:@default=check@:>@])],
                  [],
                  [with_opengl=check])
AS_IF([test "x$with_opengl" != "xno"],
      [PKG_CHECK_MODULES([GTKGLEXT], [gtkglext-3.0], [], [echo "Using software rendering"])
       if $PKG_CONFIG --exists gtkglext-3.0 ; then
         USEOPENGL="yes"
         AC_DEFINE([OPENGL], [1], [Define to 1 if GtkGLExt/OpenGL dependent code is to be compiled.])
         LIBS="$LIBS $GTKGLEXT_LIBS"
         CPPFLAGS="$CPPFLAGS $GTKGLEXT_CFLAGS"
         AC_SUBST(OPENGL_COND,[gtkglext-3.0])
       elif test "x$with_opengl" != "xcheck" ; then
         AC_MSG_ERROR([--with-opengl was given, but pkgconfig check for gtkglext-3.0 library failed])
       else
         AC_SUBST(OPENGL_COND,"")
       fi]
)
AM_CONDITIONAL([OPENGL], [test x$USEOPENGL = xyes])


CHECK_UMFPACK()
AM_CONDITIONAL([UMFPACK], [test "x$USEUMFPACK" = "xyes"])


dnl CSG
dnl Check by default, can be disabled
dnl ---------------------------------
AC_ARG_WITH([csg],
  [AS_HELP_STRING([--with-csg],
                  [enable support for csg solid @<:@default=check@:>@])],
                  [],
                  [with_csg=check])
AS_IF([test "x$with_csg" != "xno"],
      [PKG_CHECK_MODULES([CSG], [csg], [], [echo "Disabling functionality"])
       if $PKG_CONFIG --exists csg ; then
         USECSG="yes"
         AC_DEFINE([CSG], [1], [Define to 1 if CSG solid is to be compiled.])
         LIBS="$LIBS $CSG_LIBS"
         CPPFLAGS="$CPPFLAGS $CSG_CFLAGS"
         AC_SUBST(CSG_COND,[csg])
       elif test "x$with_csg" != "xcheck" ; then
         AC_MSG_ERROR([--with-csg was given, but pkgconfig check for csg library failed])
       else
         AC_SUBST(CSG_COND,"")
       fi]
)
AM_CONDITIONAL([CSG], [test x$USECSG = xyes])

dnl ---------------------------------
dnl Check for libraries
dnl ---------------------------------

dnl -------------
dnl libm

AC_CHECK_LIBM
CPPFLAGS="$CFLAGS $CPPFLAGS"
LIBS="$LIBS $LIBM"

dnl -------------
dnl gtk+-3.0

dnl AM_PATH_GTK_3_0([3.6.0])

LIBS="$LIBS $GTK_LIBS"
CPPFLAGS="$CPPFLAGS $GTK_CFLAGS"


dnl -------------
dnl libpng

if $PKG_CONFIG --exists libpng ; then
  PKG_CHECK_MODULES([LIBPNG], [libpng])
else
  AC_MSG_ERROR([required library libpng not found])
fi
LIBS="$LIBS $LIBPNG_LIBS"
CPPFLAGS="$CPPFLAGS $LIBPNG_CFLAGS"


dnl -------------
dnl fontconfig

if $PKG_CONFIG --exists fontconfig ; then
  PKG_CHECK_MODULES([FONTCONFIG], [fontconfig])
else
  AC_MSG_ERROR([required library fontconfig not found])
fi
LIBS="$LIBS $FONTCONFIG_LIBS"
CPPFLAGS="$CPPFLAGS $FONTCONFIG_CFLAGS"


dnl -------------
dnl Freetype2

if $PKG_CONFIG --exists freetype2 ; then
  PKG_CHECK_MODULES([FREETYPE], [freetype2])
else
  AC_MSG_ERROR([required library freetype2 not found])
fi
LIBS="$LIBS $FREETYPE_LIBS"
CPPFLAGS="$CPPFLAGS $FREETYPE_CFLAGS"


dnl -------------
dnl Cairo
dnl
dnl No accurate information about required cairo version. 
dnl Tested to work down to 1.2.4.
dnl Postscript support may need higher version

if $PKG_CONFIG --exists cairo ; then
  PKG_CHECK_MODULES([CAIRO], [cairo >= 1.2.4])
else
  AC_MSG_ERROR([required library cairo not found])
fi
LIBS="$LIBS $CAIRO_LIBS"
CPPFLAGS="$CPPFLAGS $CAIRO_CFLAGS"

dnl -------------
dnl POSIX threads 

ACX_PTHREAD()
LIBS="$LIBS $PTHREAD_LIBS $PTHREAD_CFLAGS"
CPPFLAGS="$CPPFLAGS $PTHREAD_CFLAGS"

dnl AC_CHECK_HEADERS([pthread.h])
dnl if $HAVE_PTHREAD_H ; then
dnl   LIBS="$LIBS -lpthread"
dnl else
dnl   AC_MSG_ERROR([required library pthread not found])
dnl fi

dnl -------------
dnl GSL

if $PKG_CONFIG --exists gsl ; then
  PKG_CHECK_MODULES([GSL], [gsl >= 1.12])
else
  AC_MSG_ERROR([required version of gsl not found])
fi
LIBS="$LIBS $GSL_LIBS"
CPPFLAGS="$CPPFLAGS $GSL_CFLAGS"

dnl AC_CHECK_LIB([gslcblas],[cblas_dgemm])
dnl AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

dnl -------------
dnl getcwd

AC_CHECK_HEADERS([unistd.h])
AC_CHECK_FUNCS([getcwd])

dnl -------------
dnl siginfo_t

AC_CHECK_TYPES([siginfo_t], [], [], [[#include <signal.h>]])
AC_CHECK_TYPES([clockid_t], [], [], [[#include <time.h>]])
AC_CHECK_TYPES([struct timespec], [], [], [[#include <time.h>]])

dnl -------------
dnl Librt

AC_CHECK_LIB(rt, clock_gettime)

dnl -------------
dnl Time

AC_CHECK_HEADERS([sys/time.h])
AC_HEADER_TIME
AC_CHECK_FUNCS([gettimeofday])
AC_CHECK_FUNCS([clock_gettime])

dnl -------------
dnl strerror_r

AC_FUNC_STRERROR_R

dnl ---------------------------------
dnl Output
dnl ---------------------------------

AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile doc/Makefile doc/Doxyfile doc/Doxyfile_web ibsimu-1.0.6dev.pc])
AC_OUTPUT

echo ""
echo "Optional compilation"
echo "--------------------"
echo ""
AS_IF([test "x$USEGTK3" = "xyes"],
      [echo "Gtk plotter     = yes"],
      [echo "Gtk plotter     = no"])
AS_IF([test "x$USEOPENGL" = "xyes"],
      [echo "OpenGL renderer = yes"],
      [echo "OpenGL renderer = no"])
AS_IF([test "x$USEUMFPACK" = "xyes"],
      [echo "UMFPACK solver  = yes"],
      [echo "UMFPACK solver  = no"])
AS_IF([test "x$USECSG" = "xyes"],
      [echo "CSG solid       = yes"],
      [echo "CSG solid       = no"])
echo ""
