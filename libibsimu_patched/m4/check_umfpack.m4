dnl @synopsis CHECK_UMFPACK()
dnl
dnl This macro searches for an installed umfpack library. If nothing was
dnl specified when calling configure, it searches first in /usr/local
dnl and then in /usr. If the --with-umfpack=DIR is specified, it will try
dnl to find it in DIR/include/umfpack.h and DIR/lib64/libumfpack.so. If
dnl --without-umfpack is specified, the library is not searched at all.
dnl
dnl If either the header file (umfpack.h) or the library (libumfpack) is not
dnl found, the configuration exits on error, asking for a valid umfpack
dnl installation directory or --without-umfpack.

AC_DEFUN([CHECK_UMFPACK],

[

# Handle user hints
AC_ARG_WITH([umfpack],
  [AS_HELP_STRING([--with-umfpack],
                  [Enable umfpack solver @<:@default=check@:>@])],
  [],
  [with_umfpack=check])

# Locate umfpack
UMFPACK_HOME=
if test "$with_umfpack" != "no"; then
	if test "$with_umfpack" != "yes" && test "$with_umfpack" != "check"; then
		# Check user given path
		UMFPACK_HOME="$with_umfpack"
		AC_MSG_CHECKING(umfpack.h in ${UMFPACK_HOME}/include/suitesparse/)
		if test ! -f "${UMFPACK_HOME}/include/suitesparse/umfpack.h"; then
			AC_MSG_RESULT(no)
			UMFPACK_HOME=
		else
			AC_MSG_RESULT(found)
		fi
	fi

	if test ! -n "${UMFPACK_HOME}"; then
		# Check default paths
		UMFPACK_HOME=/usr/local
		AC_MSG_CHECKING(umfpack.h in ${UMFPACK_HOME}/include/suitesparse/)
		if test ! -f "${UMFPACK_HOME}/include/suitesparse/umfpack.h"; then
			AC_MSG_RESULT(no)
			UMFPACK_HOME=/usr
			AC_MSG_CHECKING(umfpack.h in ${UMFPACK_HOME}/include/suitesparse/)
			if test ! -f "${UMFPACK_HOME}/include/suitesparse/umfpack.h"; then
				AC_MSG_RESULT(no)
				UMFPACK_HOME=
				AC_SUBST([UMFPACK_LIBS], [])
				AC_SUBST([UMFPACK_CFLAGS], [])
			else
				AC_MSG_RESULT(found)
			fi
		else
			AC_MSG_RESULT(found)
		fi
	fi
fi


if test -n "${UMFPACK_HOME}"; then
        UMFPACK_OLD_LDFLAGS=$LDFLAGS
        UMFPACK_OLD_CPPFLAGS=$CPPFLAGS
        LDFLAGS="$LDFLAGS -L${UMFPACK_HOME}/lib64"
        CPPFLAGS="$CPPFLAGS -I${UMFPACK_HOME}/include/suitesparse"
        AC_LANG_SAVE
        AC_LANG_C
	AC_CHECK_LIB(umfpack, umfpack_di_symbolic, [umfpack_cv_libumfpack=yes], [umfpack_cv_libumfpack=no], [-lamd -lblas])
        AC_CHECK_HEADER(umfpack.h, [umfpack_cv_umfpack_h=yes], [umfpack_cv_umfpack_h=no])
	AC_MSG_CHECKING(umfpack.h)
	AC_LANG_RESTORE
        if test "$umfpack_cv_libumfpack" = "yes" -a "$umfpack_cv_umfpack_h" = "yes"; then
		AC_MSG_RESULT(ok)
                AC_CHECK_LIB(blas, dgemm_)
                AC_CHECK_LIB(amd, amd_1)
                AC_CHECK_LIB(umfpack, umfpack_di_symbolic)
		USEUMFPACK=yes
		AC_SUBST([UMFPACK_LIBS], ["-L${UMFPACK_HOME}/lib64 -lumfpack -lamd -lblas"])
		AC_SUBST([UMFPACK_CFLAGS], ["-I${UMFPACK_HOME}/include/suitesparse"])
        else
		AC_MSG_RESULT(failed)
                LDFLAGS="$UMFPACK_OLD_LDFLAGS"
                CPPFLAGS="$UMFPACK_OLD_CPPFLAGS"
                AC_MSG_ERROR(either specify a valid umfpack installation with --with-umfpack=DIR or disable umfpack usage with --without-umfpack)
        fi
else
	AC_SUBST([UMFPACK_LIBS], [""])
	AC_SUBST([UMFPACK_CFLAGS], [""])
fi

])
