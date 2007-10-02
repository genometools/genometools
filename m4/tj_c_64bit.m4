dnl @synopsis TJ_C_64BIT
dnl
dnl Provides a test for the existance of the long long int type and
dnl defines HAVE_LONG_LONG if it is found.
dnl Modified by Thomas.Jahns <Thomas.Jahns@epost.de> to include 3rd
dnl description parameter
dnl
dnl @version $Id: tj_c_64bit.m4,v 1.1 2007/07/03 18:05:31 tjahns Exp $
dnl @author Thomas Jahns <Thomas.Jahns@gmx.net>
dnl
AC_DEFUN([TJ_C_64BIT],
[AC_REQUIRE([AC_TYPE_LONG_LONG_INT])dnl
AC_REQUIRE([AC_TYPE_LONG_LONG_INT])dnl
AC_CACHE_CHECK([for 64 bit integral type], [tj_cv_c_64bit],
dnl first try to get uint64_t from standard headers
[tj_cv_c_64bit=no
  AC_LANG_PUSH([C])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif  /* HAVE_STDINT_H */
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif /* HAVE_SYS_TYPES_H */
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif /* HAVE_INTTYPES_H */
], [uint64_t a; int64_t b; a = b = 1;])],dnl
tj_cv_c_64bit=yes,dnl
[AC_MSG_NOTICE([no native (u)int64_t define.])
  AC_MSG_CHECKING([wether long long replacement is available])
  if test "$ac_cv_type_long_long_int" = yes; then
    tj_cv_c_64bit=llreplace
  else
    tj_cv_c_64bit=no
  fi
])
AC_LANG_POP([C])])
   if test $tj_cv_c_64bit = llreplace ; then
   if test $ac_cv_type_long_long_int = yes -a $ac_cv_c_int64_t = no; then
     AC_MSG_NOTICE([using long long replacement])
     AC_DEFINE([int64_t], [long long],
       [default int64_t to long long])
   else
     tj_cv_c_64bit=noreplace
   fi
   if test $tj_cv_c_64bit = llreplace -a $ac_cv_type_unsigned_long_long_int = yes -a $ac_cv_c_uint64_t = no; then
     AC_MSG_NOTICE([using long long replacement])
     AC_DEFINE([uint64_t], [unsigned long long],
       [default uint64_t to unsigned long long])
     tj_cv_c_64bit=yes
   else
     tj_cv_c_64bit=no
   fi
   fi
   if test $tj_cv_c_64bit = yes; then
     AC_DEFINE(HAVE_64BIT_TYPE,,[System provides support for (u)int64_t])
   fi
])
