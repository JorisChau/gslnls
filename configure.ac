
AC_INIT([gslnls],[1.0.1])

AC_PROG_CC

AC_PATH_PROG([GSL_CONFIG], [gsl-config])

if test "${GSL_CONFIG}" != ""; then
    # Use gsl-config for header and linker arguments (without BLAS which we get from R)
    GSL_CFLAGS=`${GSL_CONFIG} --cflags`
    GSL_LIBS=`${GSL_CONFIG} --libs`
    CFLAGS="$CFLAGS $GSL_CFLAGS"
else
    echo "------------------------- LIBGSL ERROR ------------------------------"
    echo "Configuration failed because libgsl was not found. Try installing:"
    echo " * deb: libgsl-dev (Debian, Ubuntu, etc)"
    echo " * rpm: gsl-devel (Fedora, CentOS, RHEL)"
    echo " * csw: gsl_dev (Solaris)"
    echo "If libgsl is already installed, check that 'pkg-config' is in your"
    echo "PATH and PKG_CONFIG_PATH contains a libgsl .pc file. If pkg-config"
    echo "is unavailable you can set INCLUDE_DIR and LIB_DIR manually via:"
    echo "R CMD INSTALL --configure-vars='INCLUDE_DIR=... LIB_DIR=...'"
    echo "--------------------------------------------------------------------"
  exit 1;
fi

AC_MSG_CHECKING([if GSL version >= 2.2])
AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <string.h>
#include <gsl/gsl_version.h>
#include <stdio.h>
#include <stdlib.h>
int main() {
#ifdef GSL_VERSION
   int major, minor;
   char *gslv = GSL_VERSION;
   if ((sscanf(gslv, "%d.%d", &major, &minor)) != 2) {
     exit (1);
   }
   exit(major < 2 || minor < 2);
#else
  exit(1);
#endif
}
]])],
[gsl_version_ok=yes],
[gsl_version_ok=no],
[gsl_version_ok=yes])
if test "${gsl_version_ok}" = no; then
  AC_MSG_ERROR([Missing required GSL headers. Upgrade GSL to version >= 2.2])
else
  AC_MSG_RESULT([yes])
fi

AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
