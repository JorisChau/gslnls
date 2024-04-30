
AC_INIT([gslnls],[1.3.2],[https://github.com/JorisChau/gslnls/issues/])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`
AC_PROG_CC

FC=`"${R_HOME}/bin/R" CMD config FC`
FCLAGS=`"${R_HOME}/bin/R" CMD config FFLAGS`
AC_PROG_FC

AC_PATH_PROG([GSL_CONFIG], [gsl-config])

if test -n "${GSL_CONFIG}"; then
    # Use gsl-config for header and linker arguments
    GSL_CFLAGS=`${GSL_CONFIG} --cflags`
    GSL_LIBS=`${GSL_CONFIG} --libs`
    CFLAGS="$CFLAGS $GSL_CFLAGS"
else
    echo "------------------------- LIBGSL ERROR ------------------------------"
    echo "Configuration failed because libgsl was not found. Try installing:"
    echo " * deb: libgsl-dev (Debian, Ubuntu, etc)"
    echo " * rpm: gsl-devel (Fedora, CentOS, RHEL)"
    echo " * csw: gsl_dev (Solaris)"
    echo "If libgsl is installed, check that 'gsl-config' is in your PATH. If"
    echo "gsl-config is unavailable, set INCLUDE_DIR and LIB_DIR manually via:"
    echo "R CMD INSTALL --configure-vars='INCLUDE_DIR=... LIB_DIR=...'"
    echo "--------------------------------------------------------------------"
  exit 1;
fi

# check gsl version
AC_MSG_CHECKING([if GSL version >= 2.2])
AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <string.h>
#include <gsl/gsl_version.h>
#include <stdio.h>
#include <stdlib.h>
int main(void) {
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
