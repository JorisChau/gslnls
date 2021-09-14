
AC_INIT([gslnls], 1.0.0)

AC_PROG_CC
AC_PATH_PROG([GSL_CONFIG], [gsl-config])

if test "${GSL_CONFIG}" != ""; then
    # Use gsl-config for header and linker arguments (without BLAS which we get from R)
    GSL_CFLAGS=`${GSL_CONFIG} --cflags`
    GSL_LIBS=`${GSL_CONFIG} --libs`
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

AC_CHECK_HEADERS([gsl/gsl_multifit_nlinear.h],,[AC_MSG_ERROR([Missing required GSL headers. Upgrade GSL to version >= 2.2])])

AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

AC_OUTPUT(src/Makevars)
