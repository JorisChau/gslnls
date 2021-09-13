
AC_INIT([gslnls], 1.0.0)

AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

#if test "x$ac_cv_lib_gsl_gsl_blas_dgemm" != xyes; then
#    echo "------------------------- LIBGSL ERROR ------------------------------"
#    echo "Configuration failed because libgsl was not found. Try installing:"
#    echo " * deb: libgsl-dev (Debian, Ubuntu, etc)"
#    echo " * rpm: gsl-devel (Fedora, CentOS, RHEL)"
#    echo " * csw: gsl_dev (Solaris)"
#    echo "If libgsl is already installed, check that 'pkg-config' is in your"
#    echo "PATH and PKG_CONFIG_PATH contains a libgsl .pc file. If pkg-config"
#    echo "is unavailable you can set INCLUDE_DIR and LIB_DIR manually via:"
#   echo "R CMD INSTALL --configure-vars='INCLUDE_DIR=... LIB_DIR=...'"
#   echo "--------------------------------------------------------------------"
#  exit 1;
#fi

## credits to Dirk Eddelbuettel and RcppGSL
AC_PATH_PROG([GSL_CONFIG], [gsl-config])
if test "${GSL_CONFIG}" != ""; then
    # Use gsl-config for header and linker arguments (without BLAS which we get from R)
    GSL_CFLAGS=`${GSL_CONFIG} --cflags`
    GSL_LIBS=`${GSL_CONFIG} --libs`
else
    AC_MSG_ERROR([gsl-config not found, verify GSL installation])
fi

AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

AC_OUTPUT(src/Makevars)
