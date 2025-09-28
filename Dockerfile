FROM rocker/r-devel

ARG PKG_VER=1.4.2
ARG GSL_VER=2.8

ADD https://ftp.gnu.org/gnu/gsl/gsl-${GSL_VER}.tar.gz gsl-${GSL_VER}.tar.gz
RUN tar -xf gsl-${GSL_VER}.tar.gz

WORKDIR /gsl-${GSL_VER}
RUN ./configure \
    && make \
    && make install
WORKDIR /

COPY gslnls_${PKG_VER}.tar.gz .
RUN Rdevel CMD INSTALL gslnls_${PKG_VER}.tar.gz

CMD ["Rdevel", "-e", "source(system.file('unit_tests', 'unit_tests_gslnls.R', package = 'gslnls'))"]
