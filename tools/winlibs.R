# Download GSL 2.7 Rtools build
# credits to Dirk Eddelbuettel and rccpziggurat
if (!file.exists("../windows/gsl-2.7/include/gsl/gsl_blas.h")) {
    download.file("https://github.com/rwinlib/gsl/archive/v2.7.zip", "lib.zip", quiet = TRUE)
    dir.create("../windows", showWarnings = FALSE)
    unzip("lib.zip", exdir = "../windows")
    unlink("lib.zip")
}
