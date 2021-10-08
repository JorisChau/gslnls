## CRAN package version 1.0.1

* System requirements: GSL (>= 2.2)

### Changes

* Fixed `gsl_linalg_QRPT_rcond` memory error detected by valgrind in CRAN check results

### Test environments

* ubuntu gcc R-release, R-devel (rhub)
* debian gcc R-release, R-devel, R-patched (rhub)
* debian clang R-devel (rhub)
* fedora clang/gcc R-devel (rhub)
* centos8-epel R-4.0.4 (rhub)
* macos-highsierra R-release (rhub)
* solaris10 gcc R-release/GSL-2.4 (local vm)
* solaris10 ods R-release/GSL-2.4 (local vm)
* win-builder R-release, R-devel, R-old-release (see below)

### Compiled code checks

* ubuntu-rchk
* ubuntu clang R-release --use-valgrind 
* ubuntu clang R-release --use-gct

### Notes

* Possibly mis-spelled words in DESCRIPTION:
    GSL (3:8, 7:101)
    Levenberg (7:154)
    Marquadt (7:164)
