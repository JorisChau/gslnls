## CRAN package version 1.4.1

* System requirements: GSL (>= 2.2)

## Author comments

* Fixed compatibility with older GSL versions (<2.5) 
  causing installation ERRORs in CRAN checks	


## Test environments

* ubuntu gcc R-release + GSL 2.2.1 (local install)
* ubuntu gcc R-oldrel, R-release, R-next, R-devel (rhub, gh-actions)
* debian gcc/clang R-devel (rhub)
* fedora gcc/clang R-devel (rhub)
* macos-14 clang R-release, R-next, R-devel (gh-actions, rhub)
* windows gcc R-oldrel, R-release, R-next, R-devel (r-winbuilder, gh-actions)

## Compiled code checks

* ubuntu-rchk R-devel (docker)
* fedora gcc R-devel --use-valgrind (rhub)
* ubuntu gcc R-4.4 --use-gct (local install)
