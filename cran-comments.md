## CRAN package version 1.3.2

* System requirements: GSL (>= 2.2)

## Comments

* Fix errors cran check results (fedora builds)
* Removed configure.win + Makevars.win.in 
* Added static Makevars.win as recommended by T. Kalibera
* Fix prototype warning in config.log

## Test environments

* ubuntu gcc R-oldrel, R-release, R-next, R-devel (rhub, gh-actions)
* debian gcc/clang R-devel (rhub)
* fedora gcc/clang R-devel (rhub)
* macos-monterey clang R-release, R-next (gh-actions)
* windows gcc R-oldrel, R-release, R-next, R-devel (r-winbuilder, gh-actions)

## Compiled code checks

* ubuntu-rchk R-devel (docker)
* fedora gcc R-devel --use-valgrind (rhub)
* ubuntu gcc R-4.3 --use-gct (local install)
