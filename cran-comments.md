## CRAN package version 1.0.0

* System requirements: GSL (>= 2.2)

### CRAN-team comments

* You could omit the single quotes in the description text.
* If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form
    authors (year) <doi:...>
    authors (year) <arXiv:...>
    authors (year, ISBN:...) 
* Missing Rd-tags in up to 11 .Rd files, e.g.:
     anova.gsl_nls.Rd: `\value`
     coef.gsl_nls.Rd: `\value`
     deviance.gsl_nls.Rd: `\value`
     df.residual.gsl_nls.Rd: `\value`
     fitted.gsl_nls.Rd: `\value`
     formula.gsl_nls.Rd: `\value`
     ... 

### Edits

* Edited DESCRIPTION file to omit single quotes.
* Included reference in DESCRIPTION File to the GNU Scientific Library.
* Edited .Rd files to include missing `\value` tags.

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
