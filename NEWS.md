# gslnls 1.3.0

* Missing starting values/ranges allowed in `gsl_nls()`
* `lower` and `upper` parameter constraints included in `gsl_nls()` 
* Added 3 regression problems from Bates & Watts (1988)
* Updated multi-start algorithm in `gsl_nls()`
* Added configure.win, cleanup.win and Makevars.win.in
* Removed old Makevars and Makevars.win
* Several minor changes

# gslnls 1.2.0

* Added multi-start algorithm to `gsl_nls()`
* Added 56 NLS regression and optimization test problems
* Added unit tests in folder `unit_tests`
* Several minor changes/fixes

# gslnls 1.1.1

* Clean exits `gsl_nls()` and `gsl_nls_large()` when interrupted
* Default algorithm in `gsl_nls_large()` set to `"lm"`

# gslnls 1.1.0

* Added large-scale NLS regression with `gsl_nls_large()`

# gslnls 1.0.2

* Added Makevars.ucrt

# gslnls 1.0.1

* Initial release
