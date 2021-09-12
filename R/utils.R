gsl_config <- function(show = TRUE) {  
	if(.Platform$OS.type == "unix") {
		if(Sys.which("gsl-config") == "") {
			if(isTRUE(show)) {
				cat("No gsl library info available\n")	
			}
			return(FALSE)
		} else 
		if(isTRUE(show)) {
			args <- c(prefix = "--prefix", libs = "--libs", cflags = "--cflags", version = "--version")
			values <- vapply(args, system2, character(1), command = "gsl-config", stdout = TRUE)
			cat("gsl-config [", Sys.which("gsl-config"), "]\n", sep = "")
			cat("    prefix:", values[["prefix"]], "\n")
			cat("    libs:", values[["libs"]], "\n")
			cat("    cflags:", values[["cflags"]], "\n")
			cat("    version:", values[["version"]], "\n")
		}
		return(TRUE)    
	} else {
		## only Rtools40
		if(Sys.which("pacman") == "" && Sys.getenv("LIB_GSL") == "") {
			if(isTRUE(show)) {
				cat("No gsl library info available\n")
			}
			return(FALSE)
		} else if(Sys.getenv("LIB_GSL") != "") {
			loc <- list.files(file.path(Sys.getenv("LIB_GSL"), "include"), pattern = "gsl_version", recursive = TRUE)
			if(length(loc) > 0) {
				if(isTRUE(show)) {
					cat("Detected gsl library files from 'LIB_GSL' variable [", Sys.getenv("LIB_GSL"), "]\n", sep = "")
				}
				return(TRUE)
			} else {
				if(isTRUE(show)) {
					cat("No gsl library info available\n")
				}
				return(FALSE)
			}
		} else {
			out <- system2("pacman", args = c("-Qi", "mingw-w64-{i686,x86_64}-gsl"), stdout = TRUE, stderr = TRUE)
			if(isTRUE(show)) {
				cat(out, sep = "\n")
			}
			return(is.null(attr(out, "status")))
		}
	}
}
