.onLoad <- function(libname, pkgname) {
  pyabf <<- reticulate::import("pyabf", delay_load = T)
  print("This package uses the python import `pyabf`. Before using, activate a Conda environment or other Python istallation with the pyabf module available. Refer to the `reticulate docs` to do this")
}