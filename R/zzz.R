.onLoad <- function(libname, pkgname) {
  reticulate::use_condaenv("analysis")
  pyabf <<- reticulate::import("pyabf")
  print("Imported pyabf with reticulate")
}