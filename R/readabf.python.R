#' Read .abf with the pyabf API
#' 
#' This is a wrapper around pyabf
#' @param abf peridic file
read.multisweep.pyth = function(file) {
  abf = pyabf$ABF(file)
  # refactor the previous python script into plain R
  
  getdf =  function(sweep, abfobj) {
    # R function to extract 1 sweep
    sweep = as.integer(sweep - 1)# need to do -1 because of typical python indexing
    abfobj$setSweep(sweep)
    
    tidyr::tibble(
      t = abfobj$sweepX,
      value = abfobj$sweepY,
      command = abfobj$sweepC,
      sweep = rep(i, abfobj$sweepPointCount)
    )
  }
  # make an ldply call to extract all
  as.list(1:abf$sweepCount) %>%
    plyr::ldply(getdf, abf) %>%
    tidyr::tibble()
}