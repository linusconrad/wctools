# This script contains functions to annotate AP shared by different protocols

#' Find likely Peak points of Action potentials
#' 
#' This function checks a timeseries of differentiated voltages for whether the sign
#' of the derivative changes or it is 0.
#' Filtering for those points combined with filtering for overshoot (Voltage > 0) gave good results in identifying AP.
#' @param dv A vector of differentiated voltage values (dV/dt)
#' @return A logical vector of length dV

findAP = function(dV)
{ #returns logical vector of length dV
  is.changepoint = base::sign(dV)*base::sign(dplyr::lag(dV)) == -1
  is.zero = dV == 0
  c(is.changepoint + is.zero) != 0
}

#' Fetch Action potentials
#' 
#' This function takes a timeseries dataframe (from abf import) with the voltage variable `Vvar`.
#' It returns a dataframe containing the timepoints of the peaks of AP (`tpeak`), their height (`Vpeak`) and position within the sweep (`APindex`).
#' It contains the time scales `tAP` and `t` in list columns that can be matched with the source data timeseries to align AP for plotting.
#' @param df Time series df from an abf
#' @param Vvar Charavter vector, name of the voltage variable within df
#' @return A data frame with variables described above, with list columns 
#' @export 
returnAPdf = function(df, Vvar) {
  copydf = df
  # calculate the relevant stuff
  copydf$dV = c(NA, (base::diff(copydf[[Vvar]]) / 1000 / (1 / 50000)))
  copydf$AP = findAP(copydf$dV)
  # write the summary AP dataframe, filter AP
  copydf %<>% filter(.data[[Vvar]] > 0, .data$AP == T)
  # Following calculations change whether this has sweeps or not
  if (is.null(copydf$sweep) == F)
    copydf %<>% group_by(.data$sweep)
  # Remove doubly identified AP (points close together)
  # Impose 1 ms "refractory period"
  # then add an index column
  copydf %<>%
    mutate(ISI = .data$t - dplyr::lag(.data$t)) %>%
    filter(.data$ISI > 0.001 | is.na(.data$ISI)) %>%
    mutate(
      APindex = base::seq_along(.data$AP),
      Vpeak = .data[[Vvar]],
      tpeak = .data$t
    ) 
  #do some cleaning then print
  copydf %<>%
    select(-starts_with("Im"), -.data$dV, -.data[[Vvar]])
  
  # For every AP create a df with an offset as well as the original timescale
  # this is to tile and plot AP side by side on the same scale
  # To this end split the time axis into segments between each event (AP peak)
  copydf %<>%
    mutate(cutpoint.post = lead(((
      .data$tpeak + lag(.data$tpeak)
    ) / 2)),
    cutpoint.pre = lag(.data$cutpoint.post))
  
  # to generate original timescale 
  add.tiled.t = function(tp, pre, post) {
    # The first and last AP will have NA as cutpoints
    # In this case replace with fixed values
    if (is.na(pre))
      pre = tp - 0.01 #default to 10 ms win
    if (is.na(post))
      post = tp + 0.01
    # return a timescale with the cutpoints
    seq(from = plyr::round_any(pre, 1/50000) + 1/50000, # remove 1 sample each to avoid overlaps
        to = plyr::round_any(post, 1/50000) - 1/50000,
        by = 1 / 50000)
  }
  # generate timescale centered on event
  add.tiled.tAP = function(tp, pre, post){
    add.tiled.t(tp,pre,post) - tp
  }
  
  copydf %>%
    dplyr::tibble() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      tAP = purrr::map(.data$tpeak, .f = add.tiled.tAP, .data$cutpoint.pre, .data$cutpoint.post),
      t = purrr::map(.data$tpeak, .f = add.tiled.t, .data$cutpoint.pre, .data$cutpoint.post)
    ) %>%
    dplyr::select(tidyr::starts_with("cutpoint"))
  
}