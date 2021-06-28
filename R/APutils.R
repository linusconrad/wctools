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
#' @param Vvar Character vector, name of the voltage variable within df
#' @param thresh1 threshold for detection of AP peaks, will default to 0 (aP is an AP if its overshoots)
#' @return A data frame with variables described above, with list columns 
#' @import tidyr
#' @importFrom magrittr %<>%
#' @export 
returnAPdf = function(df, Vvar, thresh1) {
  copydf = df
  
  # default value for AP detection threshold
  if(missing(thresh1))
    thresh1 <- 0
    
  # calculate the relevant stuff
  copydf$dV = c(NA, (base::diff(copydf[[Vvar]]) / 1000 / (1 / 50000)))
  copydf$AP = findAP(copydf$dV)
  # write the summary AP dataframe, filter AP
  copydf %<>% filter(.data[[Vvar]] > thresh1, .data$AP == T)
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
      pre = tp - 0.015 #default to 15 ms win
    if (is.na(post))
      post = tp + 0.015
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
    dplyr::select(-tidyr::starts_with("cutpoint"))
  
}


#' Fortify an DF with AP information
#'
#' This function takes an abf timeseries dataframe and detects AP within it.
#' Then it adds the index of the AP within the sweep and adds a timescale to align all AP found by their peaks.
#' @param df Dataframe with the rawdata
#' @param Vvar name of the voltage variable in `df`, string
#' @param thresh2 threshold for detection of AP peaks, will default to 0 (AP is an AP if its overshoots)
#' @return source dataset with added AP index and AP centered timescale
#' @import tidyr
#' @importFrom magrittr %<>%
#' @export
addAP = function(df, Vvar, thresh2) {
  # default value for AP detection threshold
  if(missing(thresh2))
    thresh2 <- 0
  # change behaviour depending on whether is gap free or sweep data
  columns = c("t", "tAP", "APindex", "sweep")
  
  if (is.null(df$sweep))
    columns = c("t", "tAP", "APindex")
  
  APdf =
    returnAPdf(df, Vvar, thresh1 = thresh2) %>%
    dplyr::select(dplyr::all_of(columns)) %>%
    unnest(cols = c(.data$t, .data$tAP))
  #join and return
  dplyr::left_join(df, APdf)
}

#' Generate AP summary stats
#'
#' This function takes an abf timeseries, detects AP within it and returns their waveform measurements.
#' @param df Dataframe with the rawdata
#' @param vvar name of the voltage variable in `df`, string
#' @param thresh3 threshold for detection of AP peaks, will default to 0 (AP is an AP if its overshoots)
#' @return summary stats of the AP such as Peak, time, afterhyperpolarisation, width and interspike interval
#' @import tidyr
#' @importFrom magrittr %<>%
#' @export
getAPstats = function(df, vvar, thresh3) {
  if(missing(thresh3))
    thresh3 <- 0
  data = addAP(df, vvar, thresh2 = thresh3)
  
  # The definition of Vrest might need some work..
  Vrest = pracma::Mode(data[[vvar]])
  # Make a normalised Voltage scale

  data =
    data %>%
    mutate(V0 = .data[[vvar]] - Vrest)
  
  # ensure proper grouping
  if (is.null(data$sweep) == F)
    data %<>% group_by(.data$sweep, .data$APindex)
  else
    data %<>% group_by(.data$APindex)
  
  data = 
    data %>%
    mutate(Vnorm = .data$V0 / mean(.data$V0[.data$V0 == max(.data$V0)])) %>%
    select(-.data$V0) %>%
    # get rid of regions outside of APs
    filter(!is.na(.data$tAP))
  
  APstats =
    wctools::returnAPdf(df, vvar, thresh1 = thresh3) %>%
    select(-.data$t,-.data$tAP)
  
  # extract the maximal upstroke velocity
  upstroke = 
    data %>%
    # this can find erroneous stimulus artefacts.
    # filter out the actual upstroke (time immediatly before the peak (lets say 2 ms))
    filter(.data$tAP > -0.002) %>%
    dplyr::mutate(dV = c(NA, (base::diff(.data[[vvar]]) / 1000 / (1 / 50000)))) %>% # unit is V/s
    #using max comes up with strange high values, use percentile
    # also here strange high values can crop up, just filter everything that is blatantly out of physiological range
    filter(.data$dV < 250)%>%
    dplyr::summarise(upstroke = stats::quantile(.data$dV, probs = 0.99, na.rm = T))
  
  APstats %<>%
    left_join(upstroke)
  
  # extract the afterhyperpolarisation
  afterhyp =
    data %>%
    filter(.data$tAP > 0) %>% #filter all data after AP peak
    dplyr::summarise(Vmin = min(.data[[vvar]]),
                     tmin = c(.data$tAP[.data[[vvar]] == min(.data[[vvar]])])) %>%
    #take first minimum value only
    mutate(rank = seq_along(.data$tmin)) %>%
    filter(.data$rank == 1) %>%
    #get rid of helper column
    select(-.data$rank)
  
  APstats %<>%
    left_join(afterhyp)
  
  # extract the AP width
  data %<>%
    mutate(Vnormtest = abs(.data$Vnorm - 0.5),
           APphase = sign(.data$tAP)) %>%
    filter(.data$APphase != 0)
  
  # ensure proper grouping
  # split the data into pre and post peak to find a value for V1/2 before and after!
  if (is.null(data$sweep) == F)
    data %<>% group_by(.data$sweep, .data$APindex, .data$APphase)
  else
    data %<>% group_by(.data$APindex, .data$APphase)
  
  width =
    data %>%
    # filter crossings that belong to a subsequent or previous AP/ depolarisation
    left_join(afterhyp) %>%
    dplyr::filter(.data$tAP < .data$tmin) %>%
    # extract the values of V1/2
    summarise(halfpoint = .data$tAP[.data$Vnormtest == min(.data$Vnormtest)],
              #testinghalfpoint = min(abs(.data$tAP[.data$Vnormtest == min(.data$Vnormtest)])),
              Vhalf = .data[[vvar]][.data$Vnormtest == min(.data$Vnormtest)]) 
  
  # now calculate the width proper, regrouping is required
  width %<>%  ungroup(.)
  
  if (is.null(width$sweep) == F)
    width %<>% group_by(.data$sweep, .data$APindex)
  else
    width %<>% group_by(.data$APindex)
  
  # calculate the width 
  width %<>%
    mutate(width = .data$halfpoint - lag(.data$halfpoint)) %>%
       # widthtest  = cumsum(testinghalfpoint)) %>%
    filter(!is.na(.data$width)) %>%
    # filter out impossible values (AP tiling is 15 +15 so anything larger is nonsense)
    # such spurious mistake happen rarely but I dont know in which precise context that bugs crop up
    # safest is to just discard, any cell will have many AP such that one individual parameter wont matter
    filter(.data$width < 0.02)
    #select(.data$width)
  
  APstats %<>%
    left_join(width)
  
  ##calculating ISI
  # regroup
  if (is.null(APstats$sweep) == F)
    APstats %<>%
    ungroup() %>%
    group_by(.data$sweep)
  else
    APstats %<>% ungroup()
  
  APstats %>%
    mutate(ISI = .data$tpeak - dplyr::lag(.data$tpeak))
}


