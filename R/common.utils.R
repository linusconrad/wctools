#' Load a periodic stimulation abf file into a tidy data structure
#' 
#' This function takes an .abf file and loads it into memory as a tidy tibble.
#' Information about sweep etc is written to another column.
#' This will not always be a very memory saving way to do it.
#' This works for Vclamp and Cclamp recordings with two timeseries variables.
#' CAREFUL: This takes the recorded stimulus trace, not the applied stimulus ("add stimulus waveform" of clampfit)
#' Advantage of this function is that it should work with complex waveform protocols (AP waveform, ZAP) without the stimulus file being present.
#' 
#' @param abf Path to the input .abf file (periodic abf, two timeseries)
#' @param SR The sampling rate of the recording in Hz
#' @return A tibble object of the file with variables, time and sweep columns
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @import tidyr
#' @export
read.multisweep = function(abf, SR){
data = readABF::readABF(abf)
# assign some metadata
# number of sweeps
nsweep = length(data$data)
# sample per sweep
npersweep = data$header$sweepLengthInPts
# Build a tidy tIV frame
df = as_tibble(do.call(rbind, data$data), .name_repair = "unique")
names(df) = data$channelNames

# create sweep and time signatures
df$t = rep(c(c(0:(npersweep-1))*(1/SR)), nsweep)
df$sweep = rep(1:nsweep, each = npersweep)
# reorder columns
df %<>% dplyr::relocate(., tidyr::starts_with("t"), .before = tidyr::starts_with("V"))
# return the df
df
}

#' Create a timeseries graph from an abf
#' 
#' I dont have a way of browsing and displaying .abf on my linux machine.
#' Therefore this function just creates a plot of the first sweep so I can see what protocol the file is.
#' ...Clampex misses the option of writing the protocol name to file...
#' Sampling rate defaults to 50000 Hz.
#' 
#' @param abf An .abf file to create a plot of
#' @return A ggplot object of the first sweep with voltage and current trace
#' @import ggplot2
#' @import tidyr
#' @importFrom rlang .data
#' @export
draw.first.sweep = function(abf){
  # use the python implementation to read but fallback to R if it fails for some reason
  data = try(read.multisweep.pyth(abf))
  
  # checks whether stimulus variable is empty and loads with R in that case
  if (is.na(data[3][[1]][1])) 
    data = read.multisweep(abf, 50000)
  
  # omit all but first sweep to save time
  data %<>% dplyr::filter(.data$sweep == 1)
  
  pivot_longer(data,
               cols = c(2:3),
               names_to = "variable",
               values_to = "timeseries") %>%
    ggplot(aes(
      x = .data$t,
      y = .data$timeseries,
      colour = .data$variable
    )) +
    labs(y = " ",
         x = "t, s") +
    geom_line() +
    scale_x_continuous(expand = c(0, 0)) +
    geom_hline(linetype = 3,
               colour = "grey50",
               yintercept = 0) +
    ggthemes::scale_color_few() +
    facet_wrap(
      ~ variable,
      ncol = 1,
      scales = "free_y",
      strip.position = "left",
      as.table = F
    ) +
    theme(legend.position = "none")
}
  
#' Write a timeseries plot to file
#' 
#'  This is the same function as "plot.first.sweep" but instead of returning the plot object it just writes a .png to file
#'  in the directory of the .abf file.
#'  This function is meant to be used in a l_ply call of a list of files.
#'  
#' @param abf An .abf file to create a plot of
#' @return does not return a thing, just writes a png to file
#' @import ggplot2
#' @import tidyr
#' @export
draw.first.sweep.png = function(abf){
  plot = draw.first.sweep(abf)
  
  ggsave(filename = paste0(abf, "first.png"), plot = plot, width = 6, height = 4)
}

#' Take a measurement of a timeseries at a given timewindow
#' 
#' This function takes a measurement at a given timepoint.
#' var being a vector in a dataframe with a time variable
#' tp being the timepoint, win being a window around it
#' this is meant for dplyr::summary calls so is supposed to return one datum
#' 
#' @param var the name of the variable
#' @param tvar the name of the time variable
#' @param tp the timepoint to measure
#' @param win the time window to take the measurement of, (centered on tp)
#' @export
measure.timepoint = function(var, tvar, tp, win) {
  mean(var[tvar > (tp - win / 2) &
             tvar < (tp + win / 2)])
}

# ggplot options that are nice
plot.font = "sans"   # activate as needed based on OS

#' Default setting for ggplot
#' Nice setting for ggplot
#' @export
theme_linus = 
 ggthemes::theme_few() +
    theme(
      legend.position = "bottom",
      legend.title = NULL,
      strip.placement = "outside",
      strip.background = element_blank(),
      #          plot.margin = margin(t = 0.2,r =0.1, b = 0, l = 0, "cm"),
      text = element_text(
        family = plot.font,
        size = 11,
        face = "plain"
      ),
      title = element_text(family = plot.font),
      strip.text = element_text(family = plot.font)
    )

#' Read .abf with the pyabf API
#' 
#' This is a wrapper around pyabf.
#' This function works on periodic stimulation abf files.
#' The original protocol stimulus is reconstructed and added to the file.
#' The function will fail on files with arbitrary waveform specification.
#' Here, the stimulus waveform file needs to be present in the project path or the directory of the rawdata file.
#' @param abf periodic stimulation abf file
#' @return tibble of the abf with stimulus and sweep variable
#' @export
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
      sweep = rep(sweep, abfobj$sweepPointCount)
    )
  }
  # get the names of the variables from the abf object
  # they are hidden away a bit
  nameC = abf$dacNames %>% stringi::stri_split_regex(pattern = " ") %>% unlist() %>% .[1]
  nameY = abf$adcNames %>% stringi::stri_split_regex(pattern = " ") %>% unlist() %>% .[1]
  # make an ldply call to extract all
  df =
    as.list(1:abf$sweepCount) %>%
    plyr::ldply(getdf, abf) %>%
    tidyr::tibble() %>%
    mutate(sweep = .data$sweep +1) # return to R indexing
  # clean the names
  names(df)[2] = nameY
  names(df)[3] = nameC
  # return the data
  df
}

#' Read a custom stimulus abf
#' 
#' This is a wrapper around pyabf.
#' This function works on custom waveform stimulation (file-based) abf files.
#' The original protocol stimulus is reconstructed and added to the file.
#' This works with the ZAP protocol I designed which gets truncated by pclamp in a particular way, this reading function fixes that.
#' I have not tested with other protocols.
#' @param abf periodic stimulation abf file
#' @return tibble of the abf with stimulus and sweep variable
#' @export
#' @importFrom magrittr %>%
read.multisweep.custom = function(file) {
  abf = pyabf$ABF(file)
  # refactor the previous python script into plain R
  getdf =  function(sweep, abfobj) {
    # R function to extract 1 sweep
    sweep = as.integer(sweep - 1)# need to do -1 because of typical python indexing
    abfobj$setSweep(sweep)
    
    # default stimulus deadtime
    deadtime = max(abfobj$sweepX) / 64
    deadsamples = round(length(abfobj$sweepX) / 64)
    # 0 to be filled
    missingsamples = length(abfobj$sweepX) - length(abfobj$sweepC)
    
    tidyr::tibble(
      t = abfobj$sweepX,
      value = abfobj$sweepY,
      command = c(abfobj$sweepC,
                  rep.int(0, missingsamples)),
      sweep = rep(sweep, abfobj$sweepPointCount)
    )
  }
  # get the names of the variables from the abf object
  # they are hidden away a bit
  nameC = abf$dacNames %>% stringi::stri_split_regex(pattern = " ") %>% unlist() %>% .[1]
  nameY = abf$adcNames %>% stringi::stri_split_regex(pattern = " ") %>% unlist() %>% .[1]
  # make an ldply call to extract all
  df =
    as.list(1:abf$sweepCount) %>%
    plyr::ldply(getdf, abf) %>%
    tidyr::tibble() %>%
    mutate(sweep = .data$sweep + 1) # return to R indexing
  # clean the names
  names(df)[2] = nameY
  names(df)[3] = nameC
  # return the data but filter out the deadtime
  df %>% dplyr::filter(.data$t > (max(.data$t) / 64))
  #df
}












