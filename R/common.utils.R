#' Load a periodic stimulation abf file into a tidy data structure
#' 
#' This function takes an .abf file and loads it into memory as a tidy tibble.
#' Information about sweep etc is written to another column.
#' This will not always be a very memory saving way to do it.
#' This works for Vclamp and Cclamp recordings with two timeseries variables.
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
df$t = rep(c(c(1:npersweep)*(1/SR)), nsweep)
df$sweep = rep(1:nsweep, each = npersweep)
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
  data = read.multisweep(abf, 50000)
  # omit all but first sweep to save time
  data %<>% dplyr::filter(.data$sweep == 1)
  
  pivot_longer(data,
               cols = c(1:2),
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


