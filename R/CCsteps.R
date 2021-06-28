#' Process an abf file of the classic CC step protocol
#' 
#' This function takes an .abf file of the passive protocol and returns to the folder a summary graph and a .csv.
#' These contain per-sweep summaries.
#' SR of 50000 is assumed.
#' Function is tailored to this protocol only.
#' @param abffile The File
#' @param Vjunc Junction potential to add
#' @param thresh4 Voltage threshold of the AP detector, defaults to 0
#' @return a .png and a .csv with the summarised analysis. Named with the filename trunk of 'abffile'.
#' @import ggplot2
#' @import tidyr
#' @importFrom magrittr %<>%
#' @import patchwork
#' @import dplyr
#' @importFrom rlang .data
#' @export
process.CCsteps = function(abffile, Vjunc, thresh4){
  if(missing(thresh4))
    thresh4 <- 0
  # read the data
  CCdata = read.multisweep.pyth(abffile)
  # assign stim onset
  tstim = 0.03436
  
  # standardise the variable names
  names(CCdata) = c("t", "Vm", "Istim", "sweep")
  # substract junction
  CCdata =
    CCdata %>%
    dplyr::mutate(Vm = .data$Vm + Vjunc)
  
  #get the APstats
  CCstepsummary =
    CCdata %>%
    getAPstats("Vm", thresh3 = thresh4) %>%
    mutate(latency = .data$tpeak - tstim,
           halfpoint1 = .data$width - .data$halfpoint,
           halfstart = .data$tpeak - .data$halfpoint1,
           halfend = .data$tpeak + .data$halfpoint) %>%
    select(-.data$halfpoint, -.data$halfpoint1)
  
  CCstepbysweep =
    # make a by sweep summary, this will help to extract rheobase etc
    CCdata %>% 
    group_by(.data$sweep, .data$Istim) %>%
    filter(.data$t > tstim, .data$t < 0.48436) %>%
    dplyr::summarise(plateau = pracma::Mode(.data$Vm))
  
  #get the amount of AP out
  CCstepbysweep =
    left_join(
      CCstepbysweep,
      CCstepsummary %>%
        group_by(.data$sweep) %>%
        summarise(nAP = max(.data$APindex),
                  Istim = .data$Istim[1])
    ) %>%
    replace_na(list(nAP = 0))
  
  
  # plot the AP detection for quality control
  plot.QC = CCdata %>%
    ggplot(aes(x = .data$t, y = .data$Vm)) +
    geom_hline(yintercept = c(0),
               colour = "grey50",
               linetype = 3) +
    geom_hline(yintercept = c(thresh4),
               colour = "blue",
               linetype = 3) +
    facet_wrap( ~ sweep,
                ncol = 5,
                as.table = T) +
    geom_vline(
      data = CCstepsummary,
      aes(xintercept = .data$tpeak),
      colour = "blue",
      alpha = 0.5
    ) +
    geom_hline(
      data = CCstepbysweep,
      aes(yintercept = .data$plateau),
      colour = "cyan",
      alpha = 0.5
    ) +
    geom_line(size = 0.2) +
    geom_segment(
      data = CCstepsummary,
      mapping = aes(
        x = .data$halfstart,
        y = .data$Vhalf,
        xend = .data$halfend,
        yend = .data$Vhalf
      ),
      colour = "red"
    ) +
    cowplot::theme_nothing()
  
  # plot the APmax
  plotAPmax =
    CCstepbysweep %>%
    ggplot(aes(x = .data$Istim, y = .data$nAP)) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(0,NA)) +
    geom_step() +
    labs(y = "AP Count",
         x = latex2exp::TeX("I$_{stim}$, pA"))
  
  # Plot all kinds of values by sweep/stim
  plot.bysweep = 
    CCstepsummary %>%
    mutate(# rescales these variables to ms
      tmin = .data$tmin * 1000,
      width = .data$width * 1000,
      latency = .data$latency * 1000,
      ISI = .data$ISI * 1000
    ) %>% 
    pivot_longer(cols = c("Vpeak", "Vmin", "width", "ISI", "tmin", "latency", "upstroke")) %>%
    mutate(sweepf = as.factor(.data$sweep)) %>%
    ggplot(aes(x = .data$Istim, y = .data$value, colour = .data$APindex, fill = .data$APindex)) +
    labs(x = latex2exp::TeX("I$_{stim}$, pA")) +
    facet_wrap(~name, strip.position = "left", scales = "free_y") +
    geom_line(linetype = 3, aes(group = .data$APindex)) +
    geom_point(alpha = 0.5, shape = 21) +
    theme(axis.title.y = element_blank()) +
    theme(legend.position = "none")
  
  plot.APsignature = 
    CCstepsummary %>%
    ggplot(aes(y = .data$Istim, x = .data$latency *1000, colour = .data$APindex)) +
    labs(x = "Latency, ms",
         y = latex2exp::TeX("I$_{stim}$, pA")) +
    geom_line(aes(group = .data$APindex), linetype = 1, alpha = 0.2)+
    geom_point(shape = 4) 
  
  
  # Plot all kinds of values by APindex
  plot.byindex = 
    CCstepsummary %>%
    mutate(# rescales these variables to ms
      tmin = .data$tmin * 1000,
      width = .data$width * 1000,
      latency = .data$latency * 1000,
      ISI = .data$ISI * 1000
    ) %>% 
    pivot_longer(cols = c("Vpeak", "Vmin","latency", "width", "ISI", "tmin", "upstroke")) %>%
    mutate(indexf = as.factor(.data$APindex)) %>%
    #filter(.data$APindex %in% c(1,2)) %>%
    ggplot(aes(x = .data$indexf, y = .data$value)) +
    labs(x = "AP Index") +
    facet_wrap(~name, strip.position = "left", scales = "free_y") +
    geom_line(aes(group = sweep), colour = "grey50") +
    geom_point(shape = 21) +
    theme(axis.title.y = element_blank())
  
  # make a multi plot page with patchwork
  summaryplot = 
    ((plot.QC / (plot.APsignature + plotAPmax))/ (plot.bysweep)/plot.byindex) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(title = "Classic CC Steps",
                               subtitle = paste0(abffile, "\n Analysed on ", Sys.Date()))
  
  ggsave(summaryplot, width = 8.25, height = 11.75, filename = paste0(abffile, ".CCsteps.png"))
  
  readr::write_csv(CCstepsummary, file = paste0(abffile, ".CCstepsAP.csv"))
  readr::write_csv(CCstepbysweep, file = paste0(abffile, ".sweeps.csv"))
  return(summaryplot)
}