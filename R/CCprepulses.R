# This contains the function to wrangle the CC protocol with prepulse
#' Process an abf file of the CC prepulse (bifurcation test) protocol
#' 
#' This function takes an .abf file of the passive protocol and returns to the folder a summary graph and a .csv.
#' These contain per-sweep summaries.
#' SR of 50000 is assumed.
#' Function is tailored to this protocol only.
#' @param abffile The File
#' @param Vjunc Junction potential to add
#' @param threshold Voltage threshold of the AP detector, defaults to 0
#' @return a .png and a .csv with the summarised analysis. Named with the filename trunk of 'abffile'.
#' @import ggplot2
#' @import tidyr
#' @import patchwork
#' @import dplyr
#' @importFrom rlang .data
#' @export

process.CCbif = function(abffile, Vjunc, threshold) {
  #default to 0 mV
  if(missing(threshold))
    threshold <- 0
  # assign time of stimulus onset
  tjump = 0.37936
  
  # get the data
  bifurcationdata =
    read.multisweep.pyth(abffile)
  
  # substract junction
  bifurcationdata = 
  bifurcationdata %>%
    mutate(Vmemb = .data$Vmemb + Vjunc)
  
  # get AP stats, calculate latency based on stimulation onset
  bifAP =
    bifurcationdata %>%
    filter(.data$t > tjump) %>%
    wctools::getAPstats("Vmemb", threshold) %>%
    mutate(
      latency = .data$tpeak - tjump,
      halfpoint1 = .data$width - .data$halfpoint,
      halfstart = .data$tpeak - .data$halfpoint1,
      halfend = .data$tpeak + .data$halfpoint
    ) %>%
    select(-.data$halfpoint,-.data$halfpoint1)
  
  
  # plot the trace
  biftrace =
    bifurcationdata %>%
    ggplot(aes(x = .data$t, y = .data$Vmemb)) +
    geom_hline(yintercept = c(0,threshold),
               colour = "grey50",
               linetype = 3) +
    facet_wrap(~ .data$sweep, ncol = 4) +
    geom_vline(
      data = bifAP,
      aes(xintercept = .data$tpeak),
      colour = "blue",
      alpha = 0.5
    ) +
    geom_line(size = 0.2) +
    cowplot::theme_nothing()
  
  # Do QC of the AP finder
  APQC = 
  bifurcationdata %>%
    filter(.data$t > tjump, .data$t < (tjump + 0.015)) %>%
    ggplot(aes(x = .data$t, y = .data$Vmemb)) +
    geom_hline(yintercept = c(0),
               colour = "grey50",
               linetype = 3) +
    geom_hline(yintercept = c(threshold),
               colour = "blue",
               linetype = 3) +
    facet_wrap( ~ sweep,
                ncol = 5,
                as.table = T) +
    geom_vline(
      data = bifAP,
      aes(xintercept = .data$tpeak),
      colour = "blue",
      alpha = 0.5
    ) +
    geom_line(size = 0.2) +
    geom_segment(
      data = bifAP,
      mapping = aes(
        x = .data$halfstart,
        y = .data$Vhalf,
        xend = .data$halfend,
        yend = .data$Vhalf
      ),
      colour = "red"
    ) +
    cowplot::theme_nothing()
  
  
  # Calculating the WF measurements
  databysweep =
    bifurcationdata %>%
    group_by(.data$sweep) %>%
    dplyr::summarise(
      Ipre = measure.timepoint(.data$Icmd, .data$t, 0.35, 0.02),
      Vpre = measure.timepoint(.data$Vmemb, .data$t, 0.35, 0.02),
      RMP = measure.timepoint(.data$Vmemb, .data$t, 0.03, 0.02),
      Vsag = min(.data$Vmemb),
      #Vscaled = 1 - (V - RMP) / (Vmin - RMP),
      tsag = min(.data$t[.data$Vmemb == .data$Vsag]),
      # the first point at which V is = Vmin, it mind find more for the garbage ones...
      sagcoef = .data$Vpre / .data$Vsag
    )
  
  RMP = mean(databysweep$RMP)
  # merge the WF measurements and AP params
  # calculate the number of AP
  databysweep =
    left_join(databysweep,
              bifAP %>%
                group_by(.data$sweep) %>%
                dplyr::summarise(nAP = max(.data$APindex)))
    
  databysweep =
    left_join(databysweep, bifAP %>%
                filter(.data$APindex == 1))
  
  # calculate a linear regression for the input resistance
  regression =
    stats::lm(Vpre ~ Ipre, data = databysweep)
  
  coefs = broom::tidy(regression)
  R.squared = generics::glance(regression)$r.squared
  
  # summary graph with linear fit of input resistance
  bifVI =
    ggplot(databysweep, aes(x = .data$Ipre, y = .data$Vpre)) +
    labs(x = latex2exp::TeX('$I_{stim}$, pA$'),
         y = latex2exp::TeX('$V_{ss}$, mV$')) +
    geom_hline(yintercept = 0,
               linetype = 3,
               colour = "grey50") +
    geom_vline(xintercept = 0,
               linetype = 3,
               colour = "grey50") +
    geom_point() +
    ggpubr::stat_regline_equation(colour = "red") +
    stat_smooth(method = "lm")
  
  # Plot the AP properties alongside the conditioning potential.
  bifsummary =
    databysweep %>%
    mutate(# rescales these variables to ms
      tmin = .data$tmin * 1000,
      width = .data$width * 1000,
      latency = .data$latency * 1000
    ) %>% 
    select(.data$sweep, .data$Vpre, .data$Vpeak, .data$Vmin, .data$tmin, .data$width, .data$upstroke, .data$latency) %>%
    pivot_longer(cols = c(.data$Vpeak, .data$Vmin, .data$tmin, .data$width, .data$upstroke, .data$latency)) %>%
    ggplot(aes(x = .data$Vpre, y = .data$value)) +
    facet_wrap( ~ name, scales = "free_y", strip.position = "left") +
    #geom_line(linetype =3, colour = "grey50")+
    geom_point() +
    geom_smooth(se = F, linetype = 2) +
    geom_hline(yintercept = 0,
               linetype = 3,
               colour = "grey50") +
    geom_vline(xintercept = RMP,
               linetype = 3,
               colour = "red") +
    theme(axis.title.y = element_blank())
  
  bifcombined = 
    (biftrace / APQC /bifsummary) + patchwork::plot_annotation(title = "CC Step w/ Prepulse",
                                                    subtitle = paste0(abffile, "\n Analysed on ", Sys.Date()))
  # Write all to file
  utils::write.csv(databysweep, file = paste0(abffile, "CCbifsummary.csv"))
  ggsave(bifcombined, file = paste0(abffile, "CCbifplot.png"),  width = 8.25, height = 11.75)
  return(bifcombined)
}