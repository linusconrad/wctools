# This contains the function to wrangle the pulsetrain protocol data (IO curves)
#' Process an abf file of the "Input - Output Curve" protocol
#' 
#' This function takes an .abf file of the protocol and returns to the folder a summary graph and a .csv.
#' These contain per-sweep summaries.
#' SR of 50000 is assumed. The protocol has 100 hz stimulations.
#' Function is tailored to this protocol only.
#' @param abffile The File
#' @param Vjunc Junction potential to add
#' @param threshold threshold for the AP detector algorithm
#' @return a .png and a .csv with the summarised analysis. Named with the filename trunk of 'abffile'.
#' @import ggplot2
#' @import tidyr
#' @import patchwork
#' @import dplyr
#' @importFrom rlang .data
#' @export
process.100hz =
  function(abffile, Vjunc, threshold) {
    #default to 0 mV AP detection
    if(missing(threshold))
      threshold <- 0
    
    stimtrain = read.multisweep.pyth(file = abffile)
    
    stimtrain$Vmemb = stimtrain$Vmemb + Vjunc
    
    # Stimulation rate and sampling in Hz
    SR = 50000
    f = 100
    # length of one cycle  in samples
    cyclength = 1 / f * SR
    # number of stimulation cycles (half a second)
    nstim = f * 0.5
    
    stimtrain %<>%
      #remove areas with no stimulation (RMP)
      filter(.data$t > 0.0593, .data$t < 0.5593) %>%
      group_by(.data$sweep) %>%
      # add cycle time
      mutate(cycle = rep(c(1:nstim), each = cyclength),
             #add dV variable
             dV = c(NA, (diff(.data$Vmemb) / 1000 / (1 / 50000)))) %>%
      group_by(.data$cycle, .data$sweep) %>%
      # add the timescale by cycle
      mutate(tcycle = seq(
        from = 0,
        by = 1 / SR,
        length.out = cyclength
      ) * 1000)
    
    # add a stimulation variable
    stimtrain =
      left_join(stimtrain,
                stimtrain %>%
                  ungroup() %>%
                  group_by(.data$sweep) %>%
                  summarise(Istim = max(.data$ICmd1nA)))
    
    #save some space...
    stimtrain$Istim = as.factor(stimtrain$Istim)
    
    Apsignatures =
      stimtrain %>%
      ggplot(aes(x = .data$tcycle, y = .data$Vmemb, colour = .data$cycle)) +
      geom_line(aes(group = .data$cycle), alpha = 0.3) +
      labs(title = "Raw data aligned to stimulus cycle") +
      geom_hline(yintercept = 0,
                 linetype = 3,
                 colour = "grey50") +
      facet_wrap( ~ .data$sweep, nrow = 1, as.table = T) +
      cowplot::theme_nothing() +
      theme(title = element_text())
    
    ApsignaturesdV =
      stimtrain %>%
      ggplot(aes(x = .data$Vmemb, y = .data$dV, colour = .data$cycle)) +
      geom_vline(xintercept = 0,
                 linetype = 3,
                 colour = "grey50") +
      geom_path(aes(group = .data$cycle), alpha = 0.3) +
      labs(y = latex2exp::TeX("$\\frac{d\\,V_m}{d\\,t}$, V/s"),
           x = latex2exp::TeX("$V_m$, mV")) +
      geom_hline(yintercept = 0,
                 linetype = 3,
                 colour = "grey50") +
      facet_wrap(
        ~ .data$Istim,
        nrow = 2,
        as.table = T,
        labeller = labeller(
          Istim = function(x)
            paste0(x, ", pA")
        )
      )
    #  theme(title = element_text())
    
    #Extract all the AP
    stimtrainAP =
      stimtrain %>%
      ungroup() %>%
      wctools::getAPstats("Vmemb", thresh3 = threshold) %>%
      rename(latency = .data$tcycle)
    
    # Plot the AP params
    APparamplot =
      stimtrainAP %>%
      mutate(
        Istim = as.numeric(as.character(.data$Istim)),
        tmin = .data$tmin * 1000,
        width = .data$width * 1000
      ) %>% # rescales these variables to mV
      group_by(.data$sweep, .data$Istim, .data$cycle) %>%
      select(.data$Vpeak, .data$upstroke, .data$latency, .data$Vmin, .data$tmin, .data$width) %>%
      pivot_longer(cols = c(.data$Vpeak, .data$upstroke, .data$latency, .data$Vmin, .data$tmin, .data$width)) %>%
      ggplot(aes(
        x = .data$cycle,
        y = .data$value,
        colour = .data$Istim / 1000
      )) +
      labs(x = "Stimulation Cycle",
           title = "Spike properties") +
      scale_colour_viridis_c(name = latex2exp::TeX("$I_{stim}$, nA"), option = "A") +
      geom_hline(yintercept = 0,
                 linetype = 3,
                 colour = "grey50") +
      geom_line(aes(group = sweep)) +
      facet_wrap( ~ name, strip.position = "left", scales = "free") +
      theme(
        axis.title.y = element_blank(),
        legend.position = "right",
        panel.background = element_rect(fill = "grey70")
      )
    
    # plot input output curve
    # prepare the data
    IOdata =
      left_join(
        stimtrain %>%
          ungroup() %>%
          group_by(.data$sweep) %>%
          summarise(Istim = max(.data$ICmd1nA)),
        stimtrainAP %>%
          group_by(.data$sweep) %>%
          summarise(nAP = max(.data$APindex))
      )
    
    IOdata %<>%
      replace_na(replace = list(nAP = 0)) %>%
      mutate(firedfraction = (.data$nAP / 50) * 100) # part of stimuli that elicited an AP in %
    
    
    Iocurve = ggplot(IOdata, aes(x = .data$Istim, y = .data$firedfraction)) +
      labs(x = latex2exp::TeX("$I_{stim}$, pA"),
           y = "AP elicited, %") +
      geom_step() +
      geom_hline(
        yintercept = c(0, 50, 100),
        linetype = 3,
        colour = "grey50"
      )
    
    # make a massive plot of all combined
    (combined = (Apsignatures/ ApsignaturesdV + theme(legend.position = "right"))/(APparamplot + Iocurve + patchwork::plot_layout(widths = c(1,0.5))))
    ggsave(combined, file = paste0(abffile, "IO-plot.png"),  width = 11, height = 11)
    utils::write.csv(stimtrainAP, file = paste0(abffile, "100hzAPstats.csv"))
    return(combined)
    
  }