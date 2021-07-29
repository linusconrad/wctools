#' Process an abf file of the fast steps VC protocol
#' 
#' This function takes an .abf file of the protocol and returns to the folder a summary graph and a .csv.
#' These contain per-sweep summaries.
#' SR of 50000 is assumed.
#' In voltage clamp I substract the junction as I record with an offset.
#' Stimulus V in the file is relative to -80 holding!
#' @param abffile The File
#' @return a .png and a .csv with the summarised analysis. Named with the filename trunk of 'abffile'.
#' @import ggplot2
#' @import tidyr
#' @import patchwork
#' @import dplyr
#' @importFrom rlang .data
#' @export
process.VCf =
  function(abffile) {
    VCf = read.multisweep.pyth(file = abffile)
    # Make the waveform measurements
    tjump = 0.05390
    
    #substract holding
    VCf$Vcmd0 =  VCf$Vcmd0 - 80
    
    VCfbysweep = 
      VCf %>%
      group_by(.data$sweep) %>%
      summarise(V = wctools::measure.timepoint(var = .data$Vcmd0, tvar = .data$t, tp = 0.07, win = 0.01),
                Ijump1 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp = tjump + 0.001, win = 0.0002),
                Ijump2 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.002, win = 0.001),
                Ijump3 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.003, win = 0.001),
                Ijump5 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.005, win = 0.001),
                Ijump7 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.007, win = 0.001),
                Ijump10 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.01, win = 0.001),
                Ijump15 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.02, win = 0.001),
                Ijump20 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.025, win = 0.001),
                Ijump30 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.03, win = 0.001),
                Ijump35 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.035, win = 0.001),
                Ijump40 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.04, win = 0.001),
                Ijump45 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.045, win = 0.001),
                Ijump49 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.049, win = 0.001),
                Itail   = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp = 0.1049, win = 0.0002),
                Itailss = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp = 0.2, win = 0.005))
    
    
    VCfbysweep =
      VCfbysweep %>%
      mutate(tailsubs = .data$Itail - mean(VCfbysweep$Itailss),
             normtail = .data$tailsubs/min(.data$tailsubs)) 
    
    #  Fit a Boltzmann for K activation
    # Make a Boltzman model
    # Added a fudge factor to account for imperfect saturation (A)
    fit = minpack.lm::nlsLM(
      formula = normtail ~ A / (1 + exp((V - Vh) / k)),
      start = list(Vh = -50, k = -20, A = 1.1),
      data = VCfbysweep
    )
    
    tidyfit = broom::tidy(fit)
    
    # waveform measurements of the NApeak
    peakdata =
      VCf %>%
      filter(.data$t < tjump + 0.03, .data$t > tjump + 0.0005) %>%
      group_by(.data$sweep) %>%
      summarise(
        Ipeak.in = min(.data$Imemb),
        tpeak.in = t[.data$Imemb == min(.data$Imemb)][1] - tjump)
    
    peakdata =
      left_join(peakdata,
                VCf %>%
      filter(.data$t < tjump + 0.03, .data$t > tjump + 0.0015) %>%
      group_by(.data$sweep) %>%
      summarise(
        Ipeak.out = max(.data$Imemb),
        tpeak.out = t[.data$Imemb == max(.data$Imemb)][1] - tjump
      ) )
      
    
    VCfbysweep =
      left_join(VCfbysweep, peakdata)
    
    # Plot tc, and IV
    
    VCftc =
      VCf %>%
      filter(.data$t > 0.04, .data$t < 0.07) %>%
      ggplot(aes(x = t, y = .data$Imemb)) +
      geom_line(aes(group = .data$sweep), size = 0.2) +
      theme(axis.title.y = element_blank(),
            legend.position = "none")
    
    
    VCfIV =
      VCfbysweep %>%
      pivot_longer(cols = starts_with("Ijump")) %>%
      mutate(tjump = readr::parse_number(.data$name)) %>%
      ggplot(aes(x = .data$V, y = .data$value, colour = .data$tjump)) +
      geom_hline(linetype = 3,
                 colour = "grey50",
                 yintercept = 0) +
      geom_vline(linetype = 3,
                 colour = "grey50",
                 xintercept = 0) +
      geom_line(
        linetype = 1,
        #se = F,
        #geom = "line",
        alpha = 0.2,
        size = 0.5,
        aes(group = .data$tjump)
      ) +
      geom_point() +
      labs(x = "Vm, mV",
           y = "I, pA") +
      scale_colour_viridis_c(
        option = "B",
        trans = "log10",
        guide = guide_colorbar(),
        name = "Time after pulse,\n ms"
      ) +
      theme(#legend.title = element_blank(),
        panel.background = element_rect(fill = "grey90"))
    
    
    plottail= VCfbysweep %>%
      ggplot(aes(x = .data$V, y = .data$normtail)) +
      geom_hline(yintercept = c(0, c(tidyfit$estimate[3]), c(tidyfit$estimate[3]/2)),
                 colour = "grey50",
                 linetype = 3) +
      labs(y = "Normalised Tail Current",
           x = "Vm, mV") +
      geom_point() +
      annotate(geom = "text",label = paste0("Vh =", plyr::round_any(tidyfit$estimate[1],0.1)), x = 0, y = 0.5, colour = "deeppink") +
      geom_line(data = generics::augment(fit), aes(y = .data$.fitted), colour = "deeppink")
    
    plot2 = VCfbysweep %>%
      filter(abs(.data$Ipeak.in) > 500) %>%
      ggplot(aes(x = .data$V, y = .data$tpeak.in * 1000)) +
      # geom_hline(yintercept = c(0, 1),
      #            colour = "grey50",
      #            linetype = 3) +
      labs(y = "Time to In peak, ms",
           x = "Vm, mV") +
      geom_point() +
      geom_smooth()
    
    (combinedplot = (VCftc + VCfIV) / (plottail + plot2) + patchwork::plot_annotation(title = "VC Fast Steps",
                                                                                   subtitle = paste0(abffile, "\n Analysed on ", Sys.Date())))
    
    # write to file
    readr::write_csv(VCfbysweep, file = paste0(abffile, "VCf.csv"))
    ggsave(combinedplot, width = 8.25, height = 5.875, filename = paste0(abffile, ".VCfsteps.png"))
    return(combinedplot)
  }