#' Process an abf file of the fast steps VC protocol (availability)
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

process.VCfavailability =
  function(abffile) {
    VCfavailability = read.multisweep.pyth(abffile)
    
    # Vcmd is offset by 80 (Vjunc offset within amp)
    VCfavailability$Vm = VCfavailability$Vcmd0 - 80
    
    
    # Make the waveform measurements
    tjump = 0.2039
    VCfbysweep = 
      VCfavailability %>%
      group_by(.data$sweep) %>%
      summarise(V = wctools::measure.timepoint(var = .data$Vm, tvar = .data$t, tp = 0.1, win = 0.1),
                Ijump1 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp = tjump + 0.001, win = 0.0002),
                Ijump5 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.005, win = 0.001),
                Ijump10 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.01, win = 0.001),
                Ijump15 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.015, win = 0.001),
                Ijump20 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.02, win = 0.001),
                Ijump25 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.025, win = 0.001),
                Ijump30 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.03, win = 0.001),
                Ijump35 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.35, win = 0.001))
    
    IK = mean(VCfbysweep$Ijump30)
    # measurements of the peak 
    VCfpeakdata =
      VCfavailability %>%
      filter(.data$t < tjump + 0.02, .data$t > tjump + 0.0002) %>%
      mutate(Isubs = .data$Imemb - IK) %>%
      group_by(.data$sweep) %>%
      summarise(
        Ipeak.in = min(.data$Imemb),
        tpeak.in = t[.data$Imemb == min(.data$Imemb)] - tjump,
        INaT = min(.data$Isubs)) %>%
      ungroup() %>%
      mutate(normINaT = .data$INaT / min(.data$INaT))
    
    VCfbysweep =
      left_join(VCfbysweep, VCfpeakdata)
    
    
    # Make a Boltzman model
    fit = minpack.lm::nlsLM(
      formula = normINaT ~ A + 1 / (1 + exp((V - Vh) / k)),
      start = list(Vh = -50, k = 20, A = 0.1),
      data = VCfbysweep
    )
    
    tidyfit = broom::tidy(fit)
    
    plottc = VCfavailability %>%
      filter(.data$t >tjump, .data$t < 0.215) %>%
      ggplot(aes(x = .data$t, y = .data$Imemb - IK)) +
      scale_y_continuous(limits = c(NA, 100)) +
      geom_line(aes(group = .data$sweep), size = 0.2) +
      ggthemes::scale_colour_solarized(accent = "blue") +
      theme(axis.title.y = element_blank(),
            legend.position = "none")
    
    
    plot0 =
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
      labs(x = "Prepulse, mV",
           y = "I NaT @ -20 mV, pA") +
      scale_colour_viridis_c(
        option = "B",
        trans = "log10",
        guide = guide_colorbar(),
        name = "Time after pulse,\n ms"
      ) +
      theme(#legend.title = element_blank(),
        panel.background = element_rect(fill = "grey90"))
    
    plot1 = VCfbysweep %>%
      ggplot(aes(x = .data$V, y = .data$normINaT)) +
      geom_hline(yintercept = c(1, c(tidyfit$estimate[3]), c((1+tidyfit$estimate[3])/2)),
                 colour = "grey50",
                 linetype = 3) +
      labs(y = "Normalised Inw. Peak",
           x = "Pre pulse, mV") +
      geom_line(data = generics::augment(fit, newdata = data.frame(V = c(seq(-120,0)))), aes(y = .data$.fitted), colour = "deeppink") +
      geom_point(shape = 22) +
      annotate(geom = "text", label = c(paste0("Vh = ", plyr::round_any(broom::tidy(fit)$estimate[1],0.1), " mV")), x =-40, y = 0.7, colour = "deeppink") 
    
    plot2 = VCfbysweep %>%
      ggplot(aes(x = .data$V, y = .data$tpeak.in * 1000)) +
      # geom_hline(yintercept = c(0, 1),
      #            colour = "grey50",
      #            linetype = 3) +
      labs(y = "Time to peak, ms",
           x = "Pre pulse, mV") +
      geom_point()
    combined =
      (plottc +  plot0)/ (plot1 + plot2) + patchwork::plot_annotation(title = "VC Fast Availability",
                                                                      subtitle = paste0(abffile, "\n Analysed on ", Sys.Date()))
    #plottc + plot1
    # write to file
    readr::write_csv(VCfbysweep, file = paste0(abffile, "VCfavailability.csv"))
    ggsave(combined, width = 8.25, height = 5.875, filename = paste0(abffile, ".VCfavail.png"))
    return(combined)
  }