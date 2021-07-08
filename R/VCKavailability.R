#' Process an abf file of the KA availability VC protocol
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

process.VCKavailability =
  function(abffile) {
    VCKavailability = read.multisweep.pyth(abffile)
    
    # Vcmd is offset by 80 (Vjunc offset within amp)
    VCKavailability$Vm = VCKavailability$Vcmd0 - 80
    
    
    # Make the waveform measurements
    tjump = 0.67342
    VCKbysweep = 
      VCKavailability %>%
      group_by(.data$sweep) %>%
      summarise(V = wctools::measure.timepoint(var = .data$Vm, tvar = .data$t, tp = 0.5, win = 0.1),
                Ijump1 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp = tjump + 0.001, win = 0.0002),
                Ijump5 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.005, win = 0.001),
                Ijump7 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.007, win = 0.001),
                Ijump10 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.01, win = 0.001),
                Ijump20 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.02, win = 0.001),
                Ijump40 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.04, win = 0.001),
                Ijump80 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.08, win = 0.001),
                Ijump160 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.16, win = 0.001),
                Ijump320 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.16, win = 0.001),
                Ijump395 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp = tjump + 0.395, win = 0.05))
    
    Iss = min(VCKbysweep$Ijump395)
    # measurements of the peak 
    VCKpeakdata =
      VCKavailability %>%
      filter(.data$t < tjump + 0.02, .data$t > tjump + 0.0002) %>%
      group_by(.data$sweep) %>%
      summarise(
        Ipeak.out = max(.data$Imemb),
        tpeak.out = .data$t[.data$Imemb == max(.data$Imemb)][1] - tjump) %>%
      ungroup() 
    
    VCKbysweep =
      left_join(VCKbysweep, VCKpeakdata)%>%
      mutate(Isubs= .data$Ijump7 - Iss,
             normKA = .data$Isubs / max(.data$Isubs))
    
    
    #  Fit a Boltzmann for K activation
    # Make a Boltzman model
    # Added a fudge factor to account for imperfect saturation (A)
    fit = minpack.lm::nlsLM(
      formula = normKA ~ A + 1 / (1 + exp((V - Vh) / k)),
      start = list(Vh = -50, k = 20, A = 0.2),
      data = VCKbysweep
    )
    
    tidyfit = broom::tidy(fit)
    
    plottc = VCKavailability %>%
      filter(.data$t > 0.5, .data$t < 1) %>%
      ggplot(aes(x = .data$t, y = .data$Imemb)) +
      geom_line(aes(group = .data$sweep), size = 0.2) +
      ggthemes::scale_colour_solarized(accent = "blue") +
      theme(axis.title.y = element_blank(),
            legend.position = "none")
    
    
    plot0 =
      VCKbysweep %>%
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
           y = "K current @ +50 mV, pA") +
      scale_colour_viridis_c(
        option = "B",
        trans = "log10",
        guide = guide_colorbar(),
        name = "Time after pulse,\n ms"
      ) +
      theme(#legend.title = element_blank(),
        panel.background = element_rect(fill = "grey90"))
    
    plot1 = VCKbysweep %>%
      ggplot(aes(x = .data$V, y = .data$normKA)) +
      labs(y = "Normalised Peak IK",
           x = "Pre pulse, mV") +
      geom_hline(yintercept = c(1, c(tidyfit$estimate[3]), c((1+tidyfit$estimate[3])/2)),
                 colour = "grey50",
                 linetype = 3) +
      geom_point() +
      annotate(geom = "text",label = paste0("Vh =", plyr::round_any(tidyfit$estimate[1],0.1)), x = 0, y = 0.5, colour = "deeppink") +
      geom_line(data = generics::augment(fit), aes(y = .data$.fitted), colour = "deeppink")
    
    plot2 = VCKbysweep %>%
      ggplot(aes(x = .data$V, y = .data$tpeak.out * 1000)) +
      # geom_hline(yintercept = c(0, 1),
      #            colour = "grey50",
      #            linetype = 3) +
      labs(y = "Time to peak, ms",
           x = "Pre pulse, mV") +
      geom_point() +
      geom_smooth()
    
    combined = 
      (plottc +  plot0)/ (plot1 + plot2) + patchwork::plot_annotation(title = "VC IKA Availability",
                                                                      subtitle = paste0(abffile, "\n Analysed on ", Sys.Date()))
    
    # write to file
    readr::write_csv(VCKbysweep, file = paste0(abffile, "VCKavailability.csv"))
    ggsave(combined, width = 8.25, height = 11.75, filename = paste0(abffile, ".VCKsteps.png"))
    return(combined)
  }