#' Process an abf file of the slow steps VC protocol
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
process.VCstep = 
  function(abffile, listobj = F){
    VCdata = read.multisweep.pyth(abffile)
    
    # Vcmd is offset by 80 (Vjunc offset within amp)
    VCdata$Vm = VCdata$Vcmd0 - 80
    #this is where the jump happens
    tjump = 0.06562
    
    # All the data of the timecourse 
    plot0 = VCdata %>%
      filter(.data$t > tjump, .data$t < 0.75) %>%
      ggplot(aes(x = .data$t, y = .data$Imemb)) +
      geom_line(aes(group = .data$sweep)) +
      theme(axis.title.y = element_blank(),
            legend.position = "none")
    
    # collect waveform measurements
    VCdatabysweep = 
      VCdata %>%
      group_by(.data$sweep) %>%
      summarise(V = wctools::measure.timepoint(var = .data$Vm, tvar = .data$t, tp = 0.3, win = 0.1),
                
                Ijump1 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp = tjump + 0.001, win = 0.0002),
                Ijump5 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.005, win = 0.001),
                Ijump10 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.01, win = 0.001),
                Ijump20 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.02, win = 0.001),
                Ijump40 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.04, win = 0.001),
                Ijump80 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.08, win = 0.001),
                Ijump160 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.16, win = 0.001),
                Ijump320 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp =  tjump + 0.16, win = 0.001),
                Ijump550 = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp = tjump + 0.55, win = 0.05),
                Itail = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp = 0.6661, win = 0.0002),
                Itailss = wctools::measure.timepoint(var = .data$Imemb, tvar = .data$t, tp = 0.8661, win = 0.005))
    
    # calculate nbormalised tails
    VCdatabysweep = 
    VCdatabysweep %>%
      mutate(Itailsubs = .data$Itail - mean(VCdatabysweep$Itailss),
             normtail = .data$Itailsubs / min(.data$Itailsubs))
    
    # measurements of the peak 
    peakdata =
      VCdata %>%
      filter(.data$t < tjump + 0.02, .data$t > tjump + 0.0002) %>%
      group_by(.data$sweep) %>%
      summarise(
        Ipeak.in = min(.data$Imemb),
        tpeak.in = mean(t[.data$Imemb == min(.data$Imemb)]) - tjump,
        Ipeak.out = max(.data$Imemb),
        tpeak.out = mean(t[.data$Imemb == max(.data$Imemb)]) - tjump
      )
    
    # merge it all
    VCdatabysweep =
      left_join(VCdatabysweep, peakdata)
    
    # write to file in case TC analysis fails  
   readr::write_csv(VCdatabysweep, file = paste0(abffile, "VCmeasurements.csv"))
    
    # Boltzmann fit
    fit = minpack.lm::nlsLM(
      formula = normtail ~ A / (1 + exp((V - Vh) / k)),
      start = list(Vh = -50, k = -20, A = 1.1),
      data = VCdatabysweep
    )
    tidyfit = broom::tidy(fit)
    
    # Time course fitting
    
    # Create dataset to fit on
    VCfitdata =
      VCdata %>%
      group_by(.data$sweep, .data$Vm) %>%
      # filter out the time after jumps
      filter(.data$t > tjump + 0.02, .data$t < tjump + 0.60) %>%
      mutate(
        Isubs = .data$Imemb - wctools::measure.timepoint(
          var = .data$Imemb,
          tvar = .data$t,
          tp = tjump + 0.55,
          win = 0.05
        ),
        Inorm = .data$Isubs / max(.data$Isubs),
        tfit = .data$t - tjump
      ) %>%
      select(-.data$Isubs) %>%
      #get rid of voltages without activation
      filter(.data$Vm > -20, .data$Inorm < 0.95, .data$Inorm > 0.03) 
    
   # Fit a biexponential to everything
    safefit = purrr::safely(minpack.lm::nlsLM)
    
    KAfits =
      VCfitdata %>%
      group_by(.data$Vm, .data$sweep) %>%
      nest() %>%
      mutate(fit = purrr::map(.data$data, ~ safefit(
        formula = Inorm ~ SSasymp(tfit, Asym, R0, lrc),
        data = .
      )),
      fit1 = purrr::map(.data$data, ~ safefit(
        formula = Imemb ~ SSasymp(tfit, Asym, R0, lrc),
        data = .
      )))

    # unpack the results
    KAfits$fit = purrr::transpose(KAfits$fit1)$result

    KAfits =
      KAfits %>%
      mutate(
        params = purrr::map(.data$fit, broom::tidy),
        KAcurves = purrr::map(.data$fit, generics::augment, newdata = tibble(tfit = seq(0.02, 0.55, 0.01)))
      )
    
   # option for offline-use of the analysis (plotting)
    if (listobj == T) {
      return(list(VCdata,
                  VCfitdata,
                  KAfits2,
                  VCdatabysweep))
      }
   # Plot of the raw tc fits
    TCfit =
      VCfitdata %>%
      ggplot(aes(x = .data$tfit, y = .data$Imemb)) +
      # geom_hline(yintercept = c(0, 1),
      #            linetype = 3,
      #            colour = "grey50") +
      geom_line(size = 0.1, aes(group = sweep)) +
      facet_wrap( ~ .data$Vm, ncol = 3, scales = "free_y") +
      geom_smooth(method = "nls",
                  aes(group = sweep),
                  formula = y ~ SSasymp(x,Asym, R0, lrc),
                  se = F,
                  linetype = 3,
                  colour = "red") 
    
    KAparams =
      KAfits %>%
      group_by(.data$sweep, .data$Vm) %>%
      select(.data$params) %>%
      unnest(.data$params) %>%
      select(.data$term, .data$estimate) %>%
      pivot_wider(names_from = .data$term,
                  values_from = .data$estimate,
                  names_prefix = "IKA.") %>%
      # rescale from log rate constant to time constant in ms
      mutate(IKA.tau1 = (1/exp(.data$IKA.lrc)) * 1000)
    
    IKAparamplot =
      KAparams %>%
      group_by(.data$sweep, .data$Vm) %>%
      dplyr::select(-IKA.lrc) |> 
      tidyr::pivot_longer(cols = tidyselect::starts_with("IKA")) |> 
      ggplot(aes(x = .data$Vm,
                 y = .data$value)) +
      geom_hline(yintercept = 0,
                 linetype = 3,
                 colour = "grey50") +
      geom_point() +
      facet_wrap(~name, scales = "free_y") +
      geom_line(linetype = 3) +
      theme(axis.title.y = element_blank(),
      legend.title = element_blank())
    
    #IKAparamplot
    
    
    # Add KA fits to sweep summary
    VCdatabysweep =
      VCdatabysweep |>
      left_join(KAparams |> dplyr::select(.data$IKA.tau1))
 
     
    # create summary plots
    # IV
    plot1 =
      VCdatabysweep %>%
      pivot_longer(., cols = starts_with("Ijump")) %>%
      mutate(tjump = readr::parse_number(.data$name)) %>%
      ggplot(., aes(x = .data$V, y = .data$value, colour = .data$tjump)) +
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
        aes(group = tjump)
      ) +
      geom_point() +
      labs(x = "Vm, mV",
           y = "Im, pA") +
      scale_colour_viridis_c(
        option = "B",
        trans = "log10",
        guide = guide_colorbar(),
        name = "Time after pulse,\n ms"
      ) +
      theme(#legend.title = element_blank(),
        panel.background = element_rect(fill = "grey90"))

    # tails
    plot2 = ggplot(VCdatabysweep, aes(x = .data$V, y = .data$normtail)) +
      labs(y = "I tail, norm.",
           x = "Vm, mV") +
      geom_hline(yintercept = c(0, c(tidyfit$estimate[3]), c(tidyfit$estimate[3]/2)),
                 colour = "grey50",
                 linetype = 3) +
      geom_point() +
      annotate(geom = "text",label = paste0("Vh =", plyr::round_any(tidyfit$estimate[1],0.1)), x = 0, y = 0.5, colour = "deeppink") +
      geom_line(data = generics::augment(fit), aes(y = .data$.fitted), colour = "deeppink")


    # peaks
    plot3 =
      VCdata %>%
      filter(.data$t < tjump + 0.03, .data$t > tjump + 0.0002) %>%
      ggplot(aes(x = .data$t - tjump, y = .data$Imemb)) +
      geom_line(aes(group = .data$sweep)) +
      geom_point(data = VCdatabysweep, aes(y = .data$Ipeak.in, x = .data$tpeak.in ), colour = "red") +
      geom_point(data = VCdatabysweep, aes(y = .data$Ipeak.out, x = .data$tpeak.out ), colour = "blue")

    plot4 =
      ggplot(VCdatabysweep, aes(x = .data$V, y = .data$Ipeak.in)) +
      geom_line(colour = "red", alpha = 0.2) +
      labs(y = "I Peak Inw., pA",
           x = "Vm, mV") +
      scale_y_continuous(limits = c(NA, 0)) +
      geom_hline(linetype = 3,
                 colour = "grey50",
                 yintercept = 0) +
      geom_vline(linetype = 3,
                 colour = "grey50",
                 xintercept = 0) +
      geom_point(colour = "red")


    plot5 =
      ggplot(VCdatabysweep, aes(x = .data$V, y = .data$Ipeak.out)) +
      geom_line(colour = "blue", alpha = 0.2) +
      labs(y = "I Peak outw., pA",
           x = "Vm, mV") +
      scale_y_continuous(limits = c(0, NA)) +
      geom_hline(linetype = 3,
                 colour = "grey50",
                 yintercept = 0) +
      geom_vline(linetype = 3,
                 colour = "grey50",
                 xintercept = 0) +
      geom_point(colour = "blue")

    # delay of peak
    plot6 =
      VCdatabysweep %>%
      select(.data$V, .data$tpeak.in, .data$tpeak.out) %>%
      pivot_longer(starts_with("t"), values_to = "tpeak") %>%
      rowwise() %>%
      mutate(name = stringr::str_split(pattern = "\\.", .data$name),
             Direction = .data$name[2]) %>%
      select(-.data$name) %>%
      filter((.data$Direction != "out" | .data$V > -70)) %>%
      ggplot(aes(x = .data$V, y = .data$tpeak * 1000, colour = .data$Direction)) +
      labs(x = "Vm, mV",
           y = "Time to Peak, ms") +
      geom_point() +
      scale_colour_manual(values = c("red", "blue")) +
      geom_line(alpha = 0.2)
    # 
    #make the layout
    (bigplot =
        (plot0 + plot1) / (TCfit + IKAparamplot) /(plot2 + plot6) + patchwork::plot_annotation(title = "VC Slow K currents",
                                                                                            subtitle = paste0(abffile, "\n Analysed on ", Sys.Date())))
    # write to file (overwrite)
    readr::write_csv(VCdatabysweep, file = paste0(abffile, "VCmeasurements.csv"))
    readr::write_csv(KAparams, file = paste0(abffile, "fit.csv"))
    
    # 
     ggsave(bigplot, width = 11, height = 11, filename = paste0(abffile, ".VCsteps.png"))
     return(bigplot)
  }