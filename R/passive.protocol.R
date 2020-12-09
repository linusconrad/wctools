# This contains the function to wrangle the "passive" protocol data
#' Process an abf file of the "passive" protocol
#' 
#' This function takes an .abf file of the passive protocol and returns to the folder a summary graph and a .csv.
#' These contain per-sweep summaries.
#' SR of 50000 is assumed.
#' Function is tailored to this protocol only.
#' @param abffile The File
#' @return a .png and a .csv with the summarised analysis. Named with the filename trunk of 'abffile'.
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @importFrom rlang .data
#' @export
process.passive = function(abffile) {
  data = read.multisweep(abffile, 50000)
  names(data) = c("V", "I", "t", "sweep")
  #
  #######
  # Calculating the simple parameters
  
  databysweep =
    data %>%
    group_by(sweep) %>%
    dplyr::summarise(
      #tscaled = t - 0.0299,
      Istim = measure.timepoint(.data$I, .data$t, 0.5, 0.02),
      Vsteady = measure.timepoint(.data$V, .data$t, 0.5, 0.02),
      RMP = measure.timepoint(.data$V, .data$t, 0.02, 0.02),
      Vmin = min(.data$V),
      #Vscaled = 1 - (V - RMP) / (Vmin - RMP),
      tpeak = min(.data$t[.data$V == .data$Vmin]),
      # the first point at which V is = Vmin, it mind find more for the garbage ones...
      sagcoef = .data$Vsteady / .data$Vmin
    )
  
  # calculate a linear regression for the input resistance
  regression =
    stats::lm(Vsteady ~ Istim, data = databysweep)
  
  coefs = generics::tidy(regression)
  R.squared = generics::glance(regression)$r.squared
  
  # make some summary graphs
  # raw timecourses
  trace = 
    ggplot(data, aes(x = .data$t, y = .data$V)) +
    labs(x = "t, s", y = "V, mV") +
    scale_x_continuous(expand = c(0, 0)) +
    geom_hline(linetype = 3,
               colour = "grey50",
               yintercept = 0) +
    # vertical coursors to show measurements taken
    geom_line(aes(group = .data$sweep))  +
    geom_vline(
      xintercept = c(0.5, 0.02),
      colour = "blue",
      linetype = 3
    )
  
    # summary graph with linear fit of input resistance
    IV =
      ggplot(databysweep, aes(x = .data$Istim, y = .data$Vsteady)) +
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

    # summary of the sag coefficient
    sag =
      ggplot(databysweep, aes(x = .data$Istim, y = .data$sagcoef)) +
      labs(x = latex2exp::TeX('$I_{stim}$, pA$'),
           y = latex2exp::TeX('$\\frac{V{ss}}{V_{peak}}$')) +
      #scale_x_continuous(expand = c(0, 0)) +
      geom_hline(linetype = 3,
                 colour = "grey50",
                 yintercept = 1) +
      geom_smooth() +
      geom_point()
    
  #######
  # tau fitting
  # make truncated normalised dataset to do the fits

  fitset =
        data %>%
        group_by(.data$sweep) %>%
        # repeat the extraction of the times to get the rising parts
        # this is not very efficient but its a fast calculation
        mutate(
          .,
          tscaled = .data$t - 0.0299,
          Istim = measure.timepoint(.data$I, .data$t, 0.5, 0.02),
          Vsteady = measure.timepoint(.data$V, .data$t, 0.5, 0.02),
          RMP = measure.timepoint(.data$V, .data$t, 0.02, 0.02),
          Vmin = min(.data$V),
          Vscaled = 1 - (.data$V - .data$RMP) / (.data$Vsteady - .data$RMP),
          tpeak = min(t[.data$V == .data$Vmin]),
          # the first point at which V is = Vmin, it mind find more for the garbage ones...
          sagcoef = .data$Vsteady / .data$Vmin
        ) %>%
       # for the fits do not take the first sweeps with sag.
       # use the first 20 ms of sweeps without
      filter(.data$sweep > 2,
             .data$t <  0.05)
    
    
    # run the fitting in failsafe mode to generate table output
    safe_nls = purrr::safely(stats::nls)
    #now do the fitting with a purr call and process the results
    fitset %<>%
      group_by(.data$sweep) %>%
      nest() %>%
      mutate(fit = purrr::map(.data$data, ~ safe_nls(
        formula = Vscaled ~ exp(-tscaled / tau),
        data = .,
        start = list(tau = 0.5)
      )))
    # # 
    # 
    # # Extract the results list from the safely output and overwrite the complicated nested column
    # fitset$fit = purrr::transpose(fitset$fit)$result
    # # 
    # # get all the coefs and tidy stuff
    # fitset %<>%
    #   mutate(.,
    #          coefs = purrr::map(.data$fit, generics::tidy),
    #          fitvalue = purrr::map(.data$fit, generics::augment))
    # 
    # # Make a plot of the fit
    (fits =
      #unnest(c(.data$data), keep_empty = T) %>%
      ggplot(fitset, aes(
        x = .data$t,
        y = .data$Vscaled,
        colour = as.factor(.data$sweep)
      )) +
      labs(x = "t, s",
           y = "V, normalised",
           title = "Time course fitting, each sweep") +
      facet_wrap(
        ~ .data$sweep,
        as.table = T,
        scales = "free",
        strip.position = "top",
      ) +
      #scale_y_continuous(limits = c(0, 1)) +
      geom_line(aes(group = .data$sweep)) +
      theme(legend.position = "none") 
      # geom_line(
      #   data = unnest(fitset, .data$fitvalue),
      #   aes(y = .data$.fitted),
      #   linetype = 3,
      #   colour = "grey50"
      # ) +
      # geom_text(
      #   data = unnest(fitset, coefs) ,
      #   aes(label = paste0("tau = ", plyr::round_any(.data$estimate * 1000, 1), " ms")),
      #   x = 0,
      #   y = 1,
      #   colour = "black",
      #   vjust = "inward",
      #   hjust = "inward"
      # )
      )
}


# 

# 


# 

# 
#   # Extract the results list from the safely output and overwrite the complicated nested column
#   fitset$fit = purrr:transpose(fitset$fit)$result
# 
#   # get all the coefs and tidy stuff
#   fitset %<>%
#     mutate(.,
#            coefs = purrr::map(fit, ~ tidy(.)),
#            fitvalue = purrr::map(fit, ~ augment(.)))
# 
#   # Make a plot of the fit
#   fits =
#     fitset %>%
#     unnest(., c(data), keep_empty = T) %>%
#     ggplot(., aes(
#       x = tscaled,
#       y = Vscaled,
#       colour = as.factor(sweep)
#     )) +
#     labs(x = "t, s",
#          y = "V, normalised",
#          title = "Time course fitting, each sweep")+
#     facet_wrap(
#       ~ sweep,
#       as.table = T,
#       scales = "free_x",
#       strip.position = "top",
#     ) +
#     scale_y_continuous(limits = c(0, 1)) +
#     geom_line(aes(group = sweep)) +
#     theme(legend.position = "none") +
#     geom_line(
#       data = unnest(fitset, fitvalue),
#       aes(y = .fitted),
#       linetype = 3,
#       colour = "grey50"
#     ) +
#     geom_text(
#       data = unnest(fitset, coefs) ,
#       aes(label = paste0("τ = ", round_any(estimate * 1000, 1), " ms")),
#       x = 0,
#       y = 1,
#       colour = "black",
#       vjust = "inward",
#       hjust = "inward"
#     )
# 
#   #############################
#   # write all to file
#   comb =   ggpubr::ggarrange(trace, IV, sag, ncol =3)
# 
#   cowplot::plot_grid(comb,
#                      fits,
#                      ncol = 1,
#                      align = "v",
#                      axis = "lr",
#                      rel_heights = c(1,4)) %>%
#     # add Text annotations
#     ggpubr::annotate_figure(
#       .,
#       top = text_grob(
#         "Hyperpolarising step protocol",
#         size = 14,
#         family  = "Palatino"
#       ),
#       bottom = text_grob(
#         paste0(
#           "Data source: \n",
#           abffile,
#           "\n RMP: ",
#           round_any(coefs$estimate[1], accuracy = 0.1) ,
#           " mV \n",
#           "Input Resistance: ",
#           round_any(coefs$estimate[2], accuracy = 0.01) * 1000 ,
#           " MOhm \n"
#         ),
#         family = "Palatino"
#       )
#     ) %>%
#     ggsave(
#       filename = paste0(abffile, ".passive.png"),
#       width = 10,
#       height = 15
#     )
# 
#   # extract the fit parameters and write all the summary to file
#   left_join(
#     databysweep,
#     fitset %>%
#       unnest(., coefs) %>%
#       select(., sweep, estimate) %>%
#       rename(., tau = estimate) %>%
#       # rescale to ms
#       mutate(., tau = tau * 1000)
#   ) %>%
#     write_csv(., path = paste0(abffile, ".sweeps.csv"))
# }