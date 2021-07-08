# This contains the function to wrangle the "passive" protocol data
#' Process an abf file of the "passive" protocol
#' 
#' This function takes an .abf file of the passive protocol and returns to the folder a summary graph and a .csv.
#' These contain per-sweep summaries.
#' SR of 50000 is assumed.
#' Function is tailored to this protocol only.
#' @param abffile The File
#' @param Vjunc Junction potential to add
#' @return a .png and a .csv with the summarised analysis. Named with the filename trunk of 'abffile'.
#' @import ggplot2
#' @import tidyr
#' @import patchwork
#' @importFrom magrittr %<>%
#' @import dplyr
#' @importFrom rlang .data
#' @export
process.passive = function(abffile, Vjunc) {
  data = read.multisweep.pyth(abffile)
  names(data) = c("t", "V", "I", "sweep")
  #
  # substract the junction potential
  data %>% mutate(V = .data$V + Vjunc)
  
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
  
  coefs = broom::tidy(regression)
  R.squared = generics::glance(regression)$r.squared
  
  # make some summary graphs
  # raw timecourses
  # set the theme
  theme_set(theme_linus)
  
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

  (fitset =
        data %>%
        group_by(.data$sweep) %>%
        # repeat the extraction of the times to get the rising parts
        # this is not very efficient but its a fast calculation
        mutate(
          .,
          tscaled = .data$t - 0.0299,
          Istim = measure.timepoint(.data$I, .data$tscaled, 0.5, 0.02),
          Vsteady = measure.timepoint(.data$V, .data$tscaled, 0.4, 0.02),
          RMP = measure.timepoint(.data$V, .data$t, 0.02, 0.02),
          Vmin = min(.data$V),
          Vscaled = 1 - (.data$V - .data$RMP) / (.data$Vsteady - .data$RMP),
          tpeak = min(t[.data$V == .data$Vmin]),
          # the first point at which V is = Vmin, it mind find more for the garbage ones...
          sagcoef = .data$Vsteady / .data$Vmin
        ) %>%
       # for the fits do not take the first sweeps with sag.
       # use the first 20 ms of sweeps without
      filter(.data$sweep > 3,
             .data$sweep < 9,
             .data$tscaled >  0,
             .data$tscaled <  0.4))
    
    
    # Alternatively use the fitset to find the crossing of the normalised response
    # with 1/e of the final value (~36% in this case)
    # update and write to the summary table
    taux = summarise(
      fitset,
      taucross = mean(.data$tscaled[.data$Vscaled > 0.3578794 &
                                      .data$Vscaled < 0.3778794]))
    
    databysweep =
      left_join(databysweep, summarise(fitset, taux))
    
    # run the fitting in failsafe mode to generate table output
    safe_nls = purrr::safely(minpack.lm::nlsLM)
    #now do the fitting with a purr call and process the results
    fitset %<>%
      group_by(.data$sweep) %>%
      nest() %>%
      mutate(
        fit = purrr::map(
          .data$data,
          ~ safe_nls(
            formula = Vscaled ~ exp(-tscaled / tau),
            data = .,
            start = list(tau = 0.5)
          )
        ),
        fitbiexp = purrr::map(
          .data$data,
          ~ safe_nls(
            formula = Vscaled ~ (A1 * exp(-tscaled / tau1) + A2 * exp(-tscaled / tau2)),
            data = .,
            start = list(
              tau1 = 0.002,
              tau2 = 0.2,
              A1 = 0.9,
              A2 = 0.1
            )
          )
        )
      )
    
    # # Extract the results list from the safely output and overwrite the complicated nested column
    fitset$fit = purrr::transpose(fitset$fit)$result
    fitset$fitbiexp = purrr::transpose(fitset$fitbiexp)$result
    
    # get all the coefs and tidy stuff
    fitset %<>%
      mutate(coefs = purrr::map(.data$fitbiexp, broom::tidy),
             fitvalue = purrr::map(.data$fitbiexp, generics::augment))
    
    # Make a plot of the fits parameters
    params = 
    unnest(fitset, coefs) %>%
      select(.data$term, .data$estimate) %>%
      mutate(var = stringr::str_sub(.data$term, start = 1, end =1)) %>%
      ggplot(aes(x = .data$sweep, y = .data$estimate, colour = .data$term)) +
      facet_wrap( ~ var,
                  scales = "free_y",
                  strip.position = "left") +
      ggthemes::scale_color_solarized(accent = "blue") +
      geom_point() +
      theme(axis.title.y = element_blank())
    
    
    fits =
      fitset %>%
      unnest(c(.data$data), keep_empty = T) %>%
      ggplot(aes(x = .data$tscaled,
                 y = .data$Vscaled)) +
      labs(x = "t, s",
           y = "V, normalised",
           title = "Time Course Analysis, per sweep") +
      facet_wrap(
        ~ .data$sweep,
        as.table = T,
        #scales = "free",
        strip.position = "top",
      ) +
      scale_y_continuous(limits = c(-0.09, 1)) +
      geom_line(aes(group = .data$sweep)) +
      geom_line(
        data = unnest(fitset, .data$fitvalue),
        aes(y = .data$.fitted),
        linetype = 1,
        colour = "red"
      ) +
      geom_hline(yintercept = c(0, 1),
                 linetype = 3,
                 colour = "grey50") +
      geom_segment(
        aes(
          x = .data$taucross,
          xend = .data$taucross,
          y = 1 / exp(1),
          yend = -Inf
        ),
        colour = "blue",
        data = taux
      ) +
      geom_segment(
        aes(
          x = .data$taucross,
          xend = -Inf,
          y = 1 / exp(1),
          yend = 1 / exp(1)
        ),
        colour = "blue",
        data = taux
      ) +
      theme(legend.position = "none")
    
    # plot layout
    combined = (((IV) + (sag))/patchwork::wrap_plots(params, fits))
    
    ggsave(
      combined ,
      filename = paste0(abffile, ".passive.png"),
      width = 8.25,
      height = 9
    )
   
    # extract the fit parameters and write all the summary to file
    left_join(
      databysweep,
      fitset %>%
        unnest(.data$coefs) %>%
        select(.data$sweep)
    ) %>%
      readr::write_csv(file = paste0(abffile, ".sweeps.csv"))
}

