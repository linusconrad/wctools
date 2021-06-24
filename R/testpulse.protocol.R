# Testpulse analysis
#' This function takes a membrane test pulse protocol and returns parameters such as
#' Input resistance capacitance and so on.
#' @param abffile The File
#' @return no object, but a .csv with the summarised analysis. Named with the filename trunk of 'abffile'. Written to the folder of the rawdata.
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @importFrom rlang .data
#' @export
process.VCtest = function(abffile) {
  tpulse =
    read.multisweep.pyth(abffile) %>%
    filter(.data$Vcmd0 != 0)
  
  tpulse %<>%
    mutate(t = .data$t - 0.01156) %>%
    filter(.data$t > 0, .data$t < 0.05) %>%
    group_by(.data$sweep)
  
  tpulse.bysweep =
    tpulse %>%
    group_by(.data$sweep) %>%
    dplyr::summarise(
      V = wctools::measure.timepoint(.data$Vcmd0, .data$t, 0.045, win = 3),
      Iss = wctools::measure.timepoint(.data$Imemb, .data$t, 0.045, win = 3)
    )
  # 
  tpulse =
    left_join(tpulse, tpulse.bysweep) %>%
    mutate(Isubs = (.data$Imemb - .data$Iss))
  
  tpulse.bysweep =
    left_join(
      tpulse.bysweep,
      tpulse %>%
        mutate(Inet = .data$Isubs * sign(.data$V)) %>%
        summarise(Ipeak = max(.data$Inet),
                  tpeak = mean(t[.data$Inet == max(.data$Inet)]))
    )
  
  plot1 =
    tpulse %>%
    filter(.data$Vcmd0 != 0, t < 0.02) %>%
    ggplot(aes(x = .data$t, y = .data$Isubs, colour = sweep)) +
    geom_line(aes(group = sweep)) +
    geom_hline(yintercept = 0,
               linetype = 3,
               colour = "grey50") +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10 ^.x)))
  
  # Calculate the capcitance with the integral
  tpulse.bysweep =
    left_join(
      tpulse.bysweep,
      tpulse %>%
        select(.data$sweep, .data$V, .data$Isubs, .data$t) %>%
        filter(.data$t < 0.02) %>%
        mutate(., Isubs.A = .data$Isubs) %>%
        dplyr::summarise(Q = pracma::trapz(.data$t, .data$Isubs))
    ) %>%
    mutate(Cm = (.data$Q /.data$V) * 1000)
  
  #return(plot1)
  # Fit the time constant
  # transform the timecourse to a 1to 0 scale
  
  safe_nls = purrr::safely(minpack.lm::nlsLM)
  
  tpulse.fit =
    tpulse %>%
    mutate(
      Inet = .data$Isubs * sign(.data$V),
      Ipeak = max(.data$Inet),
      tpeak = mean(.data$t[.data$Inet == max(.data$Inet)]),
      Inorm = .data$Inet / .data$Ipeak
    ) %>%
    filter(.data$t > .data$tpeak) %>%
    select(.data$t, .data$Inorm) %>%
    #
    filter(.data$Inorm > 0.03, .data$Inorm <0.97) %>%
    nest() %>%
    # mutate(fit = purrr::map(
    #   .data$data,
    #   ~ safe_nls(
    #     formula = Inorm ~ A1 * exp(-exp(lrc1)*t) + A2 * exp(-exp(lrc2)*t),
    #     data = .,
    #     minpack.lm::nls.lm.control(factor = 10),
    #     start = list(
    #       lrc1 = 9,
    #       lrc2 = 0.15,
    #       A1 = 10,
    #       A2 = 0.1
    #     )
    #   )
    # ))
  mutate(fit = purrr::map(
    .data$data,
    ~ safe_nls(
      formula = Inorm ~ A1 * exp(-t / tau1) + A2 * exp(-t / tau2),
      data = .,
      #minpack.lm::nls.lm.control(maxiter = 10000),
      start = list(
        tau1 = 6e-5,
        tau2 = 0.0009,
        A1 = 5,
        A2 = 0.1
      )
    )
  ))
  # 
  # # extract the params
  tpulse.fit$fit = purrr::transpose(tpulse.fit$fit)$result
  tpulse.fit
  

  # get all the coefs and tidy stuff
  tpulse.fit = 
  tpulse.fit %>%
    mutate(
      coefs = purrr::map(.data$fit, broom::tidy),
      fitvalue = purrr::map(.data$fit, generics::augment)
    )
  
  #Make a plot of the fits parameters
  params =
    unnest(tpulse.fit, .data$coefs) %>%
    select(.data$term, .data$estimate) %>%
    mutate(var = stringr::str_sub(.data$term, start = 1, end = 1)) %>%
    ggplot(aes(x = sweep, y = .data$estimate, colour = .data$term)) +
    facet_wrap(~ .data$var,
               scales = "free_y",
               strip.position = "left") +
    ggthemes::scale_color_solarized(accent = "blue") +
    geom_point() +
    theme(axis.title.y = element_blank())

  fits =
    tpulse.fit %>%
    unnest(c(.data$data), keep_empty = T) %>%
    ggplot(aes(x = .data$t,
               y = .data$Inorm)) +
    labs(x = "t, s",
         y = "I, normalised") +
    scale_x_log10(labels = scales::label_number_auto())+
    facet_wrap(
      ~ .data$sweep,
      as.table = T,
      #scales = "free",
      drop = T,
      strip.position = "top",
    ) +
    geom_line(aes(group = .data$sweep)) +
    geom_line(
      data = unnest(tpulse.fit, .data$fitvalue),
      aes(y = .data$.fitted),
      linetype = 1,
      colour = "red"
    )

  # 
  # # # Add the fits to the by-sweep sumamry and calculate params based on tau.
  # 
  tpulse.bysweep =
    left_join(
      tpulse.bysweep,
      unnest(tpulse.fit, .data$coefs) %>%
        select(.data$term, .data$estimate) %>%
        pivot_wider(names_from = .data$term,
                    values_from = .data$estimate,
                    names_prefix = "VCtest.")
    ) %>%
    # calculate currents at t0
    mutate(I01 = .data$VCtest.A1 * .data$Ipeak,
           I02 = .data$VCtest.A2 * .data$Ipeak,
           I0sum = .data$I01 + .data$I02,
           Vt.tau1 = .data$VCtest.tau1*1000,
           Vt.tau2 = .data$VCtest.tau2*1000,
           Comp.soma = .data$VCtest.A1/(.data$VCtest.A1+ .data$VCtest.A2),
           Ra = (abs(.data$V)/.data$I0sum) *1000) %>%
    select(.data$V, .data$Cm, .data$Vt.tau1, .data$Vt.tau2, .data$Comp.soma, .data$Ra)

  # collapse to one number
  utils::write.csv(tpulse.bysweep, file = paste0(abffile, "testpulse.csv"))
  (plot = (fits + params) + patchwork::plot_annotation(title = "VC testpulse",
                                           subtitle = paste0(abffile, "\n analysed on ", Sys.Date())))
  
  
  ggsave(plot, width = 6, height = 4, filename = paste0(abffile, ".testpulse.png"))
  return(plot)
   }