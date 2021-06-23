# Testpulse analysis
#' This function takes a membrane test pulse protocol and returns parameters such as
#' Input resistance capacitance and so on.
#' @param abffile The File
#' @return a .csv with the summarised analysis. Named with the filename trunk of 'abffile'.
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
    filter(Vcmd0 != 0, t < 0.02) %>%
    ggplot(aes(x = t, y = Isubs, colour = sweep)) +
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
        filter(t < 0.02) %>%
        mutate(., Isubs.A = .data$Isubs) %>%
        dplyr::summarise(Q = pracma::trapz(.data$t, .data$Isubs))
    ) %>%
    mutate(Cm = (.data$Q /.data$V) * 1000)
  
  #return(plot1)
  # Fit the time constant
  # transform the timecourse to a 1to 0 scale
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
    nest() %>%
    mutate(fit = purrr::map(
      .data$data,
      ~ nls.multstart::nls_multstart(
        formula = Inorm ~ A1 * exp(-t / tau1) + A2 * exp(-t / tau2),
        data = .,
        start_lower = c(
          tau1 = 0.01,
          tau2 = 0.2,
          A1 = 0.5,
          A2 = 0.01
        ),
        start_upper = c(
          tau1 = 0.1,
          tau2 = 2,
          A1 = 5,
          A2 = 0.1
        ),
        convergence_count = 100,
        iter = 1000
      )
    ))
  
  # # extract the params
  tpulse.fit$fit = purrr::transpose(tpulse.fit$fit)$result
  #print(tpulse.fit$fit)
  
  # get all the coefs and tidy stuff
  tpulse.fit %<>%
    mutate(
      coefs = purrr::map(.data$fit, broom::tidy),
      fitvalue = purrr::map(.data$fit, generics::augment)
    )
  
  # Make a plot of the fits parameters
  # params =
  #   unnest(tpulse.fit, coefs) %>%
  #   select(term, estimate) %>%
  #   mutate(var = stringr::str_sub(term, start = 1, end = 1)) %>%
  #   ggplot(aes(x = sweep, y = estimate, colour = term)) +
  #   facet_wrap(~ var,
  #              scales = "free_y",
  #              strip.position = "left") +
  #   ggthemes::scale_color_solarized(accent = "blue") +
  #   geom_point() +
  #   theme(axis.title.y = element_blank())
  # 
  # fits =
  #   tpulse.fit %>%
  #   unnest(c(.data$data), keep_empty = T) %>%
  #   ggplot(aes(x = .data$t,
  #              y = .data$Inorm)) +
  #   labs(x = "t, s",
  #        y = "I, normalised") +
  #   scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10 ^.x)))+
  #   facet_wrap(
  #     ~ .data$sweep,
  #     as.table = T,
  #     #scales = "free",
  #     drop = T,
  #     strip.position = "top",
  #   ) +
  #   geom_line(aes(group = .data$sweep)) + 
  #   geom_line(
  #     data = unnest(tpulse.fit, .data$fitvalue),
  #     aes(y = .data$.fitted),
  #     linetype = 1,
  #     colour = "red"
  #   )
  
  
  # Add the fits to the by-sweep sumamry and calculate params based on tau.
  
  tpulse.bysweep =
    left_join(
      tpulse.bysweep,
      unnest(tpulse.fit, coefs) %>%
        select(.data$term, .data$estimate) %>%
        pivot_wider(names_from = .data$term, values_from = .data$estimate)
    ) %>%
    # calculate currents at t0
    mutate(I01 = .data$A1 * .data$Ipeak,
           I02 = .data$A2 * .data$Ipeak,
           I0sum = .data$I01 + .data$I02,
           Vt.tau1 = .data$tau1*1000,
           Vt.tau2 = .data$tau2*1000,
           Comp.soma = .data$A1/(.data$A1+ .data$A2),
           Ra = (abs(.data$V)/.data$I0sum) *1000) %>%
    select(.data$V, .data$Cm, .data$Vt.tau1, .data$Vt.tau2, .data$Comp.soma, .data$Ra)
  
  # collapse to one number
  tpulse.bysweep %>%
    pivot_longer(cols = c(2:6)) %>%
    group_by(.data$name) %>%
    summarise(value = mean(.data$value)) %>%
    pivot_wider() %>%
    write.csv(., file = paste0(abffile, "testpulse.csv"))
}