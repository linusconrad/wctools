# This contains the function to wrangle the ZAP protocol data

#' Calculation of complex impedance spectra from Current Clamp data
#' 
#' This function takes a dataframe with columns V and I from current clamp
#' and calculates a impedance amplitude profile (ZAP).
#' For the formulae etc see [this publication](https://journals.physiology.org/doi/abs/10.1152/jn.1986.55.5.995)
#' @param ZAPinput A dataframe with timeseries data, V (voltage, mV) and I (current, pA) and sweep (works best using read multisweep before).
#' assumes 50khz sampling 
#' @return A dataframe with columns realZAP (MOhm) and imZAP, the real and imaginary complex impedance values at each frequency
#' @export
calcZAP = function(ZAPinput){
  lengthsweep = length(ZAPinput$V/max(ZAPinput$sweep))
  
  ZAPdf = ZAPinput %>%
    group_by(sweep) %>%
    mutate(
      FFT.V = stats::fft(.data$V),
      FFT.I = stats::fft(.data$I),
      n = c(seq_along(.data$FFT.V))-1, #needs to start at 0!
      freq = n * 50000 /lengthsweep,
      ZAP = Re((.data$FFT.V / .data$FFT.I)) * 1000,
      phase = Im(.data$FFT.V / .data$FFT.I)
    ) %>%
    ungroup %>%
    group_by(.data$freq) %>%
    summarise(realZAP = mean(Re(.data$ZAP)),
              imZAP = mean(.data$phase))
  
  ZAPdf
}

#' Processing function for ZAP rawdata
#' 
#' This function takes a path to a ZAP recording and writes the ZAP profile to file
#' Its intended to be used in l_ply calls
#' @param file ABF file with ZAP data (3 sweeps, 50000 hz SR)
#' @return writes to file ZAP csv file and plot png in the path of the input file
#' @export
process.ZAP = function(file){
  # 1. get data, SR is 50000 always
  data = read.multisweep(file, 50000)
  print(paste0("Done reading file ", file,"..."))
  names(data) = c("V", "I", "t", "sweep")
  # 2. Calculate it 
  print("Calculating ZAP...")
  ZAPdata = calcZAP(data)
  # 3. Write to file
  # Truncate the spectra to the upper limit of the current input frequency
  # Files get huge otherwise
  
  ZAPdata %<>%
    dplyr::filter(.data$freq < 50) %>%
    mutate(filename = file)
  
  print("...done")
  print("Writing ZAP...")
  readr::write_csv(ZAPdata, file = paste0(file, ".ZAPoutput.csv"))
  print("...done")
  
  
  # plot the rawdata
  raw =
    data %>%
    filter(.data$t > 0.5 , .data$t < max(30), .data$sweep == 1) %>%
    ggplot(., aes(x = .data$t, y = .data$V)) +
    labs(y = "V, mV",
         x = "t, s",
         title = "Raw Traces") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(breaks = c(-100, 0)) +
    ggthemes::scale_color_few() +
    facet_wrap(~ .data$sweep,
               ncol = 1) +
    geom_line() +
    geom_hline(linetype = 3,
               colour = "grey50",
               yintercept = c(0, -100)) +
    cowplot::theme_nothing() +
    theme(axis.text.y = element_text(),
          panel.border = element_rect(colour = "grey50", fill = NA),
          plot.title = element_text(),
          text = element_text())
  # plot the zap spectra
  # First add an SG filter to plot
  ZAPdata %<>%
    mutate(smoothR = signal::filter(signal::sgolay(p=1, n=67, m=0), .data$realZAP),
           smoothIM = signal::filter(signal::sgolay(p=1, n=67, m=0), .data$imZAP))
  
  theme_set(theme_linus)
  plot1 = ZAPdata %>%
    filter(.data$freq > 1, .data$freq < 50, .data$realZAP > 0) %>%
    ggplot(aes(x = .data$freq, y = .data$realZAP)) +
    labs(x = "f, Hz",
         title = "Unprocessed ZAP",
         y = latex2exp::TeX("Z, M$\\Omega$"))+
    geom_line() +
    geom_line(aes(y = .data$smoothR), colour = "blue")
  
  plot2 =  ZAPdata %>%
    pivot_longer(all_of(c("smoothR", "smoothIM")),
                 names_to = "component",
                 values_to = "ZAP") %>%
    filter(.data$freq > 1, .data$freq < 20) %>%
    ggplot(aes(x = .data$freq, y = .data$ZAP)) +
    labs(title = "Smoothed Spectra",
         y = " ",
         x = "f, Hz") +
    geom_line(colour = "blue") +
    facet_wrap(
      ~ component,
      labeller = labeller(
        component = c('smoothR' = "Z, MOhm",
                      'smoothIM' = "Phase")
      ),
      scales = "free",
      ncol = 1,
      strip.position = "left"
    )
  print("Creating summary page...")
  cowplot::plot_grid(
    raw,
    egg::ggarrange(plot1,
                   plot2,
                   ncol = 2),
    ncol = 1,
    align = "h",
    axis = "lr",
    rel_heights = c(0.2,1)
  ) %>%
    ggpubr::annotate_figure(top = ggpubr::text_grob(
        "Impedance Amplitude Profile",
        size = 14),
      bottom = ggpubr::text_grob(paste0("Data source: \n",
                                file))
    ) %>%
    ggsave(filename = paste0(file, "ZAPplot.png"),
      width = 6,
      height = 4
    )
}



