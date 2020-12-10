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
      n = seq_along(.data$FFT.V),
      freq = n * (50000 /lengthsweep),
      ZAP = Re((.data$FFT.V / .data$FFT.I)) * 1000,
      phase = Im(.data$FFT.V / .data$FFT.I)
    ) %>%
    ungroup %>%
    group_by(.data$freq) %>%
    summarise(realZAP = mean(Re(.data$ZAP)),
              imZAP = mean(.data$phase))
  
  ZAPdf
}