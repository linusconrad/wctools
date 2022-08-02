# This contains the function to wrangle the "passive" protocol data
#' Process an abf file of the "passive" protocol (updated version with rebound spike)
#' 
#' This function takes an .abf file of the passive protocol and returns to the folder a summary graph and a .csv.
#' These contain per-sweep summaries.
#' SR of 50000 is assumed.
#' Function is tailored to this protocol only.
#' @param abffile The File
#' @param Vjunc Junction potential to add
#' @param threshold Threshold for AP detection
#' @return a .png and a .csv with the summarised analysis. Named with the filename trunk of 'abffile'.
#' @import ggplot2
#' @import tidyr
#' @import patchwork
#' @import dplyr
#' @importFrom rlang .data
#' @export
process.passiveV2 = function(abffile, Vjunc, threshold, listobj = F) {
data = read.multisweep.pyth(abffile)
names(data) = c("t", "V", "I", "sweep")

#default to 0 mV AP detection
if(missing(threshold))
  threshold <- 0

# substract the junction potential
data %>% mutate(V = .data$V + Vjunc)

####### 
# Calculating the simple parameters
tjump = 0.02936
tstim = 0.32936

databysweep =
  data %>%
  group_by(.data$sweep) %>%
  dplyr::summarise(
    #tscaled = t - tjump,
    Istim = measure.timepoint(.data$I, .data$t, 0.3, 0.02),
    Vsteady = measure.timepoint(.data$V, .data$t, 0.32, 0.02),
    RMP = measure.timepoint(.data$V, .data$t, 0.02, 0.02),
    Vsag = min(.data$V),
    #Vscaled = 1 - (V - RMP) / (Vmin - RMP), #Vscaled to 1 to 0
    tsag = min(.data$t[.data$V == .data$Vsag]),
    # the first point at which V is = Vsag, it mind find more for the garbage ones...
    sagcoef = .data$Vsteady / .data$Vsag
  )

RMP = mean(databysweep$RMP)
# calculate the width of the peak and make a plot for QC
peaks =
  left_join(data, databysweep) %>%
  filter(.data$t > tjump, .data$t < tstim, .data$sweep < 10) %>%
  group_by(.data$sweep) %>%
  #make some handy timeseries based on the peak
  mutate(
    Vsubs = .data$V - .data$Vsteady,
    Vnorm = .data$Vsubs / mean(.data$Vsubs[.data$Vsubs == min(.data$Vsubs)]),
    #Vscaled to steady state
    Vtest = abs(.data$Vnorm - 0.5),
    tevent = t - mean(t[.data$Vsubs == min(.data$Vsubs)]),
    segment = sign(.data$tevent)
  )

# Do some TC fitting on this
peakfit = peaks
safefit = purrr::safely(minpack.lm::nlsLM)

peakfit %<>%
  group_by(.data$sweep, .data$Istim) %>%
  filter(.data$Vnorm > 0.03, .data$Vnorm <0.95) %>%
  nest() %>%
  mutate(fit = purrr::map(
    .data$data,
    ~ safefit(
      formula = Vnorm ~ SSbiexp(tevent, A1, lrc1, A2, lrc2),
      data = .,
      control = list(maxiter = 100)
    )
  ))

peakfit$fit = purrr::transpose(peakfit$fit)$result

peakfit =
  peakfit %>%
  mutate(
    params = purrr::map(.data$fit, broom::tidy),
    Sagcurves = purrr::map(fit, generics::augment, newdata = tibble(tevent = seq(0.01, 0.3, 0.01)))
  )

#number of succsessfully fit sweeps
nfits =
  peakfit |>
   mutate(test = class(.data$fit[[1]])) |> 
  rowwise() |>
  dplyr::mutate(fit.done = test != "NULL")

nfits = 
  sum(nfits$fit.done)

# Flow conrol: only do the following if the fits worked, else it will give errors
if (nfits != 0){

# get the fitted curves out
sagfitted = unnest(peakfit, .data$Sagcurves)
#substract the offset from the timeaxis again
sagfitted =
  left_join(sagfitted,
            peaks %>%
              group_by(.data$sweep) %>%
              dplyr::summarise(toffset = .data$t[.data$Vsubs == min(.data$Vsubs)]))

#return(peakfit)

sagfitted = 
  sagfitted %>%
  mutate(t = .data$tevent + .data$toffset)


peakfit %<>%
  select(.data$params) %>%
  unnest(.data$params) %>%
  select(.data$term, .data$estimate) %>%
  pivot_wider(names_from = .data$term,
              values_from = .data$estimate, names_prefix = "Sag.")
}
peak.params =
  peaks %>%
  group_by(.data$sweep, .data$segment) %>%
  summarise(thalf = mean(.data$t[.data$Vtest == min(.data$Vtest)])) %>%
  filter(.data$segment != 0) %>%
  mutate(segment = .data$segment + 1) %>%
  pivot_wider(
    values_from = .data$thalf,
    names_from = .data$segment,
    names_prefix = "saghalf"
  ) %>%
  mutate(sagwidth = .data$saghalf2 - .data$saghalf0)


peakQC =
  ggplot(peaks, aes(x = .data$t, y = .data$Vnorm)) +
  geom_hline(yintercept = c(0, 1),
             linetype = 3,
             colour = "grey50") +
  labs(title = "Membrane-Sag Peaks",
       y = "Vm, normalised",
       x = "t, s") +
  geom_line() +
  scale_y_continuous(limits = c(-0.2, 1)) +
  geom_segment(
    data = peak.params,
    aes(x = .data$saghalf0, xend = .data$saghalf2),
    y = 0.5,
    yend = 0.5,
    colour = "deeppink"
  ) +
  # geom_line(data = sagfitted, aes(y = .data$.fitted), colour = "deeppink")+
  facet_wrap( ~ sweep) +
  cowplot::theme_nothing() +
  theme(plot.title = element_text())

if (nfits != 0){
  peakQC = 
    peakQC +
    geom_line(data = sagfitted, aes(y = .data$.fitted), colour = "deeppink")
}

databysweep =
  left_join(databysweep, peak.params)

#calculate a linear regression for the input resistance
regression =
  stats::lm(Vsteady ~ Istim, data = databysweep)

coefs = broom::tidy(regression)
R.squared = generics::glance(regression)$r.squared

#make some summary graphs
trace =
  ggplot(data, aes(x = .data$t, y = .data$V)) +
  labs(x = "t, s", y = "V, mV") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse()+
  geom_hline(linetype = 3,
             colour = "grey50",
             yintercept = 0) +
  # vertical coursors to show measurements taken
  geom_line(aes(group = .data$sweep))  +
  geom_vline(
    xintercept = c(0.32, 0.02),
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
  ggpubr::stat_regline_equation(colour = "deeppink") +
  stat_smooth(method = "lm")

# summary of the sag coefficient
sag =
  ggplot(databysweep, aes(x = .data$Istim, y = .data$sagcoef)) +
  labs(
    x = latex2exp::TeX('$I_{stim}$, pA$'),
    y = latex2exp::TeX('$\\frac{V{ss}}{V_{peak}}$')
  ) +
  #scale_x_continuous(expand = c(0, 0)) +
  geom_hline(linetype = 3,
             colour = "grey50",
             yintercept = 1) +
  geom_smooth() +
  geom_point()

# plot of the peak widths
peakwidth = 
  ggplot(databysweep, aes(x = .data$Istim, y = .data$sagwidth*1000)) +
  labs(x = latex2exp::TeX("$I_{stim}$, pA"),
       y = "Sag width, ms") +
  geom_hline(linetype = 3,
             colour = "grey50",
             yintercept = 0) +
  geom_smooth() +
  geom_point()

rebAP = 
  data %>%
  filter(.data$t > tstim) %>%
  wctools::getAPstats("V", threshold) %>%
  mutate(
    latency = .data$tpeak - tstim,
    halfpoint1 = .data$width - .data$halfpoint,
    halfstart = .data$tpeak - .data$halfpoint1,
    halfend = .data$tpeak + .data$halfpoint)
 
# Analysis of the rebound spike and potential
databysweep =
  left_join(
    databysweep,
    data %>%
      group_by(.data$sweep) %>%
      summarise(Vreb = pracma::Mode(.data$V[.data$t > 0.35 & .data$t < 0.36]))
  ) %>%
  mutate(relreb = .data$Vreb - RMP)


databysweep =
  left_join(
    databysweep,
    rebAP)
# 
# Plot the properties of the rebound spike
AP =
  databysweep %>%
  select(sweep, .data$Vsteady, .data$Vpeak, .data$Vmin, .data$tmin, .data$width, .data$upstroke, .data$latency) %>%
  mutate(latency = .data$latency*1000,
         width = .data$width *1000,
         tmin = .data$tmin *1000) %>%
  pivot_longer(cols = c(.data$Vpeak, .data$Vmin, .data$tmin, .data$width, .data$upstroke, .data$latency)) %>%
  ggplot(aes(x = .data$Vsteady, y = .data$value)) +
  facet_wrap( ~ name, scales = "free_y", strip.position = "left") +
  #geom_line(linetype =3, colour = "grey50")+
  geom_point() +
  geom_smooth(se = T, linetype = 1, colour = "blue") +
  geom_hline(yintercept = 0,
             linetype = 3,
             colour = "grey50") +
  geom_vline(xintercept = RMP,
             linetype = 3,
             colour = "deeppink") +
  theme(axis.title.y = element_blank())


Vreb = 
  ggplot(databysweep, aes(x = .data$Istim, y = .data$relreb)) +
  labs(x = latex2exp::TeX("$I_{stim}$, pA"),
       y = "Rebound, mV") +
  geom_hline(linetype = 3,
             colour = "grey50",
             yintercept = 0) +
  geom_smooth() +
  geom_point() 

# extract the relaxation after the jump
tracereb =  data %>%
  filter(.data$t > tstim, .data$t < (tstim+0.02)) %>%
  ggplot(aes(x = .data$t, y = .data$V)) +
  labs(title = "Rebound Spike") +
  geom_hline(
    data = databysweep,
    aes(yintercept = RMP),
    colour = "deeppink",
    alpha = 0.6
  ) +
  geom_hline(
    data = databysweep,
    aes(yintercept = Vreb),
    colour = "deeppink",
    alpha = 0.6
  ) +
  geom_text(
    data = databysweep,
    aes(label = paste0(plyr::round_any(.data$relreb, 0.01))),
    x = 0.38,
    y = -20,
    colour = "deeppink"
  ) +
  geom_line(size = 0.2) +
  geom_segment(
    data = databysweep,
    mapping = aes(
      x = .data$halfstart,
      y = .data$Vhalf,
      xend = .data$halfend,
      yend = .data$Vhalf
    ),
    colour = "red"
  ) +
  geom_line() +
  scale_y_continuous(limits = c(-70, 20)) +
  facet_wrap( ~ sweep, nrow =2) +
  cowplot::theme_nothing() +
  theme(plot.title = element_text())

#plot layout
# combined_port = (peakQC + tracereb)/ (AP) / (IV + peakwidth) / (sag  + Vreb)
# combined_land = (((peakQC + tracereb)/AP)| (IV / peakwidth) | (sag  / Vreb)) +
#   patchwork::plot_layout(widths = c(2, 1, 1))


combined_land = ((tracereb / AP) |
  (peakQC / ((IV / peakwidth) |
               (sag  / Vreb))) ) +
  patchwork::plot_annotation(title = "CC Passive Response + Rebound",
                             subtitle = paste0(abffile, "\n Analysed on ", Sys.Date()))

# Write all to file
#Write all to file
utils::write.csv(databysweep, file = paste0(abffile, "CCpassivesummary.csv"))

if (nfits != 0) {
  utils::write.csv(peakfit, file = paste0(abffile, "Sag-relaxation-fit.csv"))
}
if (listobj == T) {
  return(list(data,
              databysweep,
              peaks,
              peakfit,
              sagfitted,
              peak.params))
}


ggsave(combined_land, file = paste0(abffile, "CCpassiveplot.png"),  width = 10, height = 6)
return(combined_land)
}