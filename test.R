library(tidyverse)
#source("./R/testpulse.protocol.R")


setwd("../SGNanalysis2/")


library(parallel)
library(doParallel)
reticulate::use_condaenv("analysis")
pyabf = reticulate::import("pyabf")
library(wctools)
library(patchwork)
ncores = parallel::detectCores()
cl <- makeCluster(ncores-1)
registerDoParallel(cl) 


## Start the analysis


new.VCtest =
  function(abffile) {
    tpulse =
      read.multisweep.pyth(abffile)
    
    # process off response data to calulate Cm
    offr =
      tpulse |>
      mutate(t = .data$t - 0.06156) %>%
      filter(.data$t > 0, .data$t < 0.05) %>%
      group_by(.data$sweep)
    
    
    offr.bysweep =
      offr %>%
      group_by(.data$sweep) %>%
      dplyr::summarise(Iss = wctools::measure.timepoint(.data$Imemb, .data$t, 0.045, win = 3))
    
    offr =
      left_join(offr, offr.bysweep) |>
      mutate(Isubs = (.data$Imemb - .data$Iss),
             Inet = abs(.data$Isubs))
    
    # calculate and merge the capacitance
    offrQ =
      offr %>%
      filter(.data$t < 0.004) |>
      group_by(.data$sweep) |>
      dplyr::summarise(Q = pracma::trapz(.data$t, .data$Inet))
    
    
    # processing of the on response
    tpulse =
      tpulse |>
      filter(.data$Vcmd0 != 0) |>
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
    
    tpulse.bysweep =
      left_join(tpulse.bysweep,
                offrQ) |>
      mutate(Cm = (.data$Q / abs(.data$V)) * 1000)
    
    
    
    plot1 =
      tpulse %>%
      filter(.data$Vcmd0 != 0, t < 0.02) %>%
      ggplot(aes(
        x = .data$t,
        y = .data$Isubs,
        colour = sweep
      )) +
      geom_line(aes(group = sweep)) +
      geom_hline(yintercept = 0,
                 linetype = 3,
                 colour = "grey50") +
      scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10 ^
                                                                                 .x)))
    
    # plot the point just after the max
    plot2 =
      tpulse |>
      left_join(tpulse.bysweep) |>
      # filter for  10 samples after peak
      filter(.data$Vcmd0 != 0, t > tpeak, t < (tpeak + 1e-04)) |>
      ggplot(aes(x = .data$t, y = .data$Isubs)) +
      geom_line(aes(group = sweep)) +
      geom_smooth(
        method = "lm",
        se = F,
        aes(group = sweep),
        formula = y ~ x
      ) +
      scale_x_log10() +
      geom_vline(xintercept = 0,
                 linetype = 3,
                 colour = "grey50") 
      geom_hline(yintercept = 0,
                 linetype = 3,
                 colour = "grey50") 
    
     Raset =
       tpulse |>
       left_join(tpulse.bysweep) |>
       # filter for  10 samples after peak
       filter(.data$Vcmd0 != 0, t > tpeak, t < (tpeak + 1e-04)) |>
       group_by(sweep) |>
       nest() |>
       summarise(fit = purrr::map(data, lm, formula = Isubs ~t)) |>
       mutate(res = purrr::map(fit, generics::tidy)) |>
       select(res, sweep) |>
       unnest(cols = res) |>
       filter(term == "(Intercept)") |>
       rename(I0 = estimate) |>
       select(sweep, I0)

     tpulse.bysweep =
       left_join(tpulse.bysweep, Raset)


     lm.Rm = lm(data = tpulse.bysweep, V ~ Iss)
     lm.Ra = lm(data = tpulse.bysweep, V ~ Ipeak)
     lm.Ra_2 = lm(data = tpulse.bysweep, V ~ I0)
     lm.Cm = lm(data = tpulse.bysweep, Q ~ V)

      estimates = list(Cm = abs(lm.Cm$coefficients[[2]])*1000,
                       Ra = abs(lm.Ra$coefficients[[2]]) *1000 ,
                       Rm = lm.Rm$coefficients[[2]])

     plotRm = ggplot(tpulse.bysweep,aes(V, Iss)) +
       geom_point() +
       labs(title = paste0("Rm estimate = ", round(lm.Rm$coefficients[[2]], digits = 2) , ", GOhm")) +
       geom_smooth(method = "lm")

     plotCm = ggplot(tpulse.bysweep,aes(V, Q)) +
       geom_point() +
       labs(title = paste0("Cm estimate = ", round(abs(lm.Cm$coefficients[[2]])*1000, digits = 2) , ", pF")) +
       geom_smooth(method = "lm")

     plotRa = ggplot(tpulse.bysweep,aes(V, Ipeak)) +
       geom_point() +
       labs(title = paste0("Ra estimate = ", round(abs(lm.Ra$coefficients[[2]]) *1000, digits = 2) , ", MOhm")) +
       geom_smooth(method = "lm")

     plotRa_2 = ggplot(tpulse.bysweep,aes(V, I0)) +
       geom_point() +
       labs(title = paste0("Ra estimate = ", round(abs(lm.Ra_2$coefficients[[2]]) *1000, digits = 2) , ", MOhm")) +
       geom_smooth(method = "lm")

    # plotA = ggplot(tpulse.bysweep)
     (plot1 + plot2) /(plotRm + plotCm + plotRa + plotRa_2)

       #list(lm.Rm, lm.Cm, lm.Ra)
      # tpulse.bysweep
  }


## process all Testpulse data
print("1/10 Processing Testpulses...")
list.files(path = "./Otof experi/VC testpulse/", pattern = "abf$", full.names = T)[108] %>%
  plyr::llply(.data = .,
              .fun = new.VCtest,
              .parallel = F)


stopCluster(cl) 
