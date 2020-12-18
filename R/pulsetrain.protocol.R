# Processing function for pulsetrain protocols


process.pulsetrain = function(abf){
  data = read.multisweep.pyth(abf)
  # Add the AP information
  data %<>% add.AP(., "Vmemb")
  
  
  # get the resting potential
  # from holding @ zero before the stimulation
  data.by.sweep =
    data %>%
    group_by(.data$sweep) %>%
    summarise(measure.timepoint(.data$Vmemb, .data$t, 0.02, 0.01))
  
  # add a variable to the rawdata that sets RMP to 0
  data %>% 
    group_by(.data$sweep) %>%
    mutate(Vm0 = .data$Vmemb - measure.timepoint(.data$Vmemb, .data$t, 0.02, 0.01))
}