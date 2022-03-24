# # experimental/ idea collection for doing a better job of isolating transient K
# 
# # read some test data
# #library(wctools)
# #library(tidyverse)
# # avialbility prot
# 
# Kavail = read.multisweep.pyth("../SGNanalysis2/Otof experi/VC availability/2021_09_22 SGNxotof P2 mouse95 middle33 cell95_0003.abf")
# 
# Kavail$V = Kavail$Vcmd0 - 80
# # time of the jump
# tjump = 0.67342
# 
# # plot to check
# # Kavail |> 
# #   filter(t > tjump) |> 
# #   ggplot(aes(x = t, y = Imemb, group = V)) +
# #   geom_line(aes(group = sweep))
# 
# Kavail =
#   left_join(Kavail,
#             Kavail |>
#               group_by(sweep) |>
#               summarise(Vpre = wctools:::measure.timepoint(V, t, 0.5, 0.01)))
# 
# # Do the substraction.
# # create a vector of the reference current
# refcurrent =
#   Kavail |>
#   filter(Vpre == -120, t > tjump) |>
#   select(Imemb) |> 
#   unlist() |> 
#   as.numeric()
# 
# datasubs = 
#   Kavail |> 
#   group_by(sweep) |> 
#   filter(t > tjump) |> 
#   mutate(Isubs = refcurrent - Imemb)
# 
# datasubs |> 
#   filter(t  <1.1) |> 
#   ggplot(aes(x = t, y = Isubs, group = Vpre)) +
#   geom_line(aes(group = sweep))
# 
# 
# datasubs |> 
#   group_by(Vpre) |> 
#   summarise(IKT = wctools::measure.timepoint(Isubs,t, tjump + 0.04, 0.05)) |> 
#   mutate(norm = 1-(IKT / max(IKT))) |> 
#   ggplot(aes(x= Vpre, y = norm)) +
#   geom_point()
# 
# 
