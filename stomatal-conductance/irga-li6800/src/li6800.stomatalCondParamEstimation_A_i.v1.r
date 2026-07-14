#'# Wheat stomatal kinetics (heat & drought)
#'## --- Ed Chaplin (USyd), *Apr--June 2025* ---
#'### Li-6800 STOMATAL CONDUCTANCE MODEL PARAMETER ESTIMATION FOR A - LIGHT OFF TO LIGHT ON (INCREASE)
#'
#' * Using the final model fitted to the data, obtain the parameters of interest
#' * These include: 
#' * A0, being the starting value of the A response (intercept/lower asymptote)
#' * A1, being the upper asymptote of the A responsE
#' * 
#'
#'### Load the libraries
library(knitr)
#knitr::spin("li6800.stomatalCondParamEstimation_A_i.v1.r")
library(ggplot2)
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, 
               "fig.width" = 7.5, "fig.height" = 10, "comment"=NA, fig.path = "figures/increase/")
library(dae)
library(dplyr)
library(brms) ##package used to fit bayesian regression model using Stan
library(tidyr)
library(emmeans)
library(tidybayes)
library(pracma)

#'### Set up constants and files
options(width=90)
load("../data/stomatalCondModel_A_i.v1.Rdata")
saveFile <- "../data/stomatalCondParamEstimation_A_i.v1.Rdata"

#'### Have a look at the final model summary and store fixed effects
summary(fit_on_A_final)

#'### Build a dataframe of model parameters
li6800.A.i.param <- expand.grid(Room_Temperature = levels(li6800.on.dat$Room_Temperature),
                             Watering = levels(li6800.on.dat$Watering))

#'### Use the emmeans package to get out the expected parameter estimates, along with HPD (credible intervals)
g0.A.emmeans <- as.data.frame(emmeans(fit_on_A_final,
                                    ~ Room_Temperature * Watering,
                                    nlpar = "g0"))
g0.A.emmeans <- g0.A.emmeans %>% rename_with(~paste0("g0_A_", .x), 
                                             -all_of(c("Room_Temperature", 
                                                       "Watering")))
g0.A.emmeans <- within(g0.A.emmeans, {
  g0_A_tabquant <- paste0(round(g0_A_emmean, 3), " [", 
                        round(g0_A_lower.HPD, 3), ", ", 
                        round(g0_A_upper.HPD, 3), "]")
  
})
li6800.A.i.param <- left_join(li6800.A.i.param, g0.A.emmeans)

g1.A.emmeans <- as.data.frame(emmeans(fit_on_A_final,
                                    ~ Room_Temperature * Watering,
                                    nlpar = "g1"))
g1.A.emmeans <- g1.A.emmeans %>% rename_with(~paste0("g1_A_", .x), 
                                             -all_of(c("Room_Temperature", 
                                                       "Watering")))
g1.A.emmeans <- within(g1.A.emmeans, {
  g1_A_tabquant <- paste0(round(g1_A_emmean, 3), " [", 
                        round(g1_A_lower.HPD, 3), ", ", 
                        round(g1_A_upper.HPD, 3), "]")
  
})
li6800.A.i.param <- left_join(li6800.A.i.param, g1.A.emmeans)

k.A.emmeans <- as.data.frame(emmeans(fit_on_A_final,
                                    ~ Room_Temperature * Watering,
                                    nlpar = "k"))
k.A.emmeans <- k.A.emmeans %>% rename_with(~paste0("k_A_", .x), 
                                       -all_of(c("Room_Temperature",
                                                 "Watering")))
k.A.emmeans <- within(k.A.emmeans, {
  k_A_tabquant <- paste0(round(k_A_emmean, 3), " [", 
                       round(k_A_lower.HPD, 3), ", ", 
                       round(k_A_upper.HPD, 3), "]")
  
})
li6800.A.i.param <- left_join(li6800.A.i.param, k.A.emmeans)

#'### Calculate the time taken to reach 95% Amax
#' t = -k(log(-log(0.95*g1-g0/g1-g0))-1)

g0_A_draws <- emmeans(fit_on_A_final,
                    ~ Room_Temperature * Watering,
                    nlpar = "g0") %>%
  gather_emmeans_draws()

g1_A_draws <- emmeans(fit_on_A_final,
                    ~ Room_Temperature * Watering,
                    nlpar = "g1") %>%
  gather_emmeans_draws()

k_A_draws <- emmeans(fit_on_A_final,
                   ~ Room_Temperature * Watering,
                   nlpar = "k") %>%
  gather_emmeans_draws()

combined_draws <- g0_A_draws %>%
  rename(g0 = .value) %>%
  left_join(g1_A_draws %>% rename(g1 = .value),
            by = c(".draw", "Room_Temperature", "Watering")) %>%
  left_join(k_A_draws %>% rename(k = .value),
            by = c(".draw", "Room_Temperature", "Watering"))

combined_draws <- combined_draws %>%
  mutate(t95 = -k*(log(-1*log((0.95*g1 - g0)/(g1-g0)))-1))

summary_t95 <- combined_draws %>%
  group_by(Room_Temperature, Watering) %>%
  median_hdi(t95) %>%
  rename(t95_emmean = t95) %>%
  rename(t95_lower.HPD = .lower) %>%
  rename(t95_upper.HPD = .upper)

summary_t95 <- within(summary_t95, {
  t95_tabquant <- paste0(round(t95_emmean, 3), " [", 
                           round(t95_lower.HPD, 3), ", ", 
                           round(t95_upper.HPD, 3), "]")
  
})
li6800.A.i.param <- left_join(li6800.A.i.param, summary_t95)

#'### Using t95 calculated above, use this to determine A95
#' Use the individual t95 values calculated above in combined_draws
#' A ~ g0 + (g1 - g0) * exp(-exp(((-tmin)/k)+1)),
combined_draws <- combined_draws %>%
  mutate(A95 = g0 + (g1 - g0) * exp(-exp(((-t95)/k)+1)),)

summary_A95 <- combined_draws %>%
  group_by(Room_Temperature, Watering) %>%
  median_hdi(A95) %>%
  rename(A95_emmean = A95) %>%
  rename(A95_lower.HPD = .lower) %>%
  rename(A95_upper.HPD = .upper)

summary_A95 <- within(summary_A95, {
  A95_tabquant <- paste0(round(A95_emmean, 3), " [", 
                         round(A95_lower.HPD, 3), ", ", 
                         round(A95_upper.HPD, 3), "]")
  
})
li6800.A.i.param <- left_join(li6800.A.i.param, summary_A95)


#'### Next, calculate the %Limitation. This involves integrating the nonlinear function
#' Can calculate this using the contents of combined_draws as well!

##define the antiderivative of our function
Aint <- function(tmin, g0, g1, k){
  g0*tmin - k*(g1-g0)*expint::expint_Ei(-exp(1 - tmin/k))
}

##get out the total integral, being the function evaluated between t=0 and t= max(tmin)
combined_draws <- combined_draws %>%
  mutate(Aint.tot.0 = Aint(0, g0, g1, k),
         Aint.tot.max = Aint(max(li6800.on.dat$tmin), g0, g1, k)) %>%
  mutate(Aint.tot.diff = Aint.tot.max - Aint.tot.0)
head(combined_draws %>% select(Aint.tot.0, Aint.tot.max, Aint.tot.diff))

##next, evaluate the integral of the difference between Amax and A, for the interval
##between 0 and t95
Aint.Amax <- function(tmin, g0, g1, k){
  g1*tmin - (g0*tmin - k*(g1-g0)*expint::expint_Ei(-exp(1 - tmin/k)))
}

##now evaluate the integral using Aint.Amax
combined_draws <- combined_draws %>%
  mutate(Aint.Amax.0 = Aint.Amax(0, g0, g1, k),
         Aint.Amax.t95 = Aint.Amax(t95, g0, g1, k)) %>%
  mutate(Aint.Amax.diff = Aint.Amax.t95 - Aint.Amax.0)

##calculate the ratio of Aint.Amax.diff and Aint.tot.diff
combined_draws <- combined_draws %>%
  mutate(Limitation = Aint.Amax.diff/Aint.tot.diff)

summary_limitation <- combined_draws %>%
  group_by(Room_Temperature, Watering) %>%
  median_hdi(Limitation) %>%
  rename(Limit_emmean = Limitation) %>%
  rename(Limit_lower.HPD = .lower) %>%
  rename(Limit_upper.HPD = .upper)

summary_limitation <- within(summary_limitation, {
  Limit_tabquant <- paste0(round(Limit_emmean*100, 3), " [", 
                         round(Limit_lower.HPD*100, 3), ", ", 
                         round(Limit_upper.HPD*100, 3), "]")
  
})
li6800.A.i.param <- left_join(li6800.A.i.param, summary_limitation)

#'### Format this in a nicer way for export to Excel
li6800.A.i.export <- li6800.A.i.param %>% select(Room_Temperature, Watering, 
                                             grep("tabquant", colnames(.))) %>%
  pivot_longer(cols = 3:ncol(.),
               names_to = "params",
               values_to = "values") %>%
  pivot_wider(names_from = c("Room_Temperature", "Watering"),
              values_from = "values")
li6800.A.i.export

#'### Save data and export
save.image(saveFile)
write.csv(li6800.A.i.export, "../data/li6800.A.i.param.export.csv", row.names = FALSE)
