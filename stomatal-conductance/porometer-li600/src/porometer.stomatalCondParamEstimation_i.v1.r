#'# Wheat stomatal kinetics (heat & drought)
#'## --- Ed Chaplin (USyd), *Apr--June 2025* ---
#'### Li-600 STOMATAL CONDUCTANCE MODEL PARAMETER ESTIMATION - LIGHT OFF TO LIGHT ON (INCREASE)
#'
#' * Using the final model fitted to the data, obtain the parameters of interest
#' * These include: 
#' * k_i, being the time constant for gs to increase to Gsmax
#' * lambda, being the lag in the response time of gs to a step increase in light
#' * g0, being the starting value of the gs response (intercept/lower asymptote)
#' * g1, being the upper asymptote of the gs response (this will be used to calculate iWUE)
#' * Slmax, a derived trait from McAusland paper
#' * 
#'
#'### Load the libraries
library(knitr)
#knitr::spin("porometer.stomatalCondParamEstimation_i.v1.r")
library(ggplot2)
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, 
               "fig.width" = 7.5, "fig.height" = 10, "comment"=NA, fig.path = "figures/increase/")
library(dae)
library(dplyr)
library(brms) ##package used to fit bayesian regression model using Stan
library(tidyr)
library(emmeans)
library(tidybayes)

#'### Set up constants and files
options(width=90)
load("../data/porometer.anal_i.v1.Rdata")
saveFile <- "../data/stomatalCondParamEstimation_porometer_i.v1.Rdata"

#'### Have a look at the final model summary and store fixed effects
summary(fit_on_final)

#'### Build a dataframe of model parameters
li600.i.param <- expand.grid(Room_Temperature = levels(li600.on.dat$Room_Temperature),
                             Device_Watering = levels(li600.on.dat$Device_Watering),
                             Line = levels(li600.on.dat$Line))

#'### Use the emmeans package to get out the expected parameter estimates, along with HPD (credible intervals)
g0.emmeans <- as.data.frame(emmeans(fit_on_final,
                                    ~ Room_Temperature * Device_Watering * Line,
                                    nlpar = "g0"))
g0.emmeans <- g0.emmeans %>% rename_with(~paste0("g0_", .x), 
                                         -all_of(c("Room_Temperature", 
                                                   "Device_Watering",
                                                   "Line")))
g0.emmeans <- within(g0.emmeans, {
  g0_tabquant <- paste0("'", round(g0_emmean, 3), "\n[", 
                        round(g0_lower.HPD, 3), ", ", 
                        round(g0_upper.HPD, 3), "]'")
  
})
li600.i.param <- left_join(li600.i.param, g0.emmeans)

g1.emmeans <- as.data.frame(emmeans(fit_on_final,
                                    ~ Room_Temperature * Device_Watering * Line,
                                    nlpar = "g1"))
g1.emmeans <- g1.emmeans %>% rename_with(~paste0("g1_", .x), 
                                         -all_of(c("Room_Temperature", 
                                                   "Device_Watering",
                                                   "Line")))
g1.emmeans <- within(g1.emmeans, {
  g1_tabquant <- paste0("'", round(g1_emmean, 3), "\n[", 
                        round(g1_lower.HPD, 3), ", ", 
                        round(g1_upper.HPD, 3), "]'")
  
})
li600.i.param <- left_join(li600.i.param, g1.emmeans)

k.emmeans <- as.data.frame(emmeans(fit_on_final,
                                    ~ Line,
                                    nlpar = "k"))
k.emmeans <- k.emmeans %>% rename_with(~paste0("k_", .x), 
                                       -all_of(c("Line")))
k.emmeans <- within(k.emmeans, {
  k_tabquant <- paste0("'", round(k_emmean, 3), "\n[", 
                        round(k_lower.HPD, 3), ", ", 
                        round(k_upper.HPD, 3), "]'")
  
})
li600.i.param <- left_join(li600.i.param, k.emmeans)

lambda.emmeans <- as.data.frame(emmeans(fit_on_final,
                                   ~ 1,
                                   nlpar = "lambda"))
lambda.emmeans <- lambda.emmeans %>% rename_with(~paste0("lambda_", .x))
lambda.emmeans <- within(lambda.emmeans, {
  lambda_tabquant <- paste0("'", round(lambda_emmean, 3), "\n[", 
                       round(lambda_lower.HPD, 3), ", ", 
                       round(lambda_upper.HPD, 3), "]'")
  
})
li600.i.param <- cbind(li600.i.param, data.frame(lambda.emmeans[,2:5]))

#'### Using these parameters, calculate Slmax
#' * Slmax looks to be reported in the units mmol m^-2 s^-2
#' * This involves multiplying g1 and g0 by 1000 and also k by 60
#' * I can calculate this using the posterior draws for each parameter and then combine them together

g0_draws <- emmeans(fit_on_final,
                    ~ Room_Temperature * Device_Watering * Line,
                    nlpar = "g0") %>%
  gather_emmeans_draws()

g1_draws <- emmeans(fit_on_final,
                    ~ Room_Temperature * Device_Watering * Line,
                    nlpar = "g1") %>%
  gather_emmeans_draws()

k_draws <- emmeans(fit_on_final,
                   ~ Line,
                   nlpar = "k") %>%
  gather_emmeans_draws()

combined_draws <- g0_draws %>%
  rename(g0 = .value) %>%
  left_join(g1_draws %>% rename(g1 = .value),
            by = c(".draw", "Room_Temperature", "Device_Watering", "Line")) %>%
  left_join(k_draws %>% rename(k = .value),
            by = c(".draw", "Line"))

combined_draws <- combined_draws %>%
  mutate(Slmax = (g1*1000 - g0*1000) / (exp(1) * k*60))

summary_slopes <- combined_draws %>%
  group_by(Room_Temperature, Device_Watering, Line) %>%
  median_hdi(Slmax) %>%
  rename(Slmax_emmean = Slmax) %>%
  rename(Slmax_lower.HPD = .lower) %>%
  rename(Slmax_upper.HPD = .upper)

summary_slopes <- within(summary_slopes, {
  Slmax_tabquant <- paste0("'",round(Slmax_emmean, 3), "\n[", 
                           round(Slmax_lower.HPD, 3), ", ", 
                           round(Slmax_upper.HPD, 3), "]'")
  
})

li600.i.param <- left_join(li600.i.param,
                            summary_slopes %>% select(Room_Temperature,
                                                      Device_Watering,
                                                      Line,
                                                      Slmax_emmean,
                                                      Slmax_lower.HPD,
                                                      Slmax_upper.HPD,
                                                      Slmax_tabquant))
li600.i.param

#'### Format this in a nicer way for export to Excel
#'Keep the columns in the same order as Li6800 results and order first by parameter and then Line
li600.i.export <- li600.i.param %>% select(Room_Temperature, Device_Watering, Line, 
                                             grep("tabquant", colnames(.))) %>%
  pivot_longer(cols = 4:ncol(.),
               names_to = "params",
               values_to = "values") %>%
  pivot_wider(names_from = c("Room_Temperature", "Device_Watering"),
              values_from = "values")
li600.i.export$params <- factor(li600.i.export$params, levels = c("g0_tabquant",
                                                                  "g1_tabquant",
                                                                  "k_tabquant",
                                                                  "lambda_tabquant",
                                                                  "Slmax_tabquant"))
li600.i.export <- li600.i.export %>% arrange(params, Line) %>%
  select("params", "Line", "CER7_Heat_D2_Limited", "CER8_Ambient_D2_Limited", 
         "CER7_Heat_D1_Unlimited", "CER8_Ambient_D1_Unlimited")
li600.i.export

#'### Save data and export
save.image(saveFile)
write.csv(li600.i.export, "../data/li600.i.param.export.csv", row.names = FALSE,
          quote = TRUE)
