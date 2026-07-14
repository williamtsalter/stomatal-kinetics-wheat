#'# Wheat stomatal kinetics (heat & drought)
#'## --- Ed Chaplin (USyd), *Apr--June 2025* ---
#'### Li-6800 STOMATAL CONDUCTANCE MODEL FITTING - LIGHT ON TO LIGHT OFF (DECREASE)
#'
#' * This script focuses on the portion of the gs response when the light is swtiched off (_i script was when light switched on)
#' * This script will attempt to fit the methods of Woning et al (2026) and McAusland et al (2016) using the script supplied by Silvere Vialet-Chabrand
#' * The models fitted in the script look to be the Cumulative Distribution Weibull (CDWeibull) models described in Eqn 4.1 of Woning et al (2026)
#' * I've also adapted the script to incorporate the McAusland model
#' * Scripts supplied by Silvere can be found in the gs_fit.zip directory
#' * There are two scripts in this directory: fit_per_curve.r and fit_per_species.r
#' * Fit_per_curve looks to take each individual plant and fit a curve to it, 
#' while fit_per_species fits some parameters (tau and theta) at a species (treatment in our setting) level and others at a replicate level
#' * In an email from Silvere (shared by Chris and put in the \textit{reference-material} directory), 
#' Silvere prefers the fit_per_species approach, as do I, as this seems to align better with conventional analysis approaches
#'
#'### Changes in v2
#' * Using learnings from li600 data, returning to this script to fit the treatments as fixed effects in the same model, rather than fitting each response separately
#' * Only fitting McAusland parameterisation in this case and ditching the need to fit multiple models
#'
#'### Load the libraries
library(knitr)
#knitr::spin("li6800.stomatalCondModel_d.v2.r")
library(ggplot2)
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, 
               "fig.width" = 7.5, "fig.height" = 10, "comment"=NA, fig.path = "figures/decrease/")
library(dae)
library(dplyr)
library(growthPheno) # calcLagged
packageVersion("growthPheno")
library(brms) ##package used to fit bayesian regression model using Stan

#'### Set up constants and files
options(width=90)
load("../data/li6800.v1.Rdata")
saveFile <- "../data/stomatalCondModel_d.v2.RData"

#'### Define TW to be unique for each treatment
li6800.dat$TW <- paste(li6800.dat$Temperature, li6800.dat$Watering, sep="_")
li6800.dat$TW <- factor(li6800.dat$TW)

#'### Define th and tmin, to be time in hours and minutes respectively
li6800.dat$th <- li6800.dat$xElapsed/3600
li6800.dat$tmin <- li6800.dat$xElapsed/60

#'### Add other useful variables
li6800.dat$Room <- "CER7"
li6800.dat$Room[li6800.dat$Temperature == "Amb"] <- "CER8"
li6800.dat <- within(li6800.dat, {
  Room_Temperature <- paste(Room, Temperature, sep="_")
  Room <- factor(Room)
  Room_Temperature <- factor(Room_Temperature)
  Watering <- factor(Watering)
})

#'###find the minimum time when the light was on
min(li6800.dat$tmin[li6800.dat$Light == "On" & li6800.dat$tmin > 0]) ##60.5
max(li6800.dat$tmin[li6800.dat$Light == "On"]) ##120.0167
li6800.dat$tmin[li6800.dat$Light == "Off" &
                  li6800.dat$tmin > 60.5 &
                  li6800.dat$tmin < 120.0167] ##nothing in between these two times


#'### Explore the gs profiles of just the light on portion of the experiment
plotProfiles(li6800.dat %>% filter(Light == "Off" & tmin > 120), response = "gsw",
             individuals = "TWR", times = "tmin",
             facet.x = "Room_Temperature",
             facet.y = "Watering",
             y.title = "gsw",
             x.title = "minutes",
             title = "gsw",
             breaks.spacing.x = 5,
             alpha = 0.75)

#'### Following *fit_per_species.R*
#' * Script begins by specifying a stomatal conductance model and setting priors on the parameters
#' * I have played a little with the priors specified in the script by Silvere, to be more representative of this dataset
#' * Following discussion with Ed on 22/04, perhaps we want to select a different ``training" dataset for each treatment combination - need to explore what impact this has
#' * I haven't implemented the above dot point, as the model is refit for each treatment - it's just the compliation of the tool that fits the model under the hood that doesn't need to be re-initialised and keeps the same structure
#' 

#'#### Stomatal conductance model - McAusland from McAusland et al 2016
nlform_mc <- bf(gsw ~ g0 + (g1 - g0) * exp(-exp(((lambda - tmin)/k)+1)),
             g0 ~ 1 + Room_Temperature*Watering +
               (1 | Rep),
             g1 ~ 1 + Room_Temperature*Watering +
               (1 | Rep),
             k ~ 1+ Room_Temperature*Watering,
             lambda ~ 1+ Room_Temperature*Watering,
             nl = TRUE)

priors_mc <- c(prior(normal(0.3, 0.1), nlpar = "g0"),
            prior(normal(0, 0.02), nlpar = "g1"),
            prior(normal(0, 0.1), class = "sd", nlpar = "g0"),
            prior(normal(0, 0.02), class = "sd", nlpar = "g1"),
            prior(normal(3, 0.1), nlpar = "k", lb = 0),
            prior(normal(0.5, 0.25), nlpar = "lambda", lb = 0),
            prior(student_t(3, 0, 0.01), class = "sigma"))

plot(0.3 + (0 - 0.3) * exp(-exp(((0.5 - seq(0, 60, 0.5))/3) + 1)) ~
       seq(0, 60, 0.5))

#'#### Put finishing touches on dataframe
#' Subset for light on portion and remove th == 0 if it's present, as the light was still on at that point
li6800.off.dat <- subset(li6800.dat, Light == "Off" & tmin > 120)
li6800.off.dat %>% group_by(Rep) %>% summarise(min(tmin))
##Align each rep to make sure time 0 is the same
li6800.off.dat <- li6800.off.dat %>% group_by(Rep) %>% mutate(tmin = tmin - min(tmin)) %>% as.data.frame()
li6800.off.dat %>% group_by(Rep) %>% summarise(min(tmin)) %>% as.data.frame()

plotProfiles(li6800.off.dat, response = "gsw",
             individuals = "TWR", times = "tmin",
             y.title = "gsw",
             x.title = "minutes",
             title = "gsw",
             breaks.spacing.x = 5,
             alpha = 0.75)

#'#### Fit the model
# Configure and compile the model - McAusland
fit_off_max <- brm(
  formula = nlform_mc,
  data = li6800.off.dat,
  prior = priors_mc,
  chains = 6,
  cores = 6,
  iter = 10000,
  warmup = 6000,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 12,
    metric = "diag_e"
  )
)

print(summary(fit_off_max))
plot(fit_off_max, ask = FALSE)
pairs(fit_off_max)
mcmc_plot(fit_off_max, type = "trace")
rstan::check_hmc_diagnostics(fit_off_max$fit)
loo(fit_off_max)

pp_check(fit_off_max)  # compares density of the observed data with data generated from posterior predictive distribution
pp_check(fit_off_max, type = "error_hist", ndraws = 12)
pp_check(fit_off_max, type = "scatter_avg", ndraws = 1000)
pp_check(fit_off_max, type = "stat_2d")
pp_check(fit_off_max, type = "loo_ribbon")

plot(gsw~tmin, li6800.off.dat, pch=16)
points(fitted(fit_off_max)[,1]~tmin, li6800.off.dat, pch=16, col="red")

##make the above plot prettier
li6800.off.plot <- li6800.off.dat
li6800.off.plot$fitted <- fitted(fit_off_max)[,1]
li6800.off.plot$fitted.fixed <- fitted(fit_off_max, re_formula = NA)[,1]

ggplot(li6800.off.plot, aes(x= tmin, y = gsw, group = Rep)) +
  geom_point() +
  facet_grid(Room_Temperature ~ Watering) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red")

ggplot(li6800.off.plot, aes(x= tmin, y = gsw, group = Rep)) +
  geom_point() +
  facet_grid(Room_Temperature ~ Watering) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red") +
  geom_point(aes(y=fitted.fixed), colour = "blue") +
  geom_line(aes(y=fitted.fixed), colour = "blue")

#'#### Testing contribution of effects in the model
# nlform_mc_reduced <- bf(gsw ~ g0 + (g1 - g0) * exp(-exp(((lambda - tmin)/k)+1)),
#                         g0 ~ 1 + Room_Temperature*Watering +
#                           (1 | Rep),
#                         g1 ~ 1 + Room_Temperature*Watering +
#                           (1 | Rep),
#                         k ~ 1,
#                         lambda ~ 1 + Watering,
#                         nl = TRUE)
# 
# # Fit the model
# fit_off_reduced <- brm(
#   formula = nlform_mc_reduced,
#   data = li6800.off.dat,
#   prior = priors_mc,
#   chains = 6,
#   cores = 6,
#   iter = 10000,
#   warmup = 6000,
#   control = list(
#     adapt_delta = 0.99,
#     max_treedepth = 12,
#     metric = "diag_e"
#   )
# )
# 
# print(summary(fit_off_reduced))
# loo(fit_off_reduced)
# loo(fit_off_current)
# plot(gsw~tmin, li6800.off.dat, pch=16)
# points(fitted(fit_off_reduced)[,1]~tmin, li6800.off.dat, pch=16, col="red")
# coef(fit_off_reduced)
# 
# fit_off_current <- add_criterion(fit_off_current, "loo")
# fit_off_reduced <- add_criterion(fit_off_reduced, "loo")
# 
# loo_compare(fit_off_current, fit_off_reduced)
##Lambda: Drop Room_Temperature:Watering interaction
##k: Drop Room_Temperature:Watering interaction
#g1: Retain Room_Temperature:Watering interaction
#g0: Retain Room_Temperature:Watering interaction

##Lambda: Drop Room_Temperature main effect
##lambda: Retain Watering main effect
##k: Drop Room_Temperature main effect
##k: Retain Watering main effect
#fit_off_current <- fit_off_reduced
#print(summary(fit_off_current))

#'#### Fit the final model
#g0: Room_Temperature:Watering
#g1: Room_Temperature:Watering
#k: Watering
#lambda: Watering
nlform_mc_final <- bf(gsw ~ g0 + (g1 - g0) * exp(-exp(((lambda - tmin)/k)+1)),
                        g0 ~ 1 + Room_Temperature*Watering +
                          (1 | Rep),
                        g1 ~ 1 + Room_Temperature*Watering +
                          (1 | Rep),
                        k ~ 1 + Watering,
                        lambda ~ 1 + Watering,
                        nl = TRUE)

# Fit the final model
fit_off_final <- brm(
  formula = nlform_mc_final,
  data = li6800.off.dat,
  prior = priors_mc,
  chains = 6,
  cores = 6,
  iter = 10000,
  warmup = 6000,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 12,
    metric = "diag_e"
  )
)

print(summary(fit_off_final))
pp_check(fit_off_final)
rstan::check_hmc_diagnostics(fit_off_final$fit)

##make the above plot prettier
li6800.off.plot <- li6800.off.dat
li6800.off.plot$fitted <- fitted(fit_off_final)[,1]
li6800.off.plot$fitted.fixed <- fitted(fit_off_final, re_formula = NA)[,1]

ggplot(li6800.off.plot, aes(x= tmin, y = gsw, group = Rep)) +
  geom_point() +
  facet_grid(Room_Temperature ~ Watering) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red")

ggplot(li6800.off.plot, aes(x= tmin, y = gsw, group = Rep)) +
  geom_point() +
  facet_grid(Room_Temperature ~ Watering) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red") +
  geom_point(aes(y=fitted.fixed), colour = "blue") +
  geom_line(aes(y=fitted.fixed), colour = "blue")

##I will set up an output script to get the parameters, etc. from these models
# ##calculate Slmax for the McAusland model
# Slmax.calc.mc <- function(coef){
#   ##pull out ceofficients, as they relate to parameters in McAusland 2016
#   gi <- coef$Rep[, "Estimate", "g0_Intercept"]
#   gf <- coef$Rep[, "Estimate", "g1_Intercept"]
#   k <- coef$Rep[, "Estimate", "k_Intercept"]
#   
#   (gf - gi)/(k*exp(1))
# }
# 
# Slmax.calc.mc(coef(fits_mc[[1]]))
# Slmax.calc.mc(coef(fits_mc[[2]]))
# Slmax.calc.mc(coef(fits_mc[[3]]))
# Slmax.calc.mc(coef(fits_mc[[4]]))

#'### Save data
save.image(saveFile)
