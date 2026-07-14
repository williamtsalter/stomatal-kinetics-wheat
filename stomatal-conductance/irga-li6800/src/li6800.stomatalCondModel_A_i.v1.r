#'# Wheat stomatal kinetics (heat & drought)
#'## --- Ed Chaplin (USyd), *Apr--June 2025* ---
#'### Li-6800 STOMATAL CONDUCTANCE MODEL FITTING TO A - LIGHT OFF TO LIGHT ON (INCREASE)
#'
#' * Using approach fit to gs data to now analyse A
#' * The model fit to A is slightly simpler than for gs as there is no lag parameter (lambda)
#'
#'### Load the libraries
library(knitr)
#knitr::spin("li6800.stomatalCondModel_A_i.v1.r")
library(ggplot2)
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, 
               "fig.width" = 7.5, "fig.height" = 10, "comment"=NA, fig.path = "figures/increaseA/")
library(dae)
library(dplyr)
library(growthPheno) # calcLagged
packageVersion("growthPheno")
library(brms) ##package used to fit bayesian regression model using Stan

#'### Set up constants and files
options(width=90)
load("../data/li6800.v1.Rdata")
saveFile <- "../data/stomatalCondModel_A_i.v1.Rdata"

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

#'### Explore the gs profiles of just the light on portion of the experiment
plotProfiles(li6800.dat %>% filter(Light == "On"), response = "A",
             individuals = "TWR", times = "tmin",
             facet.x = "Room_Temperature",
             facet.y = "Watering",
             y.title = "A",
             x.title = "minutes",
             title = "A",
             breaks.spacing.x = 5,
             alpha = 0.75)

#'### Following approach implemented for analysis of gs data
#' * Script begins by specifying a model for A and setting priors on the parameters

#'#### Model for A, based on McAusland et al 2016, where we don't have a lag (lambda) parameter
nlform_A <- bf(A ~ g0 + (g1 - g0) * exp(-exp(((-tmin)/k)+1)),
             g0 ~ 1 + Room_Temperature*Watering +
               (1 | Rep),
             g1 ~ 1 + Room_Temperature*Watering +
               (1 | Rep),
             k ~ 1+ Room_Temperature*Watering,
             nl = TRUE)

priors_A <- c(prior(normal(0.5, 0.1), nlpar = "g0"),
            prior(normal(25, 1), nlpar = "g1"),
            prior(normal(0, 0.5), class = "sd", nlpar = "g0"),
            prior(normal(0, 1), class = "sd", nlpar = "g1"),
            prior(normal(3, 0.5), nlpar = "k", lb = 0),
            prior(student_t(3, 0, 0.5), class = "sigma"))

plot(1.5 + (25 - 1.5) * exp(-exp(((- seq(0, 60, 0.5))/3) + 1)) ~
       seq(0, 60, 0.5))

#'#### Put finishing touches on dataframe
#' Subset for light on portion and remove th == 0 if it's present, as the light was still on at that point
li6800.on.dat <- subset(li6800.dat, Light == "On" & th != 0)
##in this case, we can align the time access to all start from zero
li6800.on.dat$th <- li6800.on.dat$th - li6800.on.dat$th[1]
li6800.on.dat$tmin <- li6800.on.dat$tmin - li6800.on.dat$tmin[1]

plotProfiles(li6800.on.dat, response = "A",
             individuals = "TWR", times = "tmin",
             y.title = "A",
             x.title = "minutes",
             title = "A",
             breaks.spacing.x = 5,
             alpha = 0.75)

#'#### Fit the model
# Configure and compile the model - McAusland
fit_on_A_max <- brm(
  formula = nlform_A,
  data = li6800.on.dat,
  prior = priors_A,
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

print(summary(fit_on_A_max))
plot(fit_on_A_max, ask = FALSE)
pairs(fit_on_A_max)
mcmc_plot(fit_on_A_max, type = "trace")
rstan::check_hmc_diagnostics(fit_on_A_max$fit)
loo(fit_on_A_max)

pp_check(fit_on_A_max)  # compares density of the observed data with data generated from posterior predictive distribution
pp_check(fit_on_A_max, type = "error_hist", ndraws = 12)
pp_check(fit_on_A_max, type = "scatter_avg", ndraws = 1000)
pp_check(fit_on_A_max, type = "stat_2d")
pp_check(fit_on_A_max, type = "loo_ribbon")

plot(A~tmin, li6800.on.dat, pch=16)
points(fitted(fit_on_A_max)[,1]~tmin, li6800.on.dat, pch=16, col="red")

##make the above plot prettier
li6800.on.plot <- li6800.on.dat
li6800.on.plot$fitted <- fitted(fit_on_A_max)[,1]
li6800.on.plot$fitted.fixed <- fitted(fit_on_A_max, re_formula = NA)[,1]

ggplot(li6800.on.plot, aes(x= tmin, y = A, group = Rep)) +
  geom_point() +
  facet_grid(Room_Temperature ~ Watering) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red")

ggplot(li6800.on.plot, aes(x= tmin, y = A, group = Rep)) +
  geom_point() +
  facet_grid(Room_Temperature ~ Watering) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red") +
  geom_point(aes(y=fitted.fixed), colour = "blue") +
  geom_line(aes(y=fitted.fixed), colour = "blue")

#'#### Testing contribution of effects in the model
# fit_on_A_current <- fit_on_A_max
# 
# nlform_A_reduced <- bf(A ~ g0 + (g1 - g0) * exp(-exp(((-tmin)/k)+1)),
#                g0 ~ 1 + Room_Temperature+Watering +
#                  (1 | Rep),
#                g1 ~ 1 + Room_Temperature*Watering +
#                  (1 | Rep),
#                k ~ 1+  Watering, #Room_Temperature +
#                nl = TRUE)
# 
# # Fit the model
# fit_on_A_reduced <- brm(
#   formula = nlform_A_reduced,
#   data = li6800.on.dat,
#   prior = priors_A,
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
# 
# print(summary(fit_on_A_reduced))
# loo(fit_on_A_reduced)
# loo(fit_on_A_current)
# plot(A~tmin, li6800.on.dat, pch=16)
# points(fitted(fit_on_A_reduced)[,1]~tmin, li6800.on.dat, pch=16, col="red")
# coef(fit_on_A_reduced)
# 
# fit_on_A_current <- add_criterion(fit_on_A_current, "loo")
# fit_on_A_reduced <- add_criterion(fit_on_A_reduced, "loo")
# 
# loo_compare(fit_on_A_current, fit_on_A_reduced)
##TWO WAY INTERACTIONS
##k: Drop Room_Temperature:Watering interaction
##g1: Retain Room_Temperature:Watering interaction
##g0: Retain Room_Temperature:Watering interaction

##MAIN EFFECT
##k: Retain Watering main effect
##k: Retain Room_Temperature main effect
# fit_on_A_current <- fit_on_A_reduced
# print(summary(fit_on_A_current))

#'#### Fit final model
##Final model:
##g0: Room_Temperature*Watering
##g1: Room_Temperature*Watering
##k: Room_Temperature + Watering

nlform_A_final <- bf(A ~ g0 + (g1 - g0) * exp(-exp(((-tmin)/k)+1)),
               g0 ~ 1 + Room_Temperature*Watering +
                 (1 | Rep),
               g1 ~ 1 + Room_Temperature*Watering +
                 (1 | Rep),
               k ~ 1+  Room_Temperature + Watering,
               nl = TRUE)

# Fit the model
fit_on_A_final <- brm(
  formula = nlform_A_final,
  data = li6800.on.dat,
  prior = priors_A,
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

print(summary(fit_on_A_final))
pp_check(fit_on_A_final)
rstan::check_hmc_diagnostics(fit_on_A_final$fit)

##make the above plot prettier
li6800.on.plot <- li6800.on.dat
li6800.on.plot$fitted <- fitted(fit_on_A_final)[,1]
li6800.on.plot$fitted.fixed <- fitted(fit_on_A_final, re_formula = NA)[,1]

ggplot(li6800.on.plot, aes(x= tmin, y = A, group = Rep)) +
  geom_point() +
  facet_grid(Room_Temperature ~ Watering) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red")

ggplot(li6800.on.plot, aes(x= tmin, y = A, group = Rep)) +
  geom_point() +
  facet_grid(Room_Temperature ~ Watering) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red") +
  geom_point(aes(y=fitted.fixed), colour = "blue") +
  geom_line(aes(y=fitted.fixed), colour = "blue")

#'### Save data
save.image(saveFile)
