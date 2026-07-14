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
#'### Load the libraries
library(knitr)
#knitr::spin("li6800.stomatalCondModel_d.v1.r")
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
saveFile <- "../data/stomatalCondModel_d.v1.Rdata"

#'### Define TW to be unique for each treatment
li6800.dat$TW <- paste(li6800.dat$Temperature, li6800.dat$Watering, sep="_")
li6800.dat$TW <- factor(li6800.dat$TW)

#'### Define th and tmin, to be time in hours and minutes respectively
li6800.dat$th <- li6800.dat$xElapsed/3600
li6800.dat$tmin <- li6800.dat$xElapsed/60

#'### Explore the gs profiles of just the second light off portion of the experiment
#' Most responses show a slight increase before the decreasing response begins
plotProfiles(li6800.dat %>% filter(Light == "Off" & th > 2), response = "gsw",
             individuals = "TWR", times = "tmin",
             facet.x = "Temperature",
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

#'#### Stomatal conductance model - CDWeibull from Woning et al 2026
nlform <- bf(gsw ~ g1 + (g0 - g1) * exp(-((th / tau)^theta)),
             g0 ~ 1 + (1 | Rep),
             g1 ~ 1 + (1 | Rep),
             tau ~ 1,
             theta ~ 1,
             nl = TRUE)

#'#### Stomatal conductance model - McAusland from McAusland et al 2016
nlform_mc <- bf(gsw ~ g0 + (g1 - g0) * exp(-exp(((lambda - tmin)/k)+1)),
             g0 ~ 1 + (1 | Rep),
             g1 ~ 1 + (1 | Rep),
             k ~ 1,
             lambda ~ 1,
             nl = TRUE)

priors <- c(prior(normal(0.1, 0.2), nlpar = "g0"),
            prior(normal(0.01, 0.02), nlpar = "g1"),
            prior(normal(0, 0.05), class = "sd", nlpar = "g0"),
            prior(normal(0, 0.01), class = "sd", nlpar = "g1"),
            prior(normal(0.25, 1), nlpar = "tau", lb = 0),
            prior(normal(2, 2), nlpar = "theta", lb = 0),
            prior(student_t(3, 0, 0.01), class = "sigma"))

priors_mc <- c(prior(normal(0.1, 0.2), nlpar = "g0"),
            prior(normal(0.01, 0.02), nlpar = "g1"),
            prior(normal(0, 0.05), class = "sd", nlpar = "g0"),
            prior(normal(0, 0.01), class = "sd", nlpar = "g1"),
            prior(normal(3, 2), nlpar = "k", lb = 0),
            prior(normal(0.5, 2), nlpar = "lambda", lb = 0),
            prior(student_t(3, 0, 0.01), class = "sigma"))

#'#### Initialise the model
#' Looks like to initialise the model, the dataframe is subset for the first treatment and in this case, for the second period where the light is turned off
# pick a small but valid dataset
dd_init <- subset(li6800.dat, Light == "Off" & TW == li6800.dat$TW[1] & th != 0 & th > 2)
##did all reps experience the light off treatment at the same time?
dd_init %>% group_by(Rep) %>% summarise(min(th)) %>% as.data.frame()
##No, not quite the same

##Align each rep to make sure time 0 is the same
dd_init <- dd_init %>% group_by(Rep) %>% mutate(th = th - th[1]) %>% as.data.frame()
min(dd_init$th)
dd_init <- dd_init %>% group_by(Rep) %>% mutate(tmin = tmin - tmin[1]) %>% as.data.frame()
dd_init %>% group_by(Rep) %>% summarise(min(tmin)) %>% as.data.frame()

##remove the first observation from the analysis - it's still impacted by the change in lighting
dd_init <- subset(dd_init, th != 0)
head(dd_init)

plotProfiles(dd_init, response = "gsw",
             individuals = "TWR", times = "tmin",
             y.title = "gsw",
             x.title = "hours",
             title = "gsw",
             breaks.spacing.x = 5,
             alpha = 0.75)

# Configure and compile the model
fit_init <- brm(
  formula = nlform,
  data = dd_init,
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 3000,
  warmup = 2000,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 12,
    metric = "diag_e"
  )
)

print(summary(fit_init))
plot(fit_init)
pairs(fit_init)
mcmc_plot(fit_init, type = "trace")
rstan::check_hmc_diagnostics(fit_init$fit)
loo(fit_init)

pp_check(fit_init)  # compares density of the observed data with data generated from posterior predictive distribution
pp_check(fit_init, type = "error_hist", ndraws = 12)
pp_check(fit_init, type = "scatter_avg", ndraws = 1000)
pp_check(fit_init, type = "stat_2d")
pp_check(fit_init, type = "loo_ribbon")

plot(gsw~th, dd_init, pch=16, main = unique(dd_init$TW))
points(fitted(fit_init)[,1]~th, dd_init, pch=16, col="red")

# Configure and compile the model - McAusland
fit_init_mc <- brm(
  formula = nlform_mc,
  data = dd_init,
  prior = priors_mc,
  chains = 4,
  cores = 4,
  iter = 3000,
  warmup = 2000,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 12,
    metric = "diag_e"
  )
)

print(summary(fit_init_mc))
plot(fit_init_mc)
pairs(fit_init_mc)
mcmc_plot(fit_init_mc, type = "trace")
rstan::check_hmc_diagnostics(fit_init_mc$fit)
loo(fit_init_mc)

pp_check(fit_init_mc)  # compares density of the observed data with data generated from posterior predictive distribution
pp_check(fit_init_mc, type = "error_hist", ndraws = 12)
pp_check(fit_init_mc, type = "scatter_avg", ndraws = 1000)
pp_check(fit_init_mc, type = "stat_2d")
pp_check(fit_init_mc, type = "loo_ribbon")

plot(gsw~tmin, dd_init, pch=16, main = unique(dd_init$TW))
points(fitted(fit_init_mc)[,1]~tmin, dd_init, pch=16, col="red")

##look at the effects from the model
fixef(fit_init_mc)
ranef(fit_init_mc)
coef(fit_init_mc)

#'#### Write a function to process each individual treatment
process <- function(d1, fit_template) {
  dd <- subset(d1, Light == "Off" & th != 0 & !is.na(gsw) & th > 2)
  ##align each replicate to have the same t = 0 value
  dd <- dd %>% group_by(Rep) %>% mutate(th = th - th[1]) %>% as.data.frame()
  dd <- dd %>% group_by(Rep) %>% mutate(tmin = tmin - tmin[1]) %>% as.data.frame()
  
  ##remove the first time point, as seems to be impacted by change of lighting still
  dd <- subset(dd, th != 0)
  
  fit <- update(
    fit_template,
    newdata = dd,
    recompile = NULL,
    refresh = 0
  )
  
  # Check diagnostics, Important!
  print(summary(fit))
  rstan::check_hmc_diagnostics(fit$fit)
  #plot(fit)
  #pp_check(fit)
  # plot(gsw~th, dd, pch=16, main = unique(dd$TW))
  # points(fitted(fit)[,1]~th, dd, pch=16, col="red")
  
  return(fit)
}

#'#### Implement loop
# Main loop that select each plant per species
fits <- list()

for (i in unique(li6800.dat$TW)) {
  d1 <- subset(li6800.dat, TW == i)
  
  cat("TW:", i, "\n")
  
  fits[[i]] <- process(d1, fit_init)
}

#'#### Explore model output
pp_check(fits[[1]])  # compares density of the observed data with data generated from posterior predictive distribution
pp_check(fits[[2]])  # compares density of the observed data with data generated from posterior predictive distribution
pp_check(fits[[3]])  # compares density of the observed data with data generated from posterior predictive distribution
pp_check(fits[[4]])  # compares density of the observed data with data generated from posterior predictive distribution

##look at the effects from the model
coef(fits[[1]])
coef(fits[[2]])
coef(fits[[3]])
coef(fits[[4]])

#'#### Implement loop - McAusland
# Main loop that select each plant per species
fits_mc <- list()

for (i in unique(li6800.dat$TW)) {
  d1 <- subset(li6800.dat, TW == i)
  
  cat("TW:", i, "\n")
  
  fits_mc[[i]] <- process(d1, fit_init_mc)
}

#'#### Explore model output
pp_check(fits_mc[[1]])  # compares density of the observed data with data generated from posterior predictive distribution
pp_check(fits_mc[[2]])  # compares density of the observed data with data generated from posterior predictive distribution
pp_check(fits_mc[[3]])  # compares density of the observed data with data generated from posterior predictive distribution
pp_check(fits_mc[[4]])  # compares density of the observed data with data generated from posterior predictive distribution

##look at the effects from the model
coef(fits_mc[[1]])
coef(fits_mc[[2]])
coef(fits_mc[[3]])
coef(fits_mc[[4]])

##calculate Slmax
Slmax.calc <- function(coef){
  ##pull out ceofficients, as they relate to parameters in Woning et al 2026
  gi <- coef$Rep[, "Estimate", "g0_Intercept"]
  gf <- coef$Rep[, "Estimate", "g1_Intercept"]
  lambda_c <- coef$Rep[, "Estimate", "theta_Intercept"]
  tau_c <- coef$Rep[, "Estimate", "tau_Intercept"]
  
  numerator_calc <- -1*(gi - gf) * lambda_c * exp((lambda_c - 1)/lambda_c)
  denominator_calc <- ((lambda_c - 1)/lambda_c)^(1/lambda_c) * tau_c
  
  numerator_calc/denominator_calc
}

Slmax.calc(coef(fits[[1]]))
Slmax.calc(coef(fits[[2]]))
Slmax.calc(coef(fits[[3]]))
Slmax.calc(coef(fits[[4]]))
##these values seem off - go thing we're more interested in the McAusland model!

##calculate Slmax for the McAusland model
Slmax.calc.mc <- function(coef){
  ##pull out ceofficients, as they relate to parameters in McAusland 2016
  gi <- coef$Rep[, "Estimate", "g0_Intercept"]
  gf <- coef$Rep[, "Estimate", "g1_Intercept"]
  k <- coef$Rep[, "Estimate", "k_Intercept"]
  
  (gf - gi)/(k*exp(1))
}

Slmax.calc.mc(coef(fits_mc[[1]]))
Slmax.calc.mc(coef(fits_mc[[2]]))
Slmax.calc.mc(coef(fits_mc[[3]]))
Slmax.calc.mc(coef(fits_mc[[4]]))

#'### Save data
save.image(saveFile)