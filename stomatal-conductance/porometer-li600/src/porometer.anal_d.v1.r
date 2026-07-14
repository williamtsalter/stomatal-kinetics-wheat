#'# Wheat stomatal kinetics (heat & drought)
#'## --- Ed Chaplin (USyd), *Apr--June 2025* ---
#'### Porometer DATA - ANALYSIS OF LIGHT OFF PORTION
#'
#' * Analysis of gs decrease, when light switched off for the second time
#'
#'### Load the libraries
library(knitr)
#knitr::spin("porometer.anal_d.v1.r")
library(ggplot2)
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, 
               "fig.width" = 7.5, "fig.height" = 10, "comment"=NA, fig.path = "figures/increase/")
library(dplyr)
library(asreml)
library(brms) ##package used to fit bayesian regression model using Stan

#'### Set up constants and files
options(width=90)
load("../data/porometer.v1.Rdata")
saveFile <- "../data/porometer.anal_d.v1.Rdata"

#'### Fit a model to the entire dataset, for the light off portion, following light on
##find the minimum time when the light was on
min(li600.dat$time.elapsed[li600.dat$Light == "On"]) ##38.46667
max(li600.dat$time.elapsed[li600.dat$Light == "On"]) ##101.5667
hist(li600.dat$time.elapsed[li600.dat$Light == "Off" &
                             li600.dat$time.elapsed > 38.46667 &
                             li600.dat$time.elapsed < 101.5667])

li600.off.dat <- subset(li600.dat, Light == "Off" &
                          time.elapsed > 60)

##using full dataset, figure out the earliest time at which the light was turned off for each day
li600.off.dat <- li600.off.dat %>% group_by(Date) %>% mutate(min.time = min(time.elapsed))
##calculate adjusted time by subtracting the earliest time the light was on each day
li600.off.dat$time.adj <- li600.off.dat$time.elapsed - li600.off.dat$min.time

#head(li600.dat)
# unique(li600.off.dat$CEBRoom)
# unique(li600.off.dat$CEBLane)
# unique(li600.off.dat$CEBPosition)
# unique(li600.off.dat$CEBSide)
# unique(li600.off.dat$CEBLocation)
# unique(li600.off.dat$CEBBlock)
# unique(li600.off.dat$CEBSubUnit)
# unique(li600.off.dat$CEBCell)
# ggplot(li600.off.dat, aes(x=CEBLane, y = CEBPosition)) +
#   geom_tile(fill = "white") +
#   geom_text(aes(label = CEBSide)) +
#   facet_wrap(~CEBRoom)

##Fit a variance component model to the time course data
# li600.off.varcomp <- asreml(gsw ~ (Room_Temperature*Device_Watering*Line) + 
#                           time.elapsed +
#                           Room_Temperature:time.elapsed +
#                           Device_Watering:time.elapsed +
#                           Line:time.elapsed +
#                           Room_Temperature:Device_Watering:time.elapsed +
#                           Room_Temperature:Line:time.elapsed +
#                           Device_Watering:Line:time.elapsed,# +
#                           #Room_Temperature:Device_Watering:Line:time.elapsed,
#                         random = ~ spl(time.elapsed) +
#                           Room_Temperature:spl(time.elapsed) +
#                           Device_Watering:spl(time.elapsed) +
#                           Line:spl(time.elapsed) +
#                           Room_Temperature:Device_Watering:spl(time.elapsed) +
#                           Room_Temperature:Line:spl(time.elapsed) +
#                           Device_Watering:Line:spl(time.elapsed) +
#                           #Room_Temperature:Device_Watering:Line:spl(time.elapsed) +
#                           CEBRoom:CEBSide + 
#                           CEBRoom:CEBSide:CEBLocation +
#                           CEBRoom:CEBSide:CEBLocation:CEBSubUnit +
#                           CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell,
#                         residual = ~ units,
#                         data = li600.off.dat,
#                         workspace = "5gb")
# 
# li600.off.varcomp <- update(li600.off.varcomp)
# summary(li600.off.varcomp)$varcomp
# wald(li600.off.varcomp, denDF = "algebraic", ssType = "conditional")

# plot(li600.off.varcomp)

# temp.asr <- asreml(gsw ~ (Room_Temperature*Device_Watering*Line) + 
#                               time.elapsed +
#                               Room_Temperature:time.elapsed +
#                               Device_Watering:time.elapsed +
#                               Line:time.elapsed +
#                               Room_Temperature:Device_Watering:time.elapsed +
#                               Room_Temperature:Line:time.elapsed +
#                               Device_Watering:Line:time.elapsed,# +
#                               #Room_Temperature:Device_Watering:Line:time.elapsed,
#                             random = ~ spl(time.elapsed) +
#                               Room_Temperature:spl(time.elapsed) +
#                               Device_Watering:spl(time.elapsed) +
#                               Line:spl(time.elapsed) +
#                               Room_Temperature:Device_Watering:spl(time.elapsed) +
#                               Room_Temperature:Line:spl(time.elapsed) +
#                               Device_Watering:Line:spl(time.elapsed) +
#                               #Room_Temperature:Device_Watering:Line:spl(time.elapsed) +
#                               CEBRoom:CEBSide + 
#                               CEBRoom:CEBSide:CEBLocation +
#                               CEBRoom:CEBSide:CEBLocation:CEBSubUnit +
#                               CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell,
#                             residual = ~ units,
#                             data = li600.off.dat,
#                             workspace = "5gb")
# temp.asr <- update(temp.asr)
# lrt(li600.off.varcomp, temp.asr)

# li600.off.varcomp.pred_d1c8 <- predict(li600.off.varcomp, 
#                                          "Room_Temperature:Device_Watering:Line:time.elapsed",
#                                          levels = list(Device_Watering = "D1_Unlimited",
#                                                        Room_Temperature = "CER8_Ambient",
#                                                        time.elapsed = seq(100, 165, 1)),
#                                  workspace = "5gb",
#                                  pworkspace = "5gb")
# head(li600.off.varcomp.pred_d1c8$pvals, 30)
# li600.off.varcomp.pred_d1c7 <- predict(li600.off.varcomp, 
#                                    "Room_Temperature:Device_Watering:Line:time.elapsed",
#                                    levels = list(Device_Watering = "D1_Unlimited",
#                                                  Room_Temperature = "CER7_Heat",
#                                                  time.elapsed = seq(100, 165, 1)),
#                                    workspace = "5gb",
#                                    pworkspace = "5gb")
# head(li600.off.varcomp.pred_d1c7$pvals, 30)
# li600.off.varcomp.pred_d2c8 <- predict(li600.off.varcomp, 
#                                    "Room_Temperature:Device_Watering:Line:time.elapsed",
#                                    levels = list(Device_Watering = "D2_Limited",
#                                                  Room_Temperature = "CER8_Ambient",
#                                                  time.elapsed = seq(100, 165, 1)),
#                                    workspace = "5gb",
#                                    pworkspace = "5gb")
# head(li600.off.varcomp.pred_d2c8$pvals, 30)
# li600.off.varcomp.pred_d2c7 <- predict(li600.off.varcomp, 
#                                    "Room_Temperature:Device_Watering:Line:time.elapsed",
#                                    levels = list(Device_Watering = "D2_Limited",
#                                                  Room_Temperature = "CER7_Heat",
#                                                  time.elapsed = seq(100, 165, 1)),
#                                    workspace = "5gb",
#                                    pworkspace = "5gb")
# head(li600.off.varcomp.pred_d1c8$pvals, 30)
# li600.off.varcomp.pred <- bind_rows(li600.off.varcomp.pred_d1c8$pvals,
#                                    li600.off.varcomp.pred_d1c7$pvals,
#                                    li600.off.varcomp.pred_d2c8$pvals,
#                                    li600.off.varcomp.pred_d2c7$pvals)
# 
# ggplot(li600.off.varcomp.pred, aes(x=time.elapsed, y = predicted.value)) +
#   facet_grid(Room_Temperature + Device_Watering ~ Line) +
#   geom_point(data = li600.off.dat, aes(x=time.elapsed, y = gsw, 
#                                            colour = Replicate)) +
#   geom_point(colour = "red") 

#'#### Stomatal conductance model - McAusland from McAusland et al 2016
#' Now let's try the Baesian approach that was fitted to the li-6800 data

nlform <- bf(gsw ~ g0 + (g1 - g0) * exp(-exp(((lambda - time.adj)/k)+1)),
             g0 ~ 1 + Line*Device_Watering*Room_Temperature + 
               (1 | CEBRoom:CEBSide + CEBRoom:CEBSide:CEBLocation +
                  CEBRoom:CEBSide:CEBLocation:CEBSubUnit + 
                  CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell),
             g1 ~ 1 + Line*Device_Watering*Room_Temperature + 
               (1 | CEBRoom:CEBSide + CEBRoom:CEBSide:CEBLocation +
                  CEBRoom:CEBSide:CEBLocation:CEBSubUnit + 
                  CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell),
             k ~ 1 + Line*Device_Watering*Room_Temperature,
             lambda ~ 1 + Line*Device_Watering*Room_Temperature,
             nl = TRUE)

priors <- c(prior(normal(0.2, 0.5), nlpar = "g0"),
            prior(normal(0, 0.02), nlpar = "g1"),
            prior(normal(3, 0.1), nlpar = "k", lb = 0),
            prior(normal(0.5, 0.25), nlpar = "lambda", lb = 0),
            prior(normal(0, 0.15), class = "sd", nlpar = "g0"),
            prior(normal(0, 0.01), class = "sd", nlpar = "g1"),
            prior(student_t(3, 0, 0.005), class = "sigma"))

plot(0.2 + (0 - 0.2) * exp(-exp(((0.5 - seq(0, 60, 0.5))/3) + 1)) ~
       seq(0, 60, 0.5))


# Fit the model
fit_off_max <- brm(
  formula = nlform,
  data = li600.off.dat,
  prior = priors,
  chains = 6,
  cores = 6,
  iter = 3000,
  warmup = 2000,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 12,
    metric = "diag_e"
  )
)

print(summary(fit_off_max))
plot(fit_off_max, ask = FALSE)
#pairs(fit_off_max)
#mcmc_plot(fit_off_max, type = "trace")
rstan::check_hmc_diagnostics(fit_off_max$fit)
loo(fit_off_max)

pp_check(fit_off_max)  # compares density of the observed data with data generated from posterior predictive distribution
pp_check(fit_off_max, type = "error_hist", ndraws = 12)
pp_check(fit_off_max, type = "scatter_avg", ndraws = 1000)
pp_check(fit_off_max, type = "stat_2d")
pp_check(fit_off_max, type = "loo_ribbon")

plot(gsw~time.adj, li600.off.dat, pch=16)
points(fitted(fit_off_max)[,1]~time.adj, li600.off.dat, pch=16, col="red")

##make the above plot prettier
li600.off.plot <- li600.off.dat
li600.off.plot$fitted <- fitted(fit_off_max)[,1]
li600.off.plot$fitted.fixed <- fitted(fit_off_max, re_formula = NA)[,1]

ggplot(li600.off.plot, aes(x= time.adj, y = gsw, group = Replicate)) +
  geom_point() +
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red")

ggplot(li600.off.plot, aes(x= time.adj, y = gsw, group = Replicate)) +
  geom_point() +
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red") +
  geom_point(aes(y=fitted.fixed), colour = "blue") +
  geom_line(aes(y=fitted.fixed), colour = "blue")

##Testing contribution of effects in the model
# nlform.reduced <- bf(gsw ~ g0 + (g1 - g0) * exp(-exp(((lambda - time.adj)/k)+1)),
#                      g0 ~ 1 + Line*Device_Watering*Room_Temperature +
#                        (1 | CEBRoom:CEBSide:CEBLocation +
#                           CEBRoom:CEBSide:CEBLocation:CEBSubUnit + 
#                           CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell),
#                      g1 ~ 1 + Line + Device_Watering + #Room_Temperature +
#                        Line:Device_Watering + 
#                         #Line:Room_Temperature + # +
#                        #Device_Watering:Room_Temperature + 
#                        (1 | CEBRoom:CEBSide:CEBLocation +
#                           CEBRoom:CEBSide:CEBLocation:CEBSubUnit + 
#                           CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell),
#                      k ~ 1  + Room_Temperature,#+ Line +  Device_Watering,# +
#                        #Line:Device_Watering,# + Line:Room_Temperature,# +
#                        #Device_Watering:Room_Temperature,
#                      lambda ~ 1 + Device_Watering,#+ Line + Room_Temperature,# +
#                        #Line:Device_Watering,# + Line:Room_Temperature,
#                      nl = TRUE)
# 
# # Fit the model
# fit_off_reduced <- brm(
#   formula = nlform.reduced,
#   data = li600.off.dat,
#   prior = priors,
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
# rstan::check_hmc_diagnostics(fit_off_reduced$fit)
# plot(gsw~time.adj, li600.on.dat, pch=16)
# points(fitted(fit_off_reduced)[,1]~time.adj, li600.on.dat, pch=16, col="red")
# length(fixef(fit_off_reduced))
# 
# fit_off_current <- add_criterion(fit_off_current, "loo")
# fit_off_reduced <- add_criterion(fit_off_reduced, "loo")
# 
# loo_compare(fit_off_current, fit_off_reduced)
##THREE WAY INTERACTION
##Lambda: Drop Room_Temperature:Device_Watering:Line interaction
##k: Drop Room_Temperature:Device_Watering:Line interaction 
##g0: Retain Room_Temperature:Device_Watering:Line interaction
##g1: Drop Room_Temperature:Device_Watering:Line interaction 

##TWO WAY INTERACTION
##lambda: Drop Room_Temperature:Device_Watering
##k: Drop Room_Temperature:Device_Watering
##g1: Drop Room_Temperature:Device_Watering
##lambda: Drop Line:Room_Temperature
##lambda: Drop Line:Device_Watering
##k: Drop Line:Room_Temperature
##k: Drop Line:Device_Watering
##g1: Retain Line:Device_Watering
##g1: Drop Line:Room_Temperature

##MAIN EFFECTS
##lambda: Drop Room_Temperature
##lambda: Retain Device_Watering
##lambda: Drop Line
##k: Retain Room_Temperature
##k: Drop Device_Watering
##k: Drop Line
##g1: Drop Room_Temperature
# fit_off_current <- fit_off_reduced
# print(summary(fit_off_current))

#'#### Fit final model
##Final model:
##g0: Room_Temperature:Device_Watering:Line interaction
##g1: Line:Device_Watering
##lambda: Device_Watering
##k: Room_Temperature

nlform_final <- bf(gsw ~ g0 + (g1 - g0) * exp(-exp(((lambda - time.adj)/k)+1)),
                   g0 ~ 1 + Line*Device_Watering*Room_Temperature + 
                     (1 | CEBRoom:CEBSide:CEBLocation +
                        CEBRoom:CEBSide:CEBLocation:CEBSubUnit + 
                        CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell),
                   g1 ~ 1 + Line + Device_Watering + Line:Device_Watering + 
                     (1 | CEBRoom:CEBSide:CEBLocation +
                        CEBRoom:CEBSide:CEBLocation:CEBSubUnit + 
                        CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell),
                   k ~ 1 + Room_Temperature,
                   lambda ~ 1+ Device_Watering,
                   nl = TRUE)

# Fit the model
fit_off_final <- brm(
  formula = nlform_final,
  data = li600.off.dat,
  prior = priors,
  chains = 6,
  cores = 6,
  iter = 10000,
  warmup = 6000,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15,
    metric = "diag_e"
  )
)

print(summary(fit_off_final))
pp_check(fit_off_final)
rstan::check_hmc_diagnostics(fit_off_final$fit)
plot(fit_off_final, ask = FALSE)

##make the above plot prettier
li600.off.plot <- li600.off.dat
li600.off.plot$fitted <- fitted(fit_off_final)[,1]
li600.off.plot$fitted.fixed <- fitted(fit_off_final, re_formula = NA)[,1]

ggplot(li600.off.plot, aes(x= time.adj, y = gsw, group = Replicate)) +
  geom_point() +
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red")

ggplot(li600.off.plot, aes(x= time.adj, y = gsw, group = Replicate)) +
  geom_point() +
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red") +
  geom_point(aes(y=fitted.fixed), colour = "blue") +
  geom_line(aes(y=fitted.fixed), colour = "blue")

#'### Save data
save.image(saveFile)
