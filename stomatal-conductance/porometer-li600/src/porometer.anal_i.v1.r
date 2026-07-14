#'# Wheat stomatal kinetics (heat & drought)
#'## --- Ed Chaplin (USyd), *Apr--June 2025* ---
#'### Porometer DATA - ANALYSIS OF LIGHT ON PORTION
#'
#' * Analysis of gs increase, when light switched on
#'
#'### Load the libraries
library(knitr)
#knitr::spin("porometer.anal_i.v1.r")
library(ggplot2)
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, 
               "fig.width" = 7.5, "fig.height" = 10, "comment"=NA, fig.path = "figures/increase/")
library(dplyr)
##library(asreml)
library(brms) ##package used to fit bayesian regression model using Stan

#'### Set up constants and files
options(width=90)
load("../data/porometer.v1.Rdata")
saveFile <- "../data/porometer.anal_i.v1.Rdata"

#'### Fit a model to an individual Date to begin with, just for the light on portion
# li600.day1.dat <- subset(li600.dat, Date_Sys == min(Date_Sys) & Light == "On")
#head(li600.day1.dat)
# unique(li600.day1.dat$CEBLane)
# unique(li600.day1.dat$CEBPosition)
# unique(li600.day1.dat$CEBSide)
# unique(li600.day1.dat$CEBLocation)
# unique(li600.day1.dat$CEBBlock)
# unique(li600.day1.dat$CEBSubUnit)
# unique(li600.day1.dat$CEBCell)
# ggplot(li600.day1.dat, aes(x=CEBLane, y = CEBPosition)) +
#   geom_tile(fill = "white") +
#   geom_text(aes(label = CEBSubUnit))

##Fit a variance component model to the time course data
# li600.day1.varcomp <- asreml(gsw ~ (Device_Watering*Line) + 
#                                time.elapsed +
#                                Device_Watering:time.elapsed +
#                                Line:time.elapsed,
#                              random = ~ spl(time.elapsed),# +
#                                #Device_Watering:spl(time.elapsed),# +
#                                #Line:spl(time.elapsed),# +
#                                #Device_Watering:Line:spl(time.elapsed),
#                              #random = ~ CEBSubUnit + CEBSubUnit:CEBCell,
#                              residual = ~ units,
#                              data = li600.day1.dat)
# 
# li600.day1.varcomp <- update(li600.day1.varcomp)
# summary(li600.day1.varcomp)$varcomp
# wald(li600.day1.varcomp, denDF = "algebraic", ssType = "conditional")
# 
# plot(li600.day1.varcomp)
# 
# temp.asr <- asreml(gsw ~ (Device_Watering*Line) + 
#                      time.elapsed +
#                      Device_Watering:time.elapsed +
#                      Line:time.elapsed,
#                    #random = ~ spl(time.elapsed) +
#                    #  Device_Watering:spl(time.elapsed),# +
#                      #Line:spl(time.elapsed),# +
#                      #Device_Watering:Line:spl(time.elapsed),
#                    #random = ~ CEBSubUnit + CEBSubUnit:CEBCell,
#                    residual = ~ units,
#                    data = li600.day1.dat)
# temp.asr <- update(temp.asr)
# lrt(li600.day1.varcomp, temp.asr)

# li600.day1.varcomp.pred <- predict(li600.day1.varcomp, "Device_Watering:Line:time.elapsed",
#                                    pworkspace = "5gb", workspace = "5gb")#,
#                                    #levels = list(Line = "PBI16C009-0C-0N-0N-010N-14N"))
# head(li600.day1.varcomp.pred$pvals, 30)
# 
# ggplot(li600.day1.varcomp.pred$pvals, aes(x=time.elapsed, y = predicted.value)) +
#   geom_point(colour = "red") +
#   facet_grid(Device_Watering ~ Line) +
#   geom_point(data = li600.day1.dat, aes(x=time.elapsed, y = gsw))
# 
# li600.day1.varcomp.spl <- predict(li600.day1.varcomp, "Device_Watering:Line:time.elapsed",
#                                    only = "spl(time.elapsed)")
# ggplot(li600.day1.varcomp.spl$pvals, aes(x=time.elapsed, y = predicted.value)) +
#   geom_point(colour = "red") +
#   facet_grid(Device_Watering ~ Line)
# li600.day1.varcomp.linear <- predict(li600.day1.varcomp, "Device_Watering:Line:time.elapsed",
#                                   ignore = "spl(time.elapsed)")
# ggplot(li600.day1.varcomp.linear$pvals, aes(x=time.elapsed, y = predicted.value)) +
#   geom_point(colour = "red") +
#   facet_grid(Device_Watering ~ Line)

#'#### Stomatal conductance model - McAusland from McAusland et al 2016
#' Now let's try the Baesian approach that was fitted to the li-6800 data

# nlform <- bf(gsw ~ g0 + (g1 - g0) * exp(-exp(((lambda - time.elapsed)/k)+1)),
#                 g0 ~ 1 + Line*Device_Watering,
#                 g1 ~ 1 + Line*Device_Watering,
#                 k ~ 1 + Line*Device_Watering,
#                 lambda ~ 1 + Line*Device_Watering,
#                 nl = TRUE)
# 
# priors <- c(prior(normal(0.025, 0.02), nlpar = "g0"),
#                prior(normal(0.3, 0.2), nlpar = "g1"),
#                # prior(normal(0, 0.01), class = "sd", nlpar = "g0"),
#                # prior(normal(0, 0.05), class = "sd", nlpar = "g1"),
#                prior(normal(3, 2), nlpar = "k", lb = 0),
#                prior(normal(4, 2), nlpar = "lambda", lb = 0),
#                prior(student_t(3, 0, 0.01), class = "sigma"))

#'#### Initialise the model
#' Looks like to initialise the model, the dataframe is subset for the first plant in the first rep and for PAR == 1, which is the period where the light is switched on
# pick a small but valid dataset

# Configure and compile the model
# fit_init <- brm(
#   formula = nlform,
#   data = li600.day1.dat,
#   prior = priors,
#   chains = 4,
#   cores = 4,
#   iter = 3000,
#   warmup = 2000,
#   control = list(
#     adapt_delta = 0.99,
#     max_treedepth = 12,
#     metric = "diag_e"
#   )
# )
# 
# print(summary(fit_init))
# plot(fit_init)
# pairs(fit_init)
# mcmc_plot(fit_init, type = "trace")
# rstan::check_hmc_diagnostics(fit_init$fit)
# loo(fit_init)
# 
# pp_check(fit_init)  # compares density of the observed data with data generated from posterior predictive distribution
# pp_check(fit_init, type = "error_hist", ndraws = 12)
# pp_check(fit_init, type = "scatter_avg", ndraws = 1000)
# pp_check(fit_init, type = "stat_2d")
# pp_check(fit_init, type = "loo_ribbon")
# 
# plot(gsw~time.elapsed, li600.day1.dat, pch=16)
# points(fitted(fit_init)[,1]~time.elapsed, li600.day1.dat, pch=16, col="red")
# 
# ##make the above plot prettier
# li600.day1.plot <- li600.day1.dat
# li600.day1.plot$fitted <- fitted(fit_init)[,1]
# 
# ggplot(li600.day1.plot, aes(x= time.elapsed, y = gsw)) +
#   geom_point() +
#   facet_grid(Device_Watering ~ Line) +
#   geom_point(aes(y=fitted), colour = "red") +
#   geom_line(aes(y=fitted), colour = "red")
# 
# ##look at the effects from the model
# fixef(fit_init)
# ranef(fit_init)
# coef(fit_init)

# fit_init <- add_criterion(fit_init, "loo")
# fit_init2 <- add_criterion(fit_init2, "loo")
# 
# loo_compare(fit_init, fit_init2)

#'### Fit a model to an individual chamber
# li600.ambient.dat <- subset(li600.dat, Temperature == "Ambient" & Light == "On")
# #head(li600.ambient.dat)
# unique(li600.ambient.dat$CEBLane)
# unique(li600.ambient.dat$CEBPosition)
# unique(li600.ambient.dat$CEBSide)
# unique(li600.ambient.dat$CEBLocation)
# unique(li600.ambient.dat$CEBBlock)
# unique(li600.ambient.dat$CEBSubUnit)
# unique(li600.ambient.dat$CEBCell)
# ggplot(li600.ambient.dat, aes(x=CEBLane, y = CEBPosition)) +
#   geom_tile(fill = "white") +
#   geom_text(aes(label = CEBSubUnit))
# 
# ##Fit a variance component model to the time course data
# li600.ambient.varcomp <- asreml(gsw ~ (Device_Watering*Line) + 
#                                   time.elapsed +
#                                   Device_Watering:time.elapsed +
#                                   Line:time.elapsed +
#                                   Device_Watering:Line:time.elapsed,
#                                 random = ~ spl(time.elapsed) +
#                                 Device_Watering:spl(time.elapsed) +
#                                 Line:spl(time.elapsed) +
#                                 Device_Watering:Line:spl(time.elapsed) +
#                                 CEBSide + CEBSide:CEBLocation +
#                                   CEBSide:CEBLocation:CEBSubUnit +
#                                   CEBSide:CEBLocation:CEBSubUnit:CEBCell,
#                                 residual = ~ units,
#                                 data = li600.ambient.dat)
# 
# li600.ambient.varcomp <- update(li600.ambient.varcomp)
# summary(li600.ambient.varcomp)$varcomp
# wald(li600.ambient.varcomp, denDF = "algebraic", ssType = "conditional")
# 
# plot(li600.ambient.varcomp)

# temp.asr <- asreml(gsw ~ (Device_Watering*Line) + 
#                                   time.elapsed +
#                                   Device_Watering:time.elapsed +
#                                   Line:time.elapsed +
#                                   Device_Watering:Line:time.elapsed,
#                                 random = ~ spl(time.elapsed) +
#                                   Device_Watering:spl(time.elapsed) +
#                                   Line:spl(time.elapsed) +
#                                   #Device_Watering:Line:spl(time.elapsed) +
#                                   CEBSide + CEBSide:CEBLocation +
#                                   CEBSide:CEBLocation:CEBSubUnit +
#                                   CEBSide:CEBLocation:CEBSubUnit:CEBCell,
#                    residual = ~ units,
#                    data = li600.ambient.dat)
# temp.asr <- update(temp.asr)
# lrt(li600.ambient.varcomp, temp.asr)

# li600.ambient.varcomp.pred_d1 <- predict(li600.ambient.varcomp, 
#                                       "Device_Watering:Line:time.elapsed",
#                                    levels = list(Device_Watering = "D1_Unlimited"))
# head(li600.ambient.varcomp.pred_d1$pvals, 30)
# li600.ambient.varcomp.pred_d2 <- predict(li600.ambient.varcomp, 
#                                          "Device_Watering:Line:time.elapsed",
#                                          levels = list(Device_Watering = "D2_Limited"))
# head(li600.ambient.varcomp.pred_d2$pvals, 30)
# li600.ambient.varcomp.pred <- bind_rows(li600.ambient.varcomp.pred_d1$pvals,
#                                         li600.ambient.varcomp.pred_d2$pvals)
# 
# ggplot(li600.ambient.varcomp.pred, aes(x=time.elapsed, y = predicted.value)) +
#   geom_point(colour = "red") +
#   facet_grid(Device_Watering ~ Line) +
#   geom_point(data = li600.ambient.dat, aes(x=time.elapsed, y = gsw, 
#                                            colour = Replicate))

#'#### Stomatal conductance model - McAusland from McAusland et al 2016
#' Now let's try the Baesian approach that was fitted to the li-6800 data

# nlform_amb <- bf(gsw ~ g0 + (g1 - g0) * exp(-exp(((lambda - time.elapsed)/k)+1)),
#              g0 ~ 1 + Line*Device_Watering + 
#                (1 + Line*Device_Watering || Replicate),
#              g1 ~ 1 + Line*Device_Watering + 
#                (1 + Line*Device_Watering || Replicate),
#              k ~ 1 + Line*Device_Watering,
#              lambda ~ 1 + Line*Device_Watering,
#              nl = TRUE)
# 
# priors_amb <- c(prior(normal(0.025, 0.02), nlpar = "g0"),
#             prior(normal(0.4, 0.15), nlpar = "g1"),
#             prior(normal(3.5, 2), nlpar = "k", lb = 0),
#             prior(normal(5, 1), nlpar = "lambda", lb = 0),
#             prior(normal(0, 0.05), class = "sd", nlpar = "g0"),
#             prior(normal(0, 0.1), class = "sd", nlpar = "g1"),
#             prior(student_t(3, 0, 0.05), class = "sigma"))

#'#### Initialise the model
#' Fit the model just sampling from the priors to begin with

##check model using priors only
# fit_prior_amb <- brm(
#   formula = nlform_amb,
#   data = li600.ambient.dat,
#   prior = priors_amb,
#   sample_prior = "only",
#   chains = 1,
#   cores = 4
# )
# 
# pp_check(fit_prior_amb, ndraws = 25)
# summary(fit_prior_amb)


# Configure and compile the model
# fit_init <- brm(
#   formula = nlform,
#   data = li600.ambient.dat,
#   prior = priors,
#   chains = 4,
#   cores = 4,
#   iter = 5000,
#   warmup = 3000,
#   control = list(
#     adapt_delta = 0.99,
#     max_treedepth = 12,
#     metric = "diag_e"
#   )
# )

# print(summary(fit_init))
# plot(fit_init)
# pairs(fit_init)
# mcmc_plot(fit_init, type = "trace")
# rstan::check_hmc_diagnostics(fit_init$fit)
# loo(fit_init)
# 
# pp_check(fit_init)  # compares density of the observed data with data generated from posterior predictive distribution
# pp_check(fit_init, type = "error_hist", ndraws = 12)
# pp_check(fit_init, type = "scatter_avg", ndraws = 1000)
# pp_check(fit_init, type = "stat_2d")
# pp_check(fit_init, type = "loo_ribbon")

# plot(gsw~time.elapsed, li600.ambient.dat, pch=16)
# points(fitted(fit_init)[,1]~time.elapsed, li600.ambient.dat, pch=16, col="red")

##make the above plot prettier
# li600.ambient.plot <- li600.ambient.dat
# li600.ambient.plot$fitted <- fitted(fit_init)[,1]
# 
# ggplot(li600.ambient.plot, aes(x= time.elapsed, y = gsw, group = Replicate)) +
#   geom_point() +
#   facet_grid(Device_Watering ~ Line) +
#   geom_point(aes(y=fitted), colour = "red") +
#   geom_line(aes(y=fitted), colour = "red")
# 
# ##look at the effects from the model
# fixef(fit_init)
# ranef(fit_init)
# coef(fit_init)

# fit_init <- add_criterion(fit_init, "loo")
# fit_init2 <- add_criterion(fit_init2, "loo")
# 
# loo_compare(fit_init, fit_init2)

#'### Fit a model to the entire dataset, for the light on portion
li600.on.dat <- subset(li600.dat, Light == "On")
##using full dataset, figure out the earliest time at which the light was turned on for each day
li600.on.dat <- li600.on.dat %>% group_by(Date) %>% mutate(min.time = min(time.elapsed))
##calculate adjusted time by subtracting the earliest time the light was on each day
li600.on.dat$time.adj <- li600.on.dat$time.elapsed - li600.on.dat$min.time

#head(li600.dat)
# unique(li600.on.dat$CEBRoom)
# unique(li600.on.dat$CEBLane)
# unique(li600.on.dat$CEBPosition)
# unique(li600.on.dat$CEBSide)
# unique(li600.on.dat$CEBLocation)
# unique(li600.on.dat$CEBBlock)
# unique(li600.on.dat$CEBSubUnit)
# unique(li600.on.dat$CEBCell)
# ggplot(li600.on.dat, aes(x=CEBLane, y = CEBPosition)) +
#   geom_tile(fill = "white") +
#   geom_text(aes(label = CEBSide)) +
#   facet_wrap(~CEBRoom)

##Fit a variance component model to the time course data
# li600.on.varcomp <- asreml(gsw ~ (Room_Temperature*Device_Watering*Line) + 
#                           time.elapsed +
#                           Room_Temperature:time.elapsed +
#                           Device_Watering:time.elapsed +
#                           Line:time.elapsed +
#                           #Room_Temperature:Device_Watering:time.elapsed +
#                           Room_Temperature:Line:time.elapsed +
#                           Device_Watering:Line:time.elapsed,
#                         random = ~ spl(time.elapsed) +
#                           Room_Temperature:spl(time.elapsed) +
#                           Device_Watering:spl(time.elapsed) +
#                           Line:spl(time.elapsed) +
#                           #Room_Temperature:Device_Watering:spl(time.elapsed) +
#                           Room_Temperature:Line:spl(time.elapsed) +
#                           Device_Watering:Line:spl(time.elapsed) +
#                           CEBRoom:CEBSide + 
#                           CEBRoom:CEBSide:CEBLocation +
#                           CEBRoom:CEBSide:CEBLocation:CEBSubUnit +
#                           CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell,
#                         residual = ~ units,
#                         data = li600.on.dat,
#                         workspace = "5gb")
# 
# li600.on.varcomp <- update(li600.on.varcomp)
# summary(li600.on.varcomp)$varcomp
# wald(li600.on.varcomp, denDF = "algebraic", ssType = "conditional")
# 
# plot(li600.on.varcomp)

# temp.asr <- asreml(gsw ~ (Room_Temperature*Device_Watering*Line) + 
#                              time.elapsed +
#                              Room_Temperature:time.elapsed +
#                              Device_Watering:time.elapsed +
#                              Line:time.elapsed +
#                              #Room_Temperature:Device_Watering:time.elapsed +
#                              Room_Temperature:Line:time.elapsed +
#                              Device_Watering:Line:time.elapsed,
#                            random = ~ spl(time.elapsed) +
#                              Room_Temperature:spl(time.elapsed) +
#                              Device_Watering:spl(time.elapsed) +
#                              Line:spl(time.elapsed) +
#                              #Room_Temperature:Device_Watering:spl(time.elapsed) +
#                              Room_Temperature:Line:spl(time.elapsed) +
#                              Device_Watering:Line:spl(time.elapsed) +
#                              CEBRoom:CEBSide + 
#                              CEBRoom:CEBSide:CEBLocation +
#                              CEBRoom:CEBSide:CEBLocation:CEBSubUnit +
#                              CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell,
#                            residual = ~ units,
#                            data = li600.on.dat,
#                            workspace = "5gb")
# temp.asr <- update(temp.asr)
# lrt(li600.on.varcomp, temp.asr)

# li600.on.varcomp.pred_d1c8 <- predict(li600.on.varcomp, 
#                                          "Room_Temperature:Device_Watering:Line:time.elapsed",
#                                          levels = list(Device_Watering = "D1_Unlimited",
#                                                        Room_Temperature = "CER8_Ambient",
#                                                        time.elapsed = seq(38.5, 101.5, 1)),
#                                  workspace = "5gb",
#                                  pworkspace = "5gb")
# head(li600.on.varcomp.pred_d1c8$pvals, 30)
# li600.on.varcomp.pred_d1c7 <- predict(li600.on.varcomp, 
#                                    "Room_Temperature:Device_Watering:Line:time.elapsed",
#                                    levels = list(Device_Watering = "D1_Unlimited",
#                                                  Room_Temperature = "CER7_Heat",
#                                                  time.elapsed = seq(38.5, 101.5, 1)),
#                                    workspace = "5gb",
#                                    pworkspace = "5gb")
# head(li600.on.varcomp.pred_d1c7$pvals, 30)
# li600.on.varcomp.pred_d2c8 <- predict(li600.on.varcomp, 
#                                    "Room_Temperature:Device_Watering:Line:time.elapsed",
#                                    levels = list(Device_Watering = "D2_Limited",
#                                                  Room_Temperature = "CER8_Ambient",
#                                                  time.elapsed = seq(38.5, 101.5, 1)),
#                                    workspace = "5gb",
#                                    pworkspace = "5gb")
# head(li600.on.varcomp.pred_d2c8$pvals, 30)
# li600.on.varcomp.pred_d2c7 <- predict(li600.on.varcomp, 
#                                    "Room_Temperature:Device_Watering:Line:time.elapsed",
#                                    levels = list(Device_Watering = "D2_Limited",
#                                                  Room_Temperature = "CER7_Heat",
#                                                  time.elapsed = seq(38.5, 101.5, 1)),
#                                    workspace = "5gb",
#                                    pworkspace = "5gb")
# head(li600.on.varcomp.pred_d1c8$pvals, 30)
# li600.on.varcomp.pred <- bind_rows(li600.on.varcomp.pred_d1c8$pvals,
#                                    li600.on.varcomp.pred_d1c7$pvals,
#                                    li600.on.varcomp.pred_d2c8$pvals,
#                                    li600.on.varcomp.pred_d2c7$pvals)
# 
# ggplot(li600.on.varcomp.pred, aes(x=time.elapsed, y = predicted.value)) +
#   facet_grid(Room_Temperature + Device_Watering ~ Line) +
#   geom_point(data = li600.on.dat, aes(x=time.elapsed, y = gsw, 
#                                            colour = Replicate)) +
#   geom_point(colour = "red") 

#'#### Stomatal conductance model - McAusland from McAusland et al 2016
#' Now let's try the Baesian approach that was fitted to the li-6800 data

nlform <- bf(gsw ~ g0 + (g1 - g0) * exp(-exp(((lambda - time.adj)/k)+1)),
             g0 ~ 1 + Line*Device_Watering*Room_Temperature + 
               (1 | CEBRoom:CEBSide:CEBLocation +
                  CEBRoom:CEBSide:CEBLocation:CEBSubUnit + 
                  CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell),
             g1 ~ 1 + Line*Device_Watering*Room_Temperature + 
               (1 | CEBRoom:CEBSide:CEBLocation +
                  CEBRoom:CEBSide:CEBLocation:CEBSubUnit + 
                  CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell),
             k ~ 1 + Line*Device_Watering*Room_Temperature,
             lambda ~ 1 + Line*Device_Watering*Room_Temperature,
             nl = TRUE)

priors <- c(prior(normal(0, 0.02), nlpar = "g0"),
            prior(normal(0.4, 0.15), nlpar = "g1"),
            prior(normal(3.5, 2), nlpar = "k", lb = 0),
            prior(normal(5, 1), nlpar = "lambda", lb = 0),
            prior(normal(0, 0.05), class = "sd", nlpar = "g0"),
            prior(normal(0, 0.1), class = "sd", nlpar = "g1"),
            prior(student_t(3, 0, 0.05), class = "sigma"))

plot(0 + (0.4 - 0) * exp(-exp(((5 - seq(0, 60, 0.5))/3.5) + 1)) ~
       seq(0, 60, 0.5))

#'#### Initialise the model
# Fit the model
fit_on_max <- brm(
  formula = nlform,
  data = li600.on.dat,
  prior = priors,
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

print(summary(fit_on_max))
plot(fit_on_max, ask = FALSE)
# pairs(fit_on_max)
# mcmc_plot(fit_on_max, type = "trace")
rstan::check_hmc_diagnostics(fit_on_max$fit)
loo(fit_on_max)

pp_check(fit_on_max)  # compares density of the observed data with data generated from posterior predictive distribution
pp_check(fit_on_max, type = "error_hist", ndraws = 12)
pp_check(fit_on_max, type = "scatter_avg", ndraws = 1000)
pp_check(fit_on_max, type = "stat_2d")
pp_check(fit_on_max, type = "loo_ribbon")

plot(gsw~time.elapsed, li600.on.dat, pch=16)
points(fitted(fit_on_max)[,1]~time.elapsed, li600.on.dat, pch=16, col="red")

##make the above plot prettier
li600.on.plot <- li600.on.dat
li600.on.plot$fitted <- fitted(fit_on_max)[,1]
li600.on.plot$fitted.fixed <- fitted(fit_on_max, re_formula = NA)[,1]

ggplot(li600.on.plot, aes(x= time.adj, y = gsw, group = Replicate)) +
  geom_point() +
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red")

ggplot(li600.on.plot, aes(x= time.adj, y = gsw, group = Replicate)) +
  geom_point() +
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red") +
  geom_point(aes(y=fitted.fixed), colour = "blue") +
  geom_line(aes(y=fitted.fixed), colour = "blue")

##Testing contribution of effects in the model
# nlform.reduced <- bf(gsw ~ g0 + (g1 - g0) * exp(-exp(((lambda - time.adj)/k)+1)),
#              g0 ~ 1 + Line*Device_Watering*Room_Temperature + 
#                (1 | CEBRoom:CEBSide:CEBLocation +
#                   CEBRoom:CEBSide:CEBLocation:CEBSubUnit + 
#                   CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell),
#              g1 ~ 1 + Line + Device_Watering + Room_Temperature +
#                Line:Device_Watering +# Line:Room_Temperature +
#                Device_Watering:Room_Temperature + 
#                  (1 | CEBRoom:CEBSide:CEBLocation +
#                   CEBRoom:CEBSide:CEBLocation:CEBSubUnit + 
#                   CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell),
#              k ~ 1 + Line,# + Device_Watering,# +  Room_Temperature
#                 #Line:Room_Temperature,#Line:Device_Watering +
#              lambda ~ 1,# + Line,#  + Room_Temperature,# + Device_Watering 
#               #Line:Room_Temperature, # Line:Device_Watering
#              nl = TRUE)
# 
# # Fit the model
# fit_on_reduced <- brm(
#   formula = nlform.reduced,
#   data = li600.on.dat,
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
# print(summary(fit_on_reduced))
# loo(fit_on_reduced)
# loo(fit_on_current)
# rstan::check_hmc_diagnostics(fit_on_reduced$fit)
# plot(gsw~time.adj, li600.on.dat, pch=16)
# points(fitted(fit_on_reduced)[,1]~time.adj, li600.on.dat, pch=16, col="red")
# length(fixef(fit_on_reduced))
# 
# fit_on_current <- add_criterion(fit_on_current, "loo")
# fit_on_reduced <- add_criterion(fit_on_reduced, "loo")
# 
# loo_compare(fit_on_current, fit_on_reduced)
##THREE WAY INTERACTION
##Lambda: Drop Room_Temperature:Device_Watering:Line interaction
##k: Drop Room_Temperature:Device_Watering:Line interaction
##g0: Retain Room_Temperature:Device_Watering:Line interaction
##g1: Drop Room_Temperature:Device_Watering:Line interaction

##TWO WAY INTERACTION
##Lambda: Drop Room_Temperature:Device_Watering
##k: Drop Room_Temperature:Device_Watering
##g1: Retain Room_Temperature:Device_Watering
##lambda: Drop Line:Device_Watering
##lambda: Drop Line:Room_Temperature
##k: Drop Line:Device_Watering
##k: Drop Line:Room_Temperature
##g1: Drop Line:Room_Temperature
##g1: Retain Line:Device_Watering

##MAIN EFFECTS
##lambda: Drop Device_Watering
##lambda: Drop Room_Temperature
##lambda: Drop Line
##k: Drop Device_Watering
##k: Drop Room_Temperature
##k: Retain Line

# fit_on_current <- fit_on_reduced
# print(summary(fit_on_current))
# pp_check(fit_on_current) 

#'#### Fit final model
##Final model:
##g0: Room_Temperature:Device_Watering:Line interaction
##g1: Room_Temperature:Device_Watering + Line:Device_Watering
##lambda: 1
##k: Line

nlform_final <- bf(gsw ~ g0 + (g1 - g0) * exp(-exp(((lambda - time.adj)/k)+1)),
                     g0 ~ 1 + Line*Device_Watering*Room_Temperature + 
                       (1 | CEBRoom:CEBSide:CEBLocation +
                          CEBRoom:CEBSide:CEBLocation:CEBSubUnit + 
                          CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell),
                     g1 ~ 1 + Line + Device_Watering + Room_Temperature +
                       Line:Device_Watering + Device_Watering:Room_Temperature + 
                       (1 | CEBRoom:CEBSide:CEBLocation +
                          CEBRoom:CEBSide:CEBLocation:CEBSubUnit + 
                          CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell),
                     k ~ 1 + Line,
                     lambda ~ 1,
                     nl = TRUE)

# Fit the model
fit_on_final <- brm(
  formula = nlform_final,
  data = li600.on.dat,
  prior = priors,
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

print(summary(fit_on_final))
pp_check(fit_on_final)
rstan::check_hmc_diagnostics(fit_on_final$fit)
plot(fit_on_final, ask = FALSE)

##make the above plot prettier
li600.on.plot <- li600.on.dat
li600.on.plot$fitted <- fitted(fit_on_final)[,1]
li600.on.plot$fitted.fixed <- fitted(fit_on_final, re_formula = NA)[,1]

ggplot(li600.on.plot, aes(x= time.adj, y = gsw, group = Replicate)) +
  geom_point() +
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red")

ggplot(li600.on.plot, aes(x= time.adj, y = gsw, group = Replicate)) +
  geom_point() +
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red") +
  geom_point(aes(y=fitted.fixed), colour = "blue") +
  geom_line(aes(y=fitted.fixed), colour = "blue")

#'### Save data
save.image(saveFile)
