---
output:
  pdf_document: default
  html_document: default
---
# Wheat stomatal kinetics (heat & drought)
## --- Ed Chaplin (USyd), *Apr--June 2025* ---
### Porometer DATA - ANALYSIS OF LIGHT ON PORTION

* Analysis of gs increase, when light switched on

### Load the libraries


``` r
library(knitr)
#knitr::spin("porometer.anal_i.v1.r")
library(ggplot2)
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, 
               "fig.width" = 7.5, "fig.height" = 10, "comment"=NA, fig.path = "figures/increase/")
library(dplyr)
##library(asreml)
library(brms) ##package used to fit bayesian regression model using Stan
```

### Set up constants and files


``` r
options(width=90)
load("../data/porometer.v1.Rdata")
saveFile <- "../data/porometer.anal_i.v1.Rdata"
```

### Fit a model to an individual Date to begin with, just for the light on portion


``` r
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
```

#### Stomatal conductance model - McAusland from McAusland et al 2016
Now let's try the Baesian approach that was fitted to the li-6800 data


``` r
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
```

#### Initialise the model
Looks like to initialise the model, the dataframe is subset for the first plant in the first rep and for PAR == 1, which is the period where the light is switched on


``` r
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
```

### Fit a model to an individual chamber


``` r
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
```

#### Stomatal conductance model - McAusland from McAusland et al 2016
Now let's try the Baesian approach that was fitted to the li-6800 data


``` r
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
```

#### Initialise the model
Fit the model just sampling from the priors to begin with


``` r
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
```

### Fit a model to the entire dataset, for the light on portion


``` r
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
```

#### Stomatal conductance model - McAusland from McAusland et al 2016
Now let's try the Baesian approach that was fitted to the li-6800 data


``` r
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
```



![plot of chunk unnamed-chunk-10](figures/increase/unnamed-chunk-10-1.png)

#### Initialise the model


``` r
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
```

```
Compiling Stan program...
```

```
Start sampling
```

``` r
print(summary(fit_on_max))
```

```
 Family: gaussian 
  Links: mu = identity 
Formula: gsw ~ g0 + (g1 - g0) * exp(-exp(((lambda - time.adj)/k) + 1)) 
         g0 ~ 1 + Line * Device_Watering * Room_Temperature + (1 | CEBRoom:CEBSide:CEBLocation + CEBRoom:CEBSide:CEBLocation:CEBSubUnit + CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell)
         g1 ~ 1 + Line * Device_Watering * Room_Temperature + (1 | CEBRoom:CEBSide:CEBLocation + CEBRoom:CEBSide:CEBLocation:CEBSubUnit + CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell)
         k ~ 1 + Line * Device_Watering * Room_Temperature
         lambda ~ 1 + Line * Device_Watering * Room_Temperature
   Data: li600.on.dat (Number of observations: 1682) 
  Draws: 6 chains, each with iter = 10000; warmup = 6000; thin = 1;
         total post-warmup draws = 24000

Multilevel Hyperparameters:
~CEBRoom:CEBSide:CEBLocation (Number of levels: 8) 
                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(g0_Intercept)     0.02      0.02     0.00     0.07 1.00     7418    13208
sd(g1_Intercept)     0.05      0.04     0.00     0.14 1.00     9766    13366

~CEBRoom:CEBSide:CEBLocation:CEBSubUnit (Number of levels: 16) 
                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(g0_Intercept)     0.07      0.02     0.04     0.11 1.00     9538    11961
sd(g1_Intercept)     0.14      0.04     0.08     0.22 1.00     9644    14718

~CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell (Number of levels: 159) 
                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(g0_Intercept)     0.04      0.01     0.03     0.06 1.00     7281     9726
sd(g1_Intercept)     0.12      0.01     0.09     0.15 1.00     5957    11029

Regression Coefficients:
                                                                                              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
g0_Intercept                                                                                      0.02      0.02    -0.01     0.05 1.00    20967    18215
g0_LinePBI17N002M0CM0NM0NM010NM16N                                                               -0.00      0.02    -0.03     0.03 1.00    29851    19447
g0_LinePBI17N006M0NM0NM89N                                                                        0.01      0.02    -0.02     0.04 1.00    33546    19089
g0_LinePBI17N007M0CM0NM0NM010NM8N                                                                 0.01      0.02    -0.03     0.04 1.00    31308    18309
g0_LinePBI17N013M0CM0NM0NM010NM8N                                                                 0.03      0.02     0.00     0.07 1.00    32436    19648
g0_LinePBI17N016M0CM0NM0NM010NM27N                                                                0.00      0.02    -0.03     0.03 1.00    32461    20068
g0_LinePBI19N007M0NM3N                                                                           -0.01      0.02    -0.04     0.02 1.00    31846    19160
g0_LinePBIC15034M0CM24NM010NM1NM0N                                                                0.01      0.02    -0.02     0.04 1.00    32757    19025
g0_LineROCKSTAR                                                                                   0.00      0.02    -0.03     0.03 1.00    32676    19420
g0_LineSUNCHASER                                                                                  0.02      0.02    -0.01     0.05 1.00    29263    17810
g0_Device_WateringD2_Limited                                                                      0.04      0.02     0.00     0.08 1.00    17215    16746
g0_Room_TemperatureCER8_Ambient                                                                   0.02      0.02    -0.02     0.05 1.00    26351    19160
g0_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited                                     -0.02      0.02    -0.06     0.01 1.00    31210    18654
g0_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited                                              0.01      0.02    -0.03     0.04 1.00    37734    19053
g0_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                      -0.02      0.02    -0.05     0.02 1.00    34625    18065
g0_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                       0.03      0.02    -0.00     0.06 1.00    32997    19862
g0_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited                                      0.00      0.02    -0.03     0.04 1.00    35785    18332
g0_LinePBI19N007M0NM3N:Device_WateringD2_Limited                                                 -0.00      0.02    -0.04     0.03 1.00    36458    19475
g0_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited                                      0.02      0.02    -0.01     0.05 1.00    32320    19000
g0_LineROCKSTAR:Device_WateringD2_Limited                                                        -0.00      0.02    -0.04     0.03 1.00    35342    18811
g0_LineSUNCHASER:Device_WateringD2_Limited                                                        0.03      0.02    -0.01     0.06 1.00    34329    18356
g0_LinePBI17N002M0CM0NM0NM010NM16N:Room_TemperatureCER8_Ambient                                  -0.00      0.02    -0.04     0.03 1.00    33268    19203
g0_LinePBI17N006M0NM0NM89N:Room_TemperatureCER8_Ambient                                           0.02      0.02    -0.02     0.05 1.00    34725    19538
g0_LinePBI17N007M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                    0.00      0.02    -0.03     0.04 1.00    35540    19117
g0_LinePBI17N013M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                    0.02      0.02    -0.01     0.06 1.00    30948    18867
g0_LinePBI17N016M0CM0NM0NM010NM27N:Room_TemperatureCER8_Ambient                                   0.01      0.02    -0.03     0.04 1.00    33074    18222
g0_LinePBI19N007M0NM3N:Room_TemperatureCER8_Ambient                                              -0.01      0.02    -0.04     0.03 1.00    33390    18587
g0_LinePBIC15034M0CM24NM010NM1NM0N:Room_TemperatureCER8_Ambient                                  -0.00      0.02    -0.03     0.03 1.00    33001    19302
g0_LineROCKSTAR:Room_TemperatureCER8_Ambient                                                      0.01      0.02    -0.03     0.04 1.00    34653    18735
g0_LineSUNCHASER:Room_TemperatureCER8_Ambient                                                    -0.02      0.02    -0.05     0.02 1.00    32791    18272
g0_Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                                         0.02      0.02    -0.02     0.06 1.00    27883    19320
g0_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient        -0.02      0.02    -0.05     0.02 1.00    36122    17829
g0_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                 0.01      0.02    -0.02     0.05 1.00    35775    17393
g0_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient         -0.00      0.02    -0.04     0.03 1.00    39634    18709
g0_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          0.03      0.02    -0.01     0.06 1.00    34695    18347
g0_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient         0.01      0.02    -0.03     0.04 1.00    40340    19133
g0_LinePBI19N007M0NM3N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                     0.00      0.02    -0.03     0.04 1.00    35919    17998
g0_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient        -0.01      0.02    -0.04     0.03 1.00    37992    19291
g0_LineROCKSTAR:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                           -0.00      0.02    -0.04     0.03 1.00    35229    16857
g0_LineSUNCHASER:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                          -0.01      0.02    -0.04     0.03 1.00    36388    19077
g1_Intercept                                                                                      0.20      0.07     0.06     0.33 1.00    10270    13910
g1_LinePBI17N002M0CM0NM0NM010NM16N                                                                0.03      0.06    -0.08     0.15 1.00    11199    14776
g1_LinePBI17N006M0NM0NM89N                                                                        0.11      0.06     0.00     0.23 1.00     9381    12977
g1_LinePBI17N007M0CM0NM0NM010NM8N                                                                 0.08      0.06    -0.04     0.19 1.00    10281    12173
g1_LinePBI17N013M0CM0NM0NM010NM8N                                                                 0.05      0.06    -0.07     0.16 1.00    10351    14525
g1_LinePBI17N016M0CM0NM0NM010NM27N                                                               -0.03      0.06    -0.14     0.09 1.00     9390    12854
g1_LinePBI19N007M0NM3N                                                                            0.15      0.06     0.04     0.27 1.00    12073    13948
g1_LinePBIC15034M0CM24NM010NM1NM0N                                                                0.16      0.06     0.04     0.28 1.00    12320    15926
g1_LineROCKSTAR                                                                                   0.11      0.06    -0.01     0.23 1.00    11151    14834
g1_LineSUNCHASER                                                                                  0.14      0.06     0.02     0.26 1.00    11288    15149
g1_Device_WateringD2_Limited                                                                     -0.13      0.08    -0.28     0.04 1.00    12232    14491
g1_Room_TemperatureCER8_Ambient                                                                   0.03      0.09    -0.13     0.20 1.00    12860    15092
g1_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited                                      0.15      0.08     0.00     0.30 1.00    14167    15637
g1_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited                                              0.25      0.08     0.10     0.40 1.00    13403    14747
g1_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                      -0.07      0.08    -0.21     0.09 1.00    12588    15214
g1_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                       0.09      0.08    -0.06     0.26 1.00    14186    15788
g1_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited                                      0.26      0.07     0.11     0.41 1.00    14432    16260
g1_LinePBI19N007M0NM3N:Device_WateringD2_Limited                                                  0.20      0.08     0.05     0.35 1.00    13931    16564
g1_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited                                      0.03      0.08    -0.12     0.19 1.00    14933    16729
g1_LineROCKSTAR:Device_WateringD2_Limited                                                         0.20      0.08     0.05     0.35 1.00    13145    15404
g1_LineSUNCHASER:Device_WateringD2_Limited                                                        0.18      0.08     0.03     0.34 1.00    13986    15648
g1_LinePBI17N002M0CM0NM0NM010NM16N:Room_TemperatureCER8_Ambient                                   0.06      0.08    -0.09     0.22 1.00    13558    16451
g1_LinePBI17N006M0NM0NM89N:Room_TemperatureCER8_Ambient                                           0.07      0.08    -0.08     0.23 1.00    12623    15036
g1_LinePBI17N007M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                    0.05      0.08    -0.10     0.20 1.00    11953    15690
g1_LinePBI17N013M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                    0.12      0.08    -0.03     0.28 1.00    13501    16285
g1_LinePBI17N016M0CM0NM0NM010NM27N:Room_TemperatureCER8_Ambient                                   0.14      0.08    -0.01     0.30 1.00    13781    17080
g1_LinePBI19N007M0NM3N:Room_TemperatureCER8_Ambient                                               0.11      0.08    -0.05     0.27 1.00    15006    17395
g1_LinePBIC15034M0CM24NM010NM1NM0N:Room_TemperatureCER8_Ambient                                   0.10      0.08    -0.06     0.27 1.00    15382    16653
g1_LineROCKSTAR:Room_TemperatureCER8_Ambient                                                      0.14      0.08    -0.02     0.30 1.00    12325    15078
g1_LineSUNCHASER:Room_TemperatureCER8_Ambient                                                     0.09      0.08    -0.07     0.25 1.00    13962    17036
g1_Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                                        -0.02      0.10    -0.21     0.18 1.00    17497    16465
g1_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient         0.15      0.11    -0.06     0.37 1.00    24687    19604
g1_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                 0.15      0.10    -0.05     0.36 1.00    20137    18009
g1_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          0.16      0.10    -0.05     0.36 1.00    21399    18623
g1_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          0.24      0.11     0.03     0.46 1.00    23985    19843
g1_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient         0.21      0.11    -0.00     0.42 1.00    22360    19288
g1_LinePBI19N007M0NM3N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                     0.14      0.11    -0.07     0.36 1.00    22127    19015
g1_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient         0.08      0.11    -0.14     0.30 1.00    25400    18237
g1_LineROCKSTAR:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                            0.28      0.11     0.06     0.50 1.00    21210    18618
g1_LineSUNCHASER:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                           0.21      0.11    -0.01     0.43 1.00    24069    19574
k_Intercept                                                                                       3.41      0.46     2.53     4.33 1.00    24751    18365
k_LinePBI17N002M0CM0NM0NM010NM16N                                                                 2.46      1.52     0.20     5.90 1.00    20467    11504
k_LinePBI17N006M0NM0NM89N                                                                         0.69      0.56     0.02     2.09 1.00    22350    13640
k_LinePBI17N007M0CM0NM0NM010NM8N                                                                  0.37      0.35     0.01     1.29 1.00    23059    12554
k_LinePBI17N013M0CM0NM0NM010NM8N                                                                  1.11      0.93     0.04     3.52 1.00    23501    14545
k_LinePBI17N016M0CM0NM0NM010NM27N                                                                 1.13      0.90     0.05     3.39 1.00    23215    13662
k_LinePBI19N007M0NM3N                                                                             1.67      0.97     0.14     3.80 1.00    20700    11339
k_LinePBIC15034M0CM24NM010NM1NM0N                                                                 7.19      1.32     4.68     9.90 1.00    28483    16824
k_LineROCKSTAR                                                                                    2.99      1.18     0.82     5.43 1.00    21364    11788
k_LineSUNCHASER                                                                                   3.97      1.37     1.40     6.80 1.00    18624    11188
k_Device_WateringD2_Limited                                                                       0.16      0.17     0.00     0.61 1.00    25933    12965
k_Room_TemperatureCER8_Ambient                                                                    3.46      1.46     1.00     6.80 1.00    10141    10072
k_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited                                       3.16      1.77     0.25     6.92 1.00    14712     7990
k_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited                                               0.55      0.53     0.02     1.95 1.00    24152    13011
k_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                        4.18      1.87     0.73     8.01 1.00    18170     9594
k_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                        4.32      1.89     0.77     8.09 1.00    18902    10100
k_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited                                       1.31      1.11     0.04     4.18 1.00    22165    13869
k_LinePBI19N007M0NM3N:Device_WateringD2_Limited                                                   0.55      0.55     0.01     2.00 1.00    23274    13265
k_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited                                       3.81      1.93     0.42     7.74 1.00    16485     8830
k_LineROCKSTAR:Device_WateringD2_Limited                                                          2.22      1.50     0.12     5.67 1.00    19905    12252
k_LineSUNCHASER:Device_WateringD2_Limited                                                         1.96      1.54     0.07     5.68 1.00    17567    11994
k_LinePBI17N002M0CM0NM0NM010NM16N:Room_TemperatureCER8_Ambient                                    3.89      1.88     0.54     7.73 1.00    17870     8850
k_LinePBI17N006M0NM0NM89N:Room_TemperatureCER8_Ambient                                            2.72      1.77     0.15     6.68 1.00    16209    11605
k_LinePBI17N007M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                     1.02      1.03     0.02     3.86 1.00    19353    13402
k_LinePBI17N013M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                     3.45      1.72     0.43     7.05 1.00    20026     9756
k_LinePBI17N016M0CM0NM0NM010NM27N:Room_TemperatureCER8_Ambient                                    3.77      1.83     0.51     7.60 1.00    16732     8252
k_LinePBI19N007M0NM3N:Room_TemperatureCER8_Ambient                                                5.15      1.68     1.93     8.51 1.00    24232    12215
k_LinePBIC15034M0CM24NM010NM1NM0N:Room_TemperatureCER8_Ambient                                    2.40      1.62     0.12     6.07 1.00    19559    11746
k_LineROCKSTAR:Room_TemperatureCER8_Ambient                                                       2.62      1.72     0.13     6.48 1.00    19963    11594
k_LineSUNCHASER:Room_TemperatureCER8_Ambient                                                      4.37      1.84     0.90     8.09 1.00    20391    10006
k_Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                                          8.18      2.38     3.03    12.57 1.00     9973     7662
k_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          4.26      1.87     0.82     8.06 1.00    19796     9897
k_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                  4.37      1.91     0.78     8.19 1.00    17073     9838
k_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient           4.01      1.84     0.65     7.75 1.00    17961     9278
k_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient           3.91      1.88     0.54     7.78 1.00    17721     8444
k_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          3.96      1.87     0.56     7.79 1.00    17920     9516
k_LinePBI19N007M0NM3N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                      4.16      1.89     0.68     7.95 1.00    18479     9511
k_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          4.16      1.86     0.70     7.91 1.00    15624     7668
k_LineROCKSTAR:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                             4.39      1.91     0.79     8.20 1.00    15751     8154
k_LineSUNCHASER:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                            4.13      1.87     0.65     7.94 1.00    18421     8947
lambda_Intercept                                                                                  0.46      0.37     0.02     1.38 1.00    23181    14370
lambda_LinePBI17N002M0CM0NM0NM010NM16N                                                            3.93      0.92     2.13     5.71 1.00    24493    11826
lambda_LinePBI17N006M0NM0NM89N                                                                    3.28      0.77     1.75     4.77 1.00    29116    15168
lambda_LinePBI17N007M0CM0NM0NM010NM8N                                                             2.39      0.78     0.85     3.95 1.00    21630    11715
lambda_LinePBI17N013M0CM0NM0NM010NM8N                                                             3.50      0.86     1.82     5.20 1.00    28622    13622
lambda_LinePBI17N016M0CM0NM0NM010NM27N                                                            3.84      0.88     2.11     5.55 1.00    29810    14286
lambda_LinePBI19N007M0NM3N                                                                        4.53      0.86     2.80     6.19 1.00    26092    15565
lambda_LinePBIC15034M0CM24NM010NM1NM0N                                                            4.97      0.93     3.13     6.82 1.00    33165    16237
lambda_LineROCKSTAR                                                                               3.52      0.88     1.78     5.22 1.00    21940    11270
lambda_LineSUNCHASER                                                                              4.25      0.89     2.48     5.98 1.00    29887    15045
lambda_Device_WateringD2_Limited                                                                  0.45      0.38     0.02     1.43 1.00    20467    12836
lambda_Room_TemperatureCER8_Ambient                                                               0.77      0.55     0.03     2.06 1.00    19831    12655
lambda_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited                                  4.83      0.99     2.89     6.76 1.00    36164    13603
lambda_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited                                          3.37      0.88     1.67     5.11 1.00    25917    12721
lambda_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                   4.99      1.00     3.05     6.96 1.00    36805    14369
lambda_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                   4.99      1.01     3.02     6.97 1.00    37057    14876
lambda_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited                                  4.04      0.96     2.17     5.91 1.00    28240    13390
lambda_LinePBI19N007M0NM3N:Device_WateringD2_Limited                                              3.48      0.89     1.75     5.23 1.00    28791    14775
lambda_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited                                  4.70      1.01     2.71     6.68 1.00    31902    14217
lambda_LineROCKSTAR:Device_WateringD2_Limited                                                     4.03      0.96     2.16     5.93 1.00    30187    12720
lambda_LineSUNCHASER:Device_WateringD2_Limited                                                    4.23      1.01     2.27     6.22 1.00    28419    14117
lambda_LinePBI17N002M0CM0NM0NM010NM16N:Room_TemperatureCER8_Ambient                               4.60      1.00     2.63     6.56 1.00    33796    13904
 [ reached 'max' / getOption("max.print") -- omitted 18 rows ]

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sigma     0.10      0.00     0.09     0.10 1.00    16427    16951

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```

``` r
plot(fit_on_max, ask = FALSE)
```



![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-1.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-2.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-3.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-4.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-5.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-6.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-7.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-8.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-9.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-10.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-11.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-12.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-13.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-14.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-15.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-16.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-17.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-18.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-19.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-20.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-21.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-22.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-23.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-24.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-25.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-26.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-27.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-28.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-29.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-30.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-31.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-32.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-33.png)

![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-34.png)

``` r
# pairs(fit_on_max)
# mcmc_plot(fit_on_max, type = "trace")
rstan::check_hmc_diagnostics(fit_on_max$fit)
```

```

Divergences:
```

```
0 of 24000 iterations ended with a divergence.
```

```

Tree depth:
```

```
0 of 24000 iterations saturated the maximum tree depth of 12.
```

```

Energy:
```

```
E-BFMI indicated no pathological behavior.
```

``` r
loo(fit_on_max)
```

```

Computed from 24000 by 1682 log-likelihood matrix.

         Estimate   SE
elpd_loo   1406.0 39.7
p_loo       210.1  9.1
looic     -2812.0 79.4
------
MCSE of elpd_loo is 0.2.
MCSE and ESS estimates assume MCMC draws (r_eff in [0.5, 1.6]).

All Pareto k estimates are good (k < 0.7).
See help('pareto-k-diagnostic') for details.
```

``` r
pp_check(fit_on_max)  # compares density of the observed data with data generated from posterior predictive distribution
```

```
Using 10 posterior draws for ppc type 'dens_overlay' by default.
```



![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-35.png)

``` r
pp_check(fit_on_max, type = "error_hist", ndraws = 12)
```

```
`stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```



![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-36.png)

``` r
pp_check(fit_on_max, type = "scatter_avg", ndraws = 1000)
```



![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-37.png)

``` r
pp_check(fit_on_max, type = "stat_2d")
```

```
Using all posterior draws for ppc type 'stat_2d' by default.
```

```
Note: in most cases the default test statistic 'mean' is too weak to detect anything of interest.
```



![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-38.png)

``` r
pp_check(fit_on_max, type = "loo_ribbon")
```

```
Using all posterior draws for ppc type 'loo_ribbon' by default.
```



![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-39.png)

``` r
plot(gsw~time.elapsed, li600.on.dat, pch=16)
points(fitted(fit_on_max)[,1]~time.elapsed, li600.on.dat, pch=16, col="red")
```



![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-40.png)

``` r
##make the above plot prettier
li600.on.plot <- li600.on.dat
li600.on.plot$fitted <- fitted(fit_on_max)[,1]
li600.on.plot$fitted.fixed <- fitted(fit_on_max, re_formula = NA)[,1]

ggplot(li600.on.plot, aes(x= time.adj, y = gsw, group = Replicate)) +
  geom_point() +
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red")
```



![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-41.png)

``` r
ggplot(li600.on.plot, aes(x= time.adj, y = gsw, group = Replicate)) +
  geom_point() +
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red") +
  geom_point(aes(y=fitted.fixed), colour = "blue") +
  geom_line(aes(y=fitted.fixed), colour = "blue")
```



![plot of chunk unnamed-chunk-11](figures/increase/unnamed-chunk-11-42.png)

``` r
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
```

#### Fit final model


``` r
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
```

```
Compiling Stan program...
```

```
Start sampling
```

``` r
print(summary(fit_on_final))
```

```
 Family: gaussian 
  Links: mu = identity 
Formula: gsw ~ g0 + (g1 - g0) * exp(-exp(((lambda - time.adj)/k) + 1)) 
         g0 ~ 1 + Line * Device_Watering * Room_Temperature + (1 | CEBRoom:CEBSide:CEBLocation + CEBRoom:CEBSide:CEBLocation:CEBSubUnit + CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell)
         g1 ~ 1 + Line + Device_Watering + Room_Temperature + Line:Device_Watering + Device_Watering:Room_Temperature + (1 | CEBRoom:CEBSide:CEBLocation + CEBRoom:CEBSide:CEBLocation:CEBSubUnit + CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell)
         k ~ 1 + Line
         lambda ~ 1
   Data: li600.on.dat (Number of observations: 1682) 
  Draws: 6 chains, each with iter = 10000; warmup = 6000; thin = 1;
         total post-warmup draws = 24000

Multilevel Hyperparameters:
~CEBRoom:CEBSide:CEBLocation (Number of levels: 8) 
                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(g0_Intercept)     0.02      0.02     0.00     0.06 1.00     6963    11013
sd(g1_Intercept)     0.03      0.03     0.00     0.10 1.00     9505    13168

~CEBRoom:CEBSide:CEBLocation:CEBSubUnit (Number of levels: 16) 
                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(g0_Intercept)     0.05      0.01     0.02     0.08 1.00     7495    10063
sd(g1_Intercept)     0.09      0.03     0.06     0.15 1.00     6857    11296

~CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell (Number of levels: 159) 
                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(g0_Intercept)     0.01      0.01     0.00     0.02 1.00     9583    11234
sd(g1_Intercept)     0.09      0.01     0.08     0.11 1.00     7107    13229

Regression Coefficients:
                                                                                          Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
g0_Intercept                                                                                  0.01      0.01    -0.02     0.03 1.00    21374    18377
g0_LinePBI17N002M0CM0NM0NM010NM16N                                                           -0.00      0.01    -0.03     0.03 1.00    33359    19700
g0_LinePBI17N006M0NM0NM89N                                                                   -0.00      0.01    -0.03     0.02 1.00    34740    17875
g0_LinePBI17N007M0CM0NM0NM010NM8N                                                             0.02      0.01    -0.01     0.05 1.00    34771    20434
g0_LinePBI17N013M0CM0NM0NM010NM8N                                                             0.02      0.01    -0.01     0.05 1.00    35655    19229
g0_LinePBI17N016M0CM0NM0NM010NM27N                                                           -0.01      0.01    -0.04     0.02 1.00    32398    18345
g0_LinePBI19N007M0NM3N                                                                       -0.01      0.01    -0.04     0.01 1.00    32076    19680
g0_LinePBIC15034M0CM24NM010NM1NM0N                                                           -0.01      0.01    -0.03     0.02 1.00    32546    19343
g0_LineROCKSTAR                                                                              -0.01      0.01    -0.03     0.02 1.00    31874    18439
g0_LineSUNCHASER                                                                              0.01      0.01    -0.02     0.04 1.00    32576    18932
g0_Device_WateringD2_Limited                                                                  0.03      0.02    -0.00     0.07 1.00    14681    15747
g0_Room_TemperatureCER8_Ambient                                                              -0.01      0.02    -0.04     0.02 1.00    23363    17530
g0_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited                                 -0.02      0.02    -0.05     0.02 1.00    38994    18962
g0_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited                                         -0.01      0.02    -0.04     0.03 1.00    42742    18068
g0_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                  -0.01      0.02    -0.04     0.03 1.00    40578    18754
g0_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                   0.01      0.02    -0.03     0.04 1.00    40674    17797
g0_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited                                 -0.00      0.02    -0.04     0.03 1.00    39445    17987
g0_LinePBI19N007M0NM3N:Device_WateringD2_Limited                                             -0.01      0.02    -0.04     0.03 1.00    36801    17772
g0_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited                                  0.04      0.02     0.00     0.07 1.00    37478    20138
g0_LineROCKSTAR:Device_WateringD2_Limited                                                    -0.01      0.02    -0.04     0.03 1.00    38905    18454
g0_LineSUNCHASER:Device_WateringD2_Limited                                                    0.03      0.02     0.00     0.07 1.00    35010    18616
g0_LinePBI17N002M0CM0NM0NM010NM16N:Room_TemperatureCER8_Ambient                               0.00      0.02    -0.03     0.04 1.00    37793    19011
g0_LinePBI17N006M0NM0NM89N:Room_TemperatureCER8_Ambient                                       0.00      0.02    -0.03     0.04 1.00    40907    18830
g0_LinePBI17N007M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                               -0.01      0.02    -0.04     0.02 1.00    38269    19380
g0_LinePBI17N013M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                0.01      0.02    -0.02     0.04 1.00    37654    19274
g0_LinePBI17N016M0CM0NM0NM010NM27N:Room_TemperatureCER8_Ambient                              -0.00      0.02    -0.03     0.03 1.00    37029    18582
g0_LinePBI19N007M0NM3N:Room_TemperatureCER8_Ambient                                          -0.01      0.02    -0.04     0.02 1.00    33867    19088
g0_LinePBIC15034M0CM24NM010NM1NM0N:Room_TemperatureCER8_Ambient                               0.01      0.02    -0.02     0.04 1.00    34421    19926
g0_LineROCKSTAR:Room_TemperatureCER8_Ambient                                                  0.00      0.02    -0.03     0.03 1.00    38740    18773
g0_LineSUNCHASER:Room_TemperatureCER8_Ambient                                                -0.02      0.02    -0.05     0.01 1.00    34969    18410
g0_Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                                     0.00      0.02    -0.03     0.04 1.00    27535    18239
g0_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient    -0.01      0.02    -0.04     0.03 1.00    40641    17788
g0_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient             0.00      0.02    -0.03     0.04 1.00    41087    18352
g0_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient     -0.00      0.02    -0.04     0.03 1.00    41760    18833
g0_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient      0.01      0.02    -0.03     0.05 1.00    38157    17839
g0_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient    -0.00      0.02    -0.04     0.03 1.00    40082    17755
g0_LinePBI19N007M0NM3N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                -0.00      0.02    -0.04     0.03 1.00    40014    17878
g0_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient     0.01      0.02    -0.03     0.05 1.00    40476    17933
g0_LineROCKSTAR:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                       -0.00      0.02    -0.04     0.03 1.00    40652    18475
g0_LineSUNCHASER:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                      -0.01      0.02    -0.04     0.03 1.00    38947    17616
g1_Intercept                                                                                  0.23      0.05     0.12     0.33 1.00     8691    12643
g1_LinePBI17N002M0CM0NM0NM010NM16N                                                            0.01      0.04    -0.07     0.09 1.00     7446    12161
g1_LinePBI17N006M0NM0NM89N                                                                    0.09      0.04     0.02     0.18 1.00     7400    10983
g1_LinePBI17N007M0CM0NM0NM010NM8N                                                             0.08      0.04    -0.00     0.16 1.00     6840    11112
g1_LinePBI17N013M0CM0NM0NM010NM8N                                                             0.05      0.04    -0.03     0.13 1.00     7300    12307
g1_LinePBI17N016M0CM0NM0NM010NM27N                                                           -0.02      0.04    -0.10     0.07 1.00     7440    12301
g1_LinePBI19N007M0NM3N                                                                        0.13      0.04     0.05     0.21 1.00     7969    12614
g1_LinePBIC15034M0CM24NM010NM1NM0N                                                            0.10      0.04     0.02     0.18 1.00     7484    12195
g1_LineROCKSTAR                                                                               0.11      0.04     0.03     0.19 1.00     7247    12034
g1_LineSUNCHASER                                                                              0.09      0.04     0.01     0.18 1.00     7062    11638
g1_Device_WateringD2_Limited                                                                 -0.14      0.07    -0.26     0.00 1.00     8662    10798
g1_Room_TemperatureCER8_Ambient                                                               0.01      0.06    -0.11     0.15 1.00    10108    13099
g1_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited                                  0.12      0.06     0.01     0.23 1.00     8697    13247
g1_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited                                          0.23      0.05     0.12     0.34 1.00     8407    12780
g1_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                  -0.03      0.06    -0.13     0.08 1.00     8494    13186
g1_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                   0.19      0.06     0.08     0.30 1.00     8925    13000
g1_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited                                  0.25      0.06     0.15     0.36 1.00     8555    11160
g1_LinePBI19N007M0NM3N:Device_WateringD2_Limited                                              0.16      0.05     0.06     0.27 1.00     8986    13450
g1_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited                                  0.08      0.06    -0.03     0.19 1.00     8070    12506
g1_LineROCKSTAR:Device_WateringD2_Limited                                                     0.16      0.06     0.05     0.27 1.00     8785    13301
g1_LineSUNCHASER:Device_WateringD2_Limited                                                    0.18      0.06     0.07     0.29 1.00     8505    13111
g1_Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                                     0.04      0.08    -0.11     0.20 1.00    12646    14805
k_Intercept                                                                                   2.03      0.32     1.41     2.69 1.00    18670    17156
k_LinePBI17N002M0CM0NM0NM010NM16N                                                             2.12      0.73     0.85     3.74 1.00    29382    14264
k_LinePBI17N006M0NM0NM89N                                                                     1.13      0.39     0.41     1.92 1.00    21077    12115
k_LinePBI17N007M0CM0NM0NM010NM8N                                                              0.50      0.34     0.03     1.28 1.00    18562    11291
k_LinePBI17N013M0CM0NM0NM010NM8N                                                              0.15      0.14     0.00     0.53 1.00    24341    12762
k_LinePBI17N016M0CM0NM0NM010NM27N                                                             2.15      0.55     1.13     3.28 1.00    27715    14969
k_LinePBI19N007M0NM3N                                                                         4.68      0.63     3.52     5.96 1.00    29291    19257
k_LinePBIC15034M0CM24NM010NM1NM0N                                                             5.22      0.96     3.46     7.20 1.00    28533    17573
k_LineROCKSTAR                                                                                3.11      0.57     2.04     4.29 1.00    29309    16698
k_LineSUNCHASER                                                                               3.98      0.86     2.44     5.79 1.00    28394    17917
lambda_Intercept                                                                              5.47      0.45     4.54     6.32 1.00    24884    15710

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sigma     0.08      0.00     0.08     0.09 1.00    31612    18569

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```

``` r
pp_check(fit_on_final)
```

```
Using 10 posterior draws for ppc type 'dens_overlay' by default.
```



![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-1.png)

``` r
rstan::check_hmc_diagnostics(fit_on_final$fit)
```

```

Divergences:
```

```
0 of 24000 iterations ended with a divergence.
```

```

Tree depth:
```

```
0 of 24000 iterations saturated the maximum tree depth of 12.
```

```

Energy:
```

```
E-BFMI indicated no pathological behavior.
```

``` r
plot(fit_on_final, ask = FALSE)
```



![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-2.png)

![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-3.png)

![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-4.png)

![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-5.png)

![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-6.png)

![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-7.png)

![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-8.png)

![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-9.png)

![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-10.png)

![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-11.png)

![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-12.png)

![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-13.png)

![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-14.png)

![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-15.png)

![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-16.png)

![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-17.png)

``` r
##make the above plot prettier
li600.on.plot <- li600.on.dat
li600.on.plot$fitted <- fitted(fit_on_final)[,1]
li600.on.plot$fitted.fixed <- fitted(fit_on_final, re_formula = NA)[,1]

ggplot(li600.on.plot, aes(x= time.adj, y = gsw, group = Replicate)) +
  geom_point() +
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red")
```



![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-18.png)

``` r
ggplot(li600.on.plot, aes(x= time.adj, y = gsw, group = Replicate)) +
  geom_point() +
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red") +
  geom_point(aes(y=fitted.fixed), colour = "blue") +
  geom_line(aes(y=fitted.fixed), colour = "blue")
```



![plot of chunk unnamed-chunk-12](figures/increase/unnamed-chunk-12-19.png)

### Save data


``` r
save.image(saveFile)
```

