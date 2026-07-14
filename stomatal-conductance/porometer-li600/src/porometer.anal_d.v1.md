---
output:
  pdf_document: default
  html_document: default
---
# Wheat stomatal kinetics (heat & drought)
## --- Ed Chaplin (USyd), *Apr--June 2025* ---
### Porometer DATA - ANALYSIS OF LIGHT OFF PORTION

* Analysis of gs decrease, when light switched off for the second time

### Load the libraries


``` r
library(knitr)
#knitr::spin("porometer.anal_d.v1.r")
library(ggplot2)
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, 
               "fig.width" = 7.5, "fig.height" = 10, "comment"=NA, fig.path = "figures/increase/")
library(dplyr)
library(asreml)
library(brms) ##package used to fit bayesian regression model using Stan
```

### Set up constants and files


``` r
options(width=90)
load("../data/porometer.v1.Rdata")
saveFile <- "../data/porometer.anal_d.v1.Rdata"
```

### Fit a model to the entire dataset, for the light off portion, following light on


``` r
##find the minimum time when the light was on
min(li600.dat$time.elapsed[li600.dat$Light == "On"]) ##38.46667
```

```
[1] 38.46667
```

``` r
max(li600.dat$time.elapsed[li600.dat$Light == "On"]) ##101.5667
```

```
[1] 101.5667
```

``` r
hist(li600.dat$time.elapsed[li600.dat$Light == "Off" &
                             li600.dat$time.elapsed > 38.46667 &
                             li600.dat$time.elapsed < 101.5667])
```



![plot of chunk unnamed-chunk-3](figures/increase/unnamed-chunk-3-1.png)

``` r
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
```

#### Stomatal conductance model - McAusland from McAusland et al 2016
Now let's try the Baesian approach that was fitted to the li-6800 data


``` r
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
```



![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-1.png)

``` r
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
```

```
Compiling Stan program...
```

```
Start sampling
```

``` r
print(summary(fit_off_max))
```

```
 Family: gaussian 
  Links: mu = identity 
Formula: gsw ~ g0 + (g1 - g0) * exp(-exp(((lambda - time.adj)/k) + 1)) 
         g0 ~ 1 + Line * Device_Watering * Room_Temperature + (1 | CEBRoom:CEBSide + CEBRoom:CEBSide:CEBLocation + CEBRoom:CEBSide:CEBLocation:CEBSubUnit + CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell)
         g1 ~ 1 + Line * Device_Watering * Room_Temperature + (1 | CEBRoom:CEBSide + CEBRoom:CEBSide:CEBLocation + CEBRoom:CEBSide:CEBLocation:CEBSubUnit + CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell)
         k ~ 1 + Line * Device_Watering * Room_Temperature
         lambda ~ 1 + Line * Device_Watering * Room_Temperature
   Data: li600.off.dat (Number of observations: 1959) 
  Draws: 6 chains, each with iter = 3000; warmup = 2000; thin = 1;
         total post-warmup draws = 6000

Multilevel Hyperparameters:
~CEBRoom:CEBSide (Number of levels: 4) 
                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(g0_Intercept)     0.04      0.04     0.00     0.16 1.00     2735     3203
sd(g1_Intercept)     0.01      0.00     0.00     0.02 1.00     3711     3126

~CEBRoom:CEBSide:CEBLocation (Number of levels: 8) 
                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(g0_Intercept)     0.02      0.02     0.00     0.06 1.00     2759     3367
sd(g1_Intercept)     0.01      0.00     0.00     0.02 1.00     3054     3359

~CEBRoom:CEBSide:CEBLocation:CEBSubUnit (Number of levels: 16) 
                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(g0_Intercept)     0.03      0.01     0.00     0.06 1.00     1307     2023
sd(g1_Intercept)     0.02      0.00     0.01     0.02 1.00     4044     4538

~CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell (Number of levels: 159) 
                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(g0_Intercept)     0.07      0.01     0.06     0.09 1.00     2033     3147
sd(g1_Intercept)     0.00      0.00     0.00     0.01 1.00     3030     3458

Regression Coefficients:
                                                                                              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
g0_Intercept                                                                                      0.38      0.06     0.26     0.51 1.00     1442
g0_LinePBI17N002M0CM0NM0NM010NM16N                                                               -0.29      0.06    -0.40    -0.18 1.01     1524
g0_LinePBI17N006M0NM0NM89N                                                                       -0.08      0.06    -0.19     0.04 1.00     1439
g0_LinePBI17N007M0CM0NM0NM010NM8N                                                                -0.12      0.06    -0.22    -0.01 1.00     1486
g0_LinePBI17N013M0CM0NM0NM010NM8N                                                                -0.11      0.06    -0.22    -0.00 1.00     1404
g0_LinePBI17N016M0CM0NM0NM010NM27N                                                               -0.29      0.06    -0.40    -0.17 1.00     1575
g0_LinePBI19N007M0NM3N                                                                           -0.02      0.06    -0.14     0.09 1.00     1416
g0_LinePBIC15034M0CM24NM010NM1NM0N                                                               -0.17      0.06    -0.28    -0.05 1.00     1394
g0_LineROCKSTAR                                                                                  -0.09      0.06    -0.20     0.02 1.00     1396
g0_LineSUNCHASER                                                                                 -0.07      0.06    -0.18     0.05 1.00     1469
g0_Device_WateringD2_Limited                                                                     -0.36      0.06    -0.47    -0.24 1.01      985
g0_Room_TemperatureCER8_Ambient                                                                   0.18      0.08     0.02     0.34 1.00     1547
g0_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited                                      0.28      0.08     0.14     0.43 1.01     1371
g0_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited                                              0.09      0.08    -0.06     0.24 1.00     1305
g0_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                       0.10      0.08    -0.04     0.25 1.00     1349
g0_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                       0.17      0.08     0.02     0.33 1.01     1355
g0_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited                                      0.28      0.08     0.13     0.42 1.00     1461
g0_LinePBI19N007M0NM3N:Device_WateringD2_Limited                                                  0.10      0.08    -0.05     0.25 1.00     1251
g0_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited                                      0.20      0.08     0.05     0.36 1.01     1171
g0_LineROCKSTAR:Device_WateringD2_Limited                                                         0.19      0.08     0.04     0.34 1.00     1296
g0_LineSUNCHASER:Device_WateringD2_Limited                                                        0.10      0.08    -0.05     0.25 1.00     1360
g0_LinePBI17N002M0CM0NM0NM010NM16N:Room_TemperatureCER8_Ambient                                  -0.18      0.07    -0.33    -0.04 1.00     1859
g0_LinePBI17N006M0NM0NM89N:Room_TemperatureCER8_Ambient                                          -0.26      0.08    -0.41    -0.11 1.00     1513
g0_LinePBI17N007M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                   -0.30      0.07    -0.45    -0.16 1.00     1592
g0_LinePBI17N013M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                   -0.22      0.08    -0.37    -0.08 1.00     1731
g0_LinePBI17N016M0CM0NM0NM010NM27N:Room_TemperatureCER8_Ambient                                  -0.18      0.08    -0.33    -0.03 1.00     1672
g0_LinePBI19N007M0NM3N:Room_TemperatureCER8_Ambient                                              -0.22      0.07    -0.37    -0.07 1.00     1524
g0_LinePBIC15034M0CM24NM010NM1NM0N:Room_TemperatureCER8_Ambient                                  -0.16      0.08    -0.30    -0.01 1.00     1546
g0_LineROCKSTAR:Room_TemperatureCER8_Ambient                                                     -0.25      0.08    -0.40    -0.10 1.00     1565
g0_LineSUNCHASER:Room_TemperatureCER8_Ambient                                                    -0.22      0.08    -0.37    -0.07 1.00     1566
g0_Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                                        -0.09      0.08    -0.25     0.06 1.00     1088
g0_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient         0.10      0.10    -0.10     0.29 1.00     1423
g0_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                 0.18      0.10    -0.02     0.39 1.00     1374
g0_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          0.22      0.10     0.02     0.42 1.00     1385
g0_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          0.12      0.11    -0.09     0.33 1.01     1471
g0_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient         0.15      0.10    -0.04     0.35 1.00     1452
g0_LinePBI19N007M0NM3N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                     0.10      0.10    -0.11     0.29 1.00     1360
g0_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient         0.06      0.10    -0.15     0.26 1.01     1304
g0_LineROCKSTAR:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                            0.17      0.10    -0.03     0.37 1.00     1385
g0_LineSUNCHASER:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                           0.15      0.10    -0.04     0.35 1.00     1526
g1_Intercept                                                                                      0.00      0.01    -0.02     0.02 1.00     4506
g1_LinePBI17N002M0CM0NM0NM010NM16N                                                               -0.01      0.01    -0.02     0.00 1.00     6290
g1_LinePBI17N006M0NM0NM89N                                                                       -0.02      0.01    -0.03    -0.00 1.00     5883
g1_LinePBI17N007M0CM0NM0NM010NM8N                                                                -0.01      0.01    -0.02     0.01 1.00     5388
g1_LinePBI17N013M0CM0NM0NM010NM8N                                                                -0.02      0.01    -0.03    -0.01 1.00     5773
g1_LinePBI17N016M0CM0NM0NM010NM27N                                                               -0.01      0.01    -0.02     0.01 1.00     5566
g1_LinePBI19N007M0NM3N                                                                           -0.01      0.01    -0.03    -0.00 1.00     5675
g1_LinePBIC15034M0CM24NM010NM1NM0N                                                               -0.02      0.01    -0.03    -0.01 1.00     5713
g1_LineROCKSTAR                                                                                  -0.01      0.01    -0.03     0.00 1.00     5545
g1_LineSUNCHASER                                                                                 -0.01      0.01    -0.03    -0.00 1.00     5388
g1_Device_WateringD2_Limited                                                                     -0.00      0.01    -0.02     0.02 1.00     4157
g1_Room_TemperatureCER8_Ambient                                                                  -0.01      0.01    -0.03     0.01 1.00     4877
g1_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited                                     -0.00      0.01    -0.02     0.02 1.00     8591
g1_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited                                              0.00      0.01    -0.02     0.02 1.00     7780
g1_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                       0.00      0.01    -0.02     0.02 1.00     7856
g1_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                      -0.00      0.01    -0.02     0.02 1.00     8610
g1_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited                                     -0.00      0.01    -0.02     0.02 1.00     7528
g1_LinePBI19N007M0NM3N:Device_WateringD2_Limited                                                 -0.01      0.01    -0.03     0.01 1.00     8802
g1_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited                                      0.00      0.01    -0.02     0.02 1.00     8544
g1_LineROCKSTAR:Device_WateringD2_Limited                                                        -0.02      0.01    -0.04     0.00 1.00     8117
g1_LineSUNCHASER:Device_WateringD2_Limited                                                       -0.00      0.01    -0.02     0.02 1.00     7357
g1_LinePBI17N002M0CM0NM0NM010NM16N:Room_TemperatureCER8_Ambient                                   0.01      0.01    -0.01     0.02 1.00     7453
g1_LinePBI17N006M0NM0NM89N:Room_TemperatureCER8_Ambient                                          -0.02      0.01    -0.03     0.00 1.00     7429
g1_LinePBI17N007M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                   -0.00      0.01    -0.02     0.02 1.00     7865
g1_LinePBI17N013M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                   -0.01      0.01    -0.03     0.01 1.00     7898
g1_LinePBI17N016M0CM0NM0NM010NM27N:Room_TemperatureCER8_Ambient                                  -0.00      0.01    -0.02     0.02 1.00     7388
g1_LinePBI19N007M0NM3N:Room_TemperatureCER8_Ambient                                              -0.02      0.01    -0.04    -0.00 1.00     7120
g1_LinePBIC15034M0CM24NM010NM1NM0N:Room_TemperatureCER8_Ambient                                  -0.01      0.01    -0.03     0.01 1.00     6186
g1_LineROCKSTAR:Room_TemperatureCER8_Ambient                                                     -0.02      0.01    -0.04    -0.00 1.00     6960
g1_LineSUNCHASER:Room_TemperatureCER8_Ambient                                                    -0.02      0.01    -0.04    -0.00 1.00     6988
g1_Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                                         0.00      0.01    -0.02     0.03 1.00     5072
g1_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient         0.00      0.02    -0.03     0.03 1.00    10291
g1_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                 0.01      0.02    -0.02     0.04 1.00    10584
g1_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          0.01      0.02    -0.02     0.04 1.00    12093
g1_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          0.01      0.02    -0.02     0.04 1.00    11460
g1_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient        -0.01      0.02    -0.04     0.02 1.00    12038
g1_LinePBI19N007M0NM3N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                     0.01      0.02    -0.02     0.04 1.00    11453
g1_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient         0.00      0.02    -0.03     0.03 1.00    10672
g1_LineROCKSTAR:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                           -0.02      0.02    -0.05     0.01 1.00    11850
g1_LineSUNCHASER:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                          -0.01      0.02    -0.04     0.02 1.00    10844
k_Intercept                                                                                       1.98      0.10     1.79     2.17 1.00     8939
k_LinePBI17N002M0CM0NM0NM010NM16N                                                                 2.97      0.10     2.77     3.16 1.00    14386
k_LinePBI17N006M0NM0NM89N                                                                         2.80      0.10     2.60     3.00 1.00    13075
k_LinePBI17N007M0CM0NM0NM010NM8N                                                                  2.93      0.10     2.74     3.12 1.00    11022
k_LinePBI17N013M0CM0NM0NM010NM8N                                                                  2.88      0.10     2.68     3.08 1.00    12637
k_LinePBI17N016M0CM0NM0NM010NM27N                                                                 2.97      0.10     2.77     3.16 1.00    12413
k_LinePBI19N007M0NM3N                                                                             2.74      0.10     2.53     2.93 1.00    12751
k_LinePBIC15034M0CM24NM010NM1NM0N                                                                 2.86      0.10     2.66     3.06 1.00    12041
k_LineROCKSTAR                                                                                    2.83      0.10     2.63     3.02 1.00    13320
k_LineSUNCHASER                                                                                   2.84      0.10     2.64     3.03 1.00    13957
k_Device_WateringD2_Limited                                                                       2.89      0.10     2.69     3.09 1.00    12998
k_Room_TemperatureCER8_Ambient                                                                    2.54      0.10     2.35     2.73 1.00    10151
k_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited                                       3.00      0.10     2.80     3.20 1.00    14185
k_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited                                               3.00      0.10     2.80     3.19 1.00    13189
k_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                        3.00      0.10     2.81     3.19 1.00    15118
k_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                        2.99      0.10     2.80     3.18 1.00    11627
k_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited                                       2.99      0.10     2.80     3.19 1.00    12210
k_LinePBI19N007M0NM3N:Device_WateringD2_Limited                                                   2.98      0.10     2.78     3.18 1.00    12214
k_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited                                       2.99      0.10     2.79     3.19 1.00    11766
k_LineROCKSTAR:Device_WateringD2_Limited                                                          2.95      0.10     2.76     3.15 1.00    12476
k_LineSUNCHASER:Device_WateringD2_Limited                                                         2.99      0.10     2.78     3.19 1.00    13662
k_LinePBI17N002M0CM0NM0NM010NM16N:Room_TemperatureCER8_Ambient                                    2.98      0.10     2.79     3.18 1.00    11938
k_LinePBI17N006M0NM0NM89N:Room_TemperatureCER8_Ambient                                            2.93      0.10     2.74     3.13 1.00    12864
k_LinePBI17N007M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                     2.98      0.10     2.78     3.17 1.00    13380
k_LinePBI17N013M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                     2.93      0.10     2.73     3.13 1.00    13623
k_LinePBI17N016M0CM0NM0NM010NM27N:Room_TemperatureCER8_Ambient                                    2.99      0.10     2.79     3.18 1.00    12727
k_LinePBI19N007M0NM3N:Room_TemperatureCER8_Ambient                                                2.90      0.10     2.71     3.09 1.00    13235
k_LinePBIC15034M0CM24NM010NM1NM0N:Room_TemperatureCER8_Ambient                                    2.92      0.10     2.72     3.11 1.00    14729
k_LineROCKSTAR:Room_TemperatureCER8_Ambient                                                       2.94      0.10     2.74     3.13 1.00    10882
k_LineSUNCHASER:Room_TemperatureCER8_Ambient                                                      2.92      0.10     2.73     3.12 1.00    12612
k_Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                                          2.96      0.10     2.77     3.15 1.00    11988
k_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          3.00      0.10     2.80     3.19 1.00    13974
k_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                  3.00      0.10     2.80     3.20 1.00    13495
k_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient           3.00      0.10     2.81     3.19 1.00    13290
k_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient           3.00      0.10     2.80     3.19 1.00    13603
k_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          2.99      0.10     2.80     3.19 1.00    13708
k_LinePBI19N007M0NM3N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                      3.00      0.10     2.79     3.21 1.00    13763
k_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          3.00      0.10     2.80     3.19 1.00    12734
k_LineROCKSTAR:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                             2.98      0.10     2.79     3.18 1.00    12494
k_LineSUNCHASER:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                            3.00      0.10     2.80     3.19 1.00    11073
lambda_Intercept                                                                                  0.02      0.02     0.00     0.09 1.00     8022
lambda_LinePBI17N002M0CM0NM0NM010NM16N                                                            0.42      0.22     0.04     0.89 1.00     3862
lambda_LinePBI17N006M0NM0NM89N                                                                    0.16      0.13     0.01     0.49 1.00     5643
lambda_LinePBI17N007M0CM0NM0NM010NM8N                                                             0.35      0.20     0.03     0.79 1.00     4333
lambda_LinePBI17N013M0CM0NM0NM010NM8N                                                             0.28      0.18     0.02     0.68 1.00     5243
lambda_LinePBI17N016M0CM0NM0NM010NM27N                                                            0.41      0.22     0.04     0.86 1.00     4501
lambda_LinePBI19N007M0NM3N                                                                        0.12      0.10     0.00     0.38 1.00     6208
lambda_LinePBIC15034M0CM24NM010NM1NM0N                                                            0.25      0.17     0.01     0.65 1.00     4997
lambda_LineROCKSTAR                                                                               0.18      0.14     0.01     0.51 1.00     5966
lambda_LineSUNCHASER                                                                              0.22      0.16     0.01     0.59 1.00     5193
lambda_Device_WateringD2_Limited                                                                  0.18      0.14     0.01     0.53 1.00     5049
lambda_Room_TemperatureCER8_Ambient                                                               0.06      0.05     0.00     0.20 1.00     6462
lambda_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited                                  0.50      0.23     0.08     0.98 1.00     3496
lambda_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited                                          0.49      0.23     0.06     0.95 1.00     4945
lambda_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                   0.50      0.23     0.07     0.98 1.00     4299
lambda_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                   0.48      0.23     0.06     0.94 1.00     4530
lambda_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited                                  0.48      0.24     0.06     0.97 1.00     4278
lambda_LinePBI19N007M0NM3N:Device_WateringD2_Limited                                              0.45      0.23     0.05     0.90 1.00     4008
lambda_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited                                  0.48      0.22     0.08     0.93 1.00     4944
lambda_LineROCKSTAR:Device_WateringD2_Limited                                                     0.37      0.21     0.03     0.82 1.00     4676
lambda_LineSUNCHASER:Device_WateringD2_Limited                                                    0.45      0.23     0.05     0.93 1.00     4436
lambda_LinePBI17N002M0CM0NM0NM010NM16N:Room_TemperatureCER8_Ambient                               0.45      0.23     0.06     0.93 1.00     4310
                                                                                              Tail_ESS
g0_Intercept                                                                                      2356
g0_LinePBI17N002M0CM0NM0NM010NM16N                                                                3040
g0_LinePBI17N006M0NM0NM89N                                                                        2855
g0_LinePBI17N007M0CM0NM0NM010NM8N                                                                 3149
g0_LinePBI17N013M0CM0NM0NM010NM8N                                                                 3030
g0_LinePBI17N016M0CM0NM0NM010NM27N                                                                2468
g0_LinePBI19N007M0NM3N                                                                            2847
g0_LinePBIC15034M0CM24NM010NM1NM0N                                                                2365
g0_LineROCKSTAR                                                                                   2662
g0_LineSUNCHASER                                                                                  2548
g0_Device_WateringD2_Limited                                                                      2555
g0_Room_TemperatureCER8_Ambient                                                                   2514
g0_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited                                      3093
g0_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited                                              2531
g0_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                       2911
g0_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                       2873
g0_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited                                      2922
g0_LinePBI19N007M0NM3N:Device_WateringD2_Limited                                                  2515
g0_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited                                      2888
g0_LineROCKSTAR:Device_WateringD2_Limited                                                         2935
g0_LineSUNCHASER:Device_WateringD2_Limited                                                        3013
g0_LinePBI17N002M0CM0NM0NM010NM16N:Room_TemperatureCER8_Ambient                                   2674
g0_LinePBI17N006M0NM0NM89N:Room_TemperatureCER8_Ambient                                           2600
g0_LinePBI17N007M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                    2715
g0_LinePBI17N013M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                    3208
g0_LinePBI17N016M0CM0NM0NM010NM27N:Room_TemperatureCER8_Ambient                                   2526
g0_LinePBI19N007M0NM3N:Room_TemperatureCER8_Ambient                                               2909
g0_LinePBIC15034M0CM24NM010NM1NM0N:Room_TemperatureCER8_Ambient                                   2560
g0_LineROCKSTAR:Room_TemperatureCER8_Ambient                                                      2710
g0_LineSUNCHASER:Room_TemperatureCER8_Ambient                                                     3088
g0_Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                                         2561
g0_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient         2251
g0_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                 2349
g0_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          3122
g0_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          3099
g0_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient         2957
g0_LinePBI19N007M0NM3N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                     2563
g0_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient         2719
g0_LineROCKSTAR:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                            2845
g0_LineSUNCHASER:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                           2581
g1_Intercept                                                                                      3948
g1_LinePBI17N002M0CM0NM0NM010NM16N                                                                4738
g1_LinePBI17N006M0NM0NM89N                                                                        4406
g1_LinePBI17N007M0CM0NM0NM010NM8N                                                                 4757
g1_LinePBI17N013M0CM0NM0NM010NM8N                                                                 4973
g1_LinePBI17N016M0CM0NM0NM010NM27N                                                                4935
g1_LinePBI19N007M0NM3N                                                                            4887
g1_LinePBIC15034M0CM24NM010NM1NM0N                                                                4910
g1_LineROCKSTAR                                                                                   4696
g1_LineSUNCHASER                                                                                  5103
g1_Device_WateringD2_Limited                                                                      3891
g1_Room_TemperatureCER8_Ambient                                                                   4250
g1_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited                                      5274
g1_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited                                              4884
g1_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                       4597
g1_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                       4977
g1_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited                                      4493
g1_LinePBI19N007M0NM3N:Device_WateringD2_Limited                                                  5347
g1_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited                                      5240
g1_LineROCKSTAR:Device_WateringD2_Limited                                                         5123
g1_LineSUNCHASER:Device_WateringD2_Limited                                                        4835
g1_LinePBI17N002M0CM0NM0NM010NM16N:Room_TemperatureCER8_Ambient                                   4932
g1_LinePBI17N006M0NM0NM89N:Room_TemperatureCER8_Ambient                                           4946
g1_LinePBI17N007M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                    5082
g1_LinePBI17N013M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                    4712
g1_LinePBI17N016M0CM0NM0NM010NM27N:Room_TemperatureCER8_Ambient                                   4497
g1_LinePBI19N007M0NM3N:Room_TemperatureCER8_Ambient                                               4796
g1_LinePBIC15034M0CM24NM010NM1NM0N:Room_TemperatureCER8_Ambient                                   4884
g1_LineROCKSTAR:Room_TemperatureCER8_Ambient                                                      4630
g1_LineSUNCHASER:Room_TemperatureCER8_Ambient                                                     5302
g1_Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                                         4436
g1_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient         4051
g1_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                 4719
g1_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          4621
g1_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          4922
g1_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient         4530
g1_LinePBI19N007M0NM3N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                     4770
g1_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient         4894
g1_LineROCKSTAR:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                            4600
g1_LineSUNCHASER:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                           4417
k_Intercept                                                                                       4896
k_LinePBI17N002M0CM0NM0NM010NM16N                                                                 4166
k_LinePBI17N006M0NM0NM89N                                                                         4797
k_LinePBI17N007M0CM0NM0NM010NM8N                                                                  4186
k_LinePBI17N013M0CM0NM0NM010NM8N                                                                  4390
k_LinePBI17N016M0CM0NM0NM010NM27N                                                                 4294
k_LinePBI19N007M0NM3N                                                                             4105
k_LinePBIC15034M0CM24NM010NM1NM0N                                                                 4327
k_LineROCKSTAR                                                                                    4484
k_LineSUNCHASER                                                                                   4144
k_Device_WateringD2_Limited                                                                       4438
k_Room_TemperatureCER8_Ambient                                                                    4381
k_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited                                       4506
k_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited                                               4314
k_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                        4305
k_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                        4195
k_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited                                       3937
k_LinePBI19N007M0NM3N:Device_WateringD2_Limited                                                   3884
k_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited                                       4034
k_LineROCKSTAR:Device_WateringD2_Limited                                                          4497
k_LineSUNCHASER:Device_WateringD2_Limited                                                         4267
k_LinePBI17N002M0CM0NM0NM010NM16N:Room_TemperatureCER8_Ambient                                    3724
k_LinePBI17N006M0NM0NM89N:Room_TemperatureCER8_Ambient                                            4290
k_LinePBI17N007M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                     4181
k_LinePBI17N013M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                                     4048
k_LinePBI17N016M0CM0NM0NM010NM27N:Room_TemperatureCER8_Ambient                                    3865
k_LinePBI19N007M0NM3N:Room_TemperatureCER8_Ambient                                                4388
k_LinePBIC15034M0CM24NM010NM1NM0N:Room_TemperatureCER8_Ambient                                    4555
k_LineROCKSTAR:Room_TemperatureCER8_Ambient                                                       4271
k_LineSUNCHASER:Room_TemperatureCER8_Ambient                                                      4075
k_Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                                          3920
k_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          4229
k_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                  4368
k_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient           4229
k_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient           4147
k_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          4429
k_LinePBI19N007M0NM3N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                      4285
k_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient          4360
k_LineROCKSTAR:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                             4032
k_LineSUNCHASER:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                            4375
lambda_Intercept                                                                                  3349
lambda_LinePBI17N002M0CM0NM0NM010NM16N                                                            1703
lambda_LinePBI17N006M0NM0NM89N                                                                    3287
lambda_LinePBI17N007M0CM0NM0NM010NM8N                                                             2588
lambda_LinePBI17N013M0CM0NM0NM010NM8N                                                             3302
lambda_LinePBI17N016M0CM0NM0NM010NM27N                                                            2363
lambda_LinePBI19N007M0NM3N                                                                        3287
lambda_LinePBIC15034M0CM24NM010NM1NM0N                                                            2699
lambda_LineROCKSTAR                                                                               3238
lambda_LineSUNCHASER                                                                              2972
lambda_Device_WateringD2_Limited                                                                  2828
lambda_Room_TemperatureCER8_Ambient                                                               3473
lambda_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited                                  1798
lambda_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited                                          2277
lambda_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                   2334
lambda_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                   1827
lambda_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited                                  1859
lambda_LinePBI19N007M0NM3N:Device_WateringD2_Limited                                              1801
lambda_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited                                  2360
lambda_LineROCKSTAR:Device_WateringD2_Limited                                                     2716
lambda_LineSUNCHASER:Device_WateringD2_Limited                                                    2406
lambda_LinePBI17N002M0CM0NM0NM010NM16N:Room_TemperatureCER8_Ambient                               2105
 [ reached 'max' / getOption("max.print") -- omitted 18 rows ]

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sigma     0.04      0.00     0.04     0.04 1.00     6566     4943

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```

``` r
plot(fit_off_max, ask = FALSE)
```



![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-2.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-3.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-4.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-5.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-6.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-7.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-8.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-9.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-10.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-11.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-12.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-13.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-14.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-15.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-16.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-17.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-18.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-19.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-20.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-21.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-22.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-23.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-24.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-25.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-26.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-27.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-28.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-29.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-30.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-31.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-32.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-33.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-34.png)

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-35.png)

``` r
#pairs(fit_off_max)
#mcmc_plot(fit_off_max, type = "trace")
rstan::check_hmc_diagnostics(fit_off_max$fit)
```

```

Divergences:
```

```
0 of 6000 iterations ended with a divergence.
```

```

Tree depth:
```

```
0 of 6000 iterations saturated the maximum tree depth of 12.
```

```

Energy:
```

```
E-BFMI indicated no pathological behavior.
```

``` r
loo(fit_off_max)
```

```
Warning: Found 61 observations with a pareto_k > 0.7 in model 'fit_off_max'. We recommend to set 'moment_match = TRUE' in order to perform moment
matching for problematic observations.
```

```

Computed from 6000 by 1959 log-likelihood matrix.

         Estimate    SE
elpd_loo   3199.2  85.0
p_loo       395.2  35.9
looic     -6398.3 169.9
------
MCSE of elpd_loo is NA.
MCSE and ESS estimates assume MCMC draws (r_eff in [0.5, 1.6]).

Pareto k diagnostic values:
                         Count Pct.    Min. ESS
(-Inf, 0.7]   (good)     1898  96.9%   133     
   (0.7, 1]   (bad)        41   2.1%   <NA>    
   (1, Inf)   (very bad)   20   1.0%   <NA>    
See help('pareto-k-diagnostic') for details.
```

``` r
pp_check(fit_off_max)  # compares density of the observed data with data generated from posterior predictive distribution
```

```
Using 10 posterior draws for ppc type 'dens_overlay' by default.
```



![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-36.png)

``` r
pp_check(fit_off_max, type = "error_hist", ndraws = 12)
```

```
`stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```



![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-37.png)

``` r
pp_check(fit_off_max, type = "scatter_avg", ndraws = 1000)
```



![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-38.png)

``` r
pp_check(fit_off_max, type = "stat_2d")
```

```
Using all posterior draws for ppc type 'stat_2d' by default.
```

```
Note: in most cases the default test statistic 'mean' is too weak to detect anything of interest.
```



![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-39.png)

``` r
pp_check(fit_off_max, type = "loo_ribbon")
```

```
Using all posterior draws for ppc type 'loo_ribbon' by default.
```

```
Warning: Found 61 observations with a pareto_k > 0.7 in model '.x1'. We recommend to set 'moment_match = TRUE' in order to perform moment matching for
problematic observations.
```



![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-40.png)

``` r
plot(gsw~time.adj, li600.off.dat, pch=16)
points(fitted(fit_off_max)[,1]~time.adj, li600.off.dat, pch=16, col="red")
```



![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-41.png)

``` r
##make the above plot prettier
li600.off.plot <- li600.off.dat
li600.off.plot$fitted <- fitted(fit_off_max)[,1]
li600.off.plot$fitted.fixed <- fitted(fit_off_max, re_formula = NA)[,1]

ggplot(li600.off.plot, aes(x= time.adj, y = gsw, group = Replicate)) +
  geom_point() +
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red")
```



![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-42.png)

``` r
ggplot(li600.off.plot, aes(x= time.adj, y = gsw, group = Replicate)) +
  geom_point() +
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red") +
  geom_point(aes(y=fitted.fixed), colour = "blue") +
  geom_line(aes(y=fitted.fixed), colour = "blue")
```



![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-43.png)

``` r
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
```

#### Fit final model


``` r
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
```

```
Compiling Stan program...
```

```
Start sampling
```

``` r
print(summary(fit_off_final))
```

```
 Family: gaussian 
  Links: mu = identity 
Formula: gsw ~ g0 + (g1 - g0) * exp(-exp(((lambda - time.adj)/k) + 1)) 
         g0 ~ 1 + Line * Device_Watering * Room_Temperature + (1 | CEBRoom:CEBSide:CEBLocation + CEBRoom:CEBSide:CEBLocation:CEBSubUnit + CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell)
         g1 ~ 1 + Line + Device_Watering + Line:Device_Watering + (1 | CEBRoom:CEBSide:CEBLocation + CEBRoom:CEBSide:CEBLocation:CEBSubUnit + CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell)
         k ~ 1 + Room_Temperature
         lambda ~ 1 + Device_Watering
   Data: li600.off.dat (Number of observations: 1959) 
  Draws: 6 chains, each with iter = 10000; warmup = 6000; thin = 1;
         total post-warmup draws = 24000

Multilevel Hyperparameters:
~CEBRoom:CEBSide:CEBLocation (Number of levels: 8) 
                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(g0_Intercept)     0.04      0.03     0.00     0.11 1.00    11691    14175
sd(g1_Intercept)     0.00      0.00     0.00     0.01 1.00     9771    14474

~CEBRoom:CEBSide:CEBLocation:CEBSubUnit (Number of levels: 16) 
                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(g0_Intercept)     0.06      0.02     0.01     0.11 1.00     6481     6006
sd(g1_Intercept)     0.01      0.00     0.01     0.02 1.00    13702    17666

~CEBRoom:CEBSide:CEBLocation:CEBSubUnit:CEBCell (Number of levels: 159) 
                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(g0_Intercept)     0.14      0.01     0.12     0.16 1.00     8078    14856
sd(g1_Intercept)     0.01      0.00     0.01     0.01 1.00    11004    15780

Regression Coefficients:
                                                                                          Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
g0_Intercept                                                                                  0.36      0.07     0.22     0.49 1.00     6443    12020
g0_LinePBI17N002M0CM0NM0NM010NM16N                                                           -0.17      0.09    -0.34    -0.00 1.00     8002    12037
g0_LinePBI17N006M0NM0NM89N                                                                    0.06      0.09    -0.11     0.23 1.00     7861    12553
g0_LinePBI17N007M0CM0NM0NM010NM8N                                                             0.02      0.09    -0.14     0.19 1.00     7841    12448
g0_LinePBI17N013M0CM0NM0NM010NM8N                                                             0.02      0.09    -0.14     0.19 1.00     7824    13393
g0_LinePBI17N016M0CM0NM0NM010NM27N                                                           -0.20      0.08    -0.36    -0.03 1.00     7531    12527
g0_LinePBI19N007M0NM3N                                                                        0.13      0.08    -0.04     0.30 1.00     8240    13140
g0_LinePBIC15034M0CM24NM010NM1NM0N                                                            0.04      0.09    -0.13     0.21 1.00     7999    12602
g0_LineROCKSTAR                                                                               0.05      0.08    -0.12     0.22 1.00     7334    12203
g0_LineSUNCHASER                                                                              0.12      0.09    -0.05     0.29 1.00     7899    12970
g0_Device_WateringD2_Limited                                                                 -0.31      0.09    -0.48    -0.14 1.00     7130    12112
g0_Room_TemperatureCER8_Ambient                                                               0.17      0.09    -0.00     0.35 1.00     7082    12651
g0_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited                                  0.22      0.12    -0.01     0.46 1.00     8644    12115
g0_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited                                          0.03      0.12    -0.20     0.26 1.00     9435    13619
g0_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                  -0.00      0.12    -0.23     0.23 1.00     8909    13634
g0_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                   0.08      0.12    -0.16     0.32 1.00     9218    13636
g0_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited                                  0.21      0.12    -0.02     0.44 1.00     9038    13611
g0_LinePBI19N007M0NM3N:Device_WateringD2_Limited                                              0.09      0.12    -0.15     0.32 1.00     9289    14349
g0_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited                                  0.15      0.12    -0.08     0.39 1.00     9501    13022
g0_LineROCKSTAR:Device_WateringD2_Limited                                                     0.15      0.12    -0.07     0.38 1.00     8403    12823
g0_LineSUNCHASER:Device_WateringD2_Limited                                                    0.03      0.12    -0.21     0.26 1.00     8829    12357
g0_LinePBI17N002M0CM0NM0NM010NM16N:Room_TemperatureCER8_Ambient                              -0.20      0.12    -0.43     0.04 1.00     8591    11609
g0_LinePBI17N006M0NM0NM89N:Room_TemperatureCER8_Ambient                                      -0.18      0.12    -0.41     0.04 1.00     8329    13248
g0_LinePBI17N007M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                               -0.29      0.12    -0.52    -0.06 1.00     7777    12299
g0_LinePBI17N013M0CM0NM0NM010NM8N:Room_TemperatureCER8_Ambient                               -0.16      0.12    -0.38     0.07 1.00     8333    13272
g0_LinePBI17N016M0CM0NM0NM010NM27N:Room_TemperatureCER8_Ambient                              -0.14      0.12    -0.37     0.09 1.00     8193    12730
g0_LinePBI19N007M0NM3N:Room_TemperatureCER8_Ambient                                          -0.12      0.12    -0.34     0.11 1.00     8646    12570
g0_LinePBIC15034M0CM24NM010NM1NM0N:Room_TemperatureCER8_Ambient                              -0.15      0.12    -0.38     0.08 1.00     8195    12352
g0_LineROCKSTAR:Room_TemperatureCER8_Ambient                                                 -0.19      0.12    -0.42     0.04 1.00     8037    12363
g0_LineSUNCHASER:Room_TemperatureCER8_Ambient                                                -0.17      0.12    -0.40     0.06 1.00     8300    12565
g0_Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                                    -0.09      0.11    -0.31     0.13 1.00     8594    11627
g0_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient     0.07      0.16    -0.24     0.38 1.00     9917    13355
g0_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient             0.06      0.16    -0.25     0.37 1.00    10280    13586
g0_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient      0.16      0.16    -0.16     0.46 1.00    10313    13794
g0_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient      0.12      0.16    -0.20     0.43 1.00    10448    13313
g0_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient     0.16      0.16    -0.14     0.47 1.00    10402    13789
g0_LinePBI19N007M0NM3N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                -0.10      0.16    -0.42     0.20 1.00    10388    13884
g0_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient    -0.05      0.16    -0.36     0.26 1.00    10418    13710
g0_LineROCKSTAR:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                        0.19      0.16    -0.12     0.50 1.00    10031    14495
g0_LineSUNCHASER:Device_WateringD2_Limited:Room_TemperatureCER8_Ambient                       0.09      0.16    -0.23     0.40 1.00    10351    14151
g1_Intercept                                                                                  0.00      0.01    -0.01     0.01 1.00    11611    14403
g1_LinePBI17N002M0CM0NM0NM010NM16N                                                           -0.00      0.00    -0.01     0.01 1.00    13692    16352
g1_LinePBI17N006M0NM0NM89N                                                                   -0.00      0.00    -0.01     0.01 1.00    14314    17789
g1_LinePBI17N007M0CM0NM0NM010NM8N                                                             0.01      0.00    -0.00     0.01 1.00    12776    16519
g1_LinePBI17N013M0CM0NM0NM010NM8N                                                            -0.00      0.00    -0.01     0.00 1.00    13997    16852
g1_LinePBI17N016M0CM0NM0NM010NM27N                                                           -0.00      0.00    -0.01     0.01 1.00    13243    17539
g1_LinePBI19N007M0NM3N                                                                        0.00      0.00    -0.01     0.01 1.00    13504    17256
g1_LinePBIC15034M0CM24NM010NM1NM0N                                                           -0.00      0.00    -0.01     0.00 1.00    14138    16610
g1_LineROCKSTAR                                                                               0.00      0.00    -0.01     0.01 1.00    13856    16070
g1_LineSUNCHASER                                                                              0.00      0.00    -0.01     0.01 1.00    13155    17009
g1_Device_WateringD2_Limited                                                                  0.00      0.01    -0.01     0.02 1.00    11125    15322
g1_LinePBI17N002M0CM0NM0NM010NM16N:Device_WateringD2_Limited                                 -0.00      0.01    -0.02     0.01 1.00    15084    18120
g1_LinePBI17N006M0NM0NM89N:Device_WateringD2_Limited                                         -0.00      0.01    -0.01     0.01 1.00    16469    18051
g1_LinePBI17N007M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                  -0.01      0.01    -0.02     0.00 1.00    15448    18302
g1_LinePBI17N013M0CM0NM0NM010NM8N:Device_WateringD2_Limited                                   0.00      0.01    -0.01     0.01 1.00    15801    17077
g1_LinePBI17N016M0CM0NM0NM010NM27N:Device_WateringD2_Limited                                 -0.01      0.01    -0.02     0.01 1.00    14718    17325
g1_LinePBI19N007M0NM3N:Device_WateringD2_Limited                                             -0.01      0.01    -0.02     0.00 1.00    16226    16851
g1_LinePBIC15034M0CM24NM010NM1NM0N:Device_WateringD2_Limited                                 -0.00      0.01    -0.01     0.01 1.00    15625    17240
g1_LineROCKSTAR:Device_WateringD2_Limited                                                    -0.01      0.01    -0.02     0.00 1.00    15659    16868
g1_LineSUNCHASER:Device_WateringD2_Limited                                                   -0.01      0.01    -0.02     0.00 1.00    15856    17897
k_Intercept                                                                                   2.62      0.04     2.54     2.69 1.00    39789    20523
k_Room_TemperatureCER8_Ambient                                                                1.97      0.06     1.85     2.09 1.00    35775    20681
lambda_Intercept                                                                              0.02      0.02     0.00     0.09 1.00    26708    13086
lambda_Device_WateringD2_Limited                                                              0.03      0.03     0.00     0.11 1.00    27002    12874

Further Distributional Parameters:
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sigma     0.02      0.00     0.02     0.02 1.00    28903    19031

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```

``` r
pp_check(fit_off_final)
```

```
Using 10 posterior draws for ppc type 'dens_overlay' by default.
```



![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-1.png)

``` r
rstan::check_hmc_diagnostics(fit_off_final$fit)
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
0 of 24000 iterations saturated the maximum tree depth of 15.
```

```

Energy:
```

```
E-BFMI indicated no pathological behavior.
```

``` r
plot(fit_off_final, ask = FALSE)
```



![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-2.png)

![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-3.png)

![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-4.png)

![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-5.png)

![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-6.png)

![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-7.png)

![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-8.png)

![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-9.png)

![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-10.png)

![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-11.png)

![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-12.png)

![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-13.png)

![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-14.png)

![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-15.png)

![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-16.png)

``` r
##make the above plot prettier
li600.off.plot <- li600.off.dat
li600.off.plot$fitted <- fitted(fit_off_final)[,1]
li600.off.plot$fitted.fixed <- fitted(fit_off_final, re_formula = NA)[,1]

ggplot(li600.off.plot, aes(x= time.adj, y = gsw, group = Replicate)) +
  geom_point() +
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red")
```



![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-17.png)

``` r
ggplot(li600.off.plot, aes(x= time.adj, y = gsw, group = Replicate)) +
  geom_point() +
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  geom_point(aes(y=fitted), colour = "red") +
  geom_line(aes(y=fitted), colour = "red") +
  geom_point(aes(y=fitted.fixed), colour = "blue") +
  geom_line(aes(y=fitted.fixed), colour = "blue")
```



![plot of chunk unnamed-chunk-5](figures/increase/unnamed-chunk-5-18.png)

### Save data


``` r
save.image(saveFile)
```

