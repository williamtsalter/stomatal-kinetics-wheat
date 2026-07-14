# Wheat stomatal kinetics (heat & drought)
## --- Ed Chaplin (USyd), *Apr--June 2025* ---
### Li-600 STOMATAL CONDUCTANCE MODEL PARAMETER ESTIMATION - LIGHT ON TO LIGHT OFF (DECREASE)

* Using the final model fitted to the data, obtain the parameters of interest
* These include: 
* k_d, being the time for gs to decrease to Gsmin
* lambda, being the lag in the time of gs to a step decrease in light
* g0, being the starting value of the gs (intercept/upper asymptote)
* g1, being the lower asymptote of the gs
* Slmax, a derived trait from McAusland paper
* 

### Load the libraries


``` r
library(knitr)
#knitr::spin("porometer.stomatalCondParamEstimation_d.v1.r")
library(ggplot2)
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, 
               "fig.width" = 7.5, "fig.height" = 10, "comment"=NA, fig.path = "figures/increase/")
library(dae)
library(dplyr)
library(brms) ##package used to fit bayesian regress.off model using Stan
library(tidyr)
library(emmeans)
library(tidybayes)
```

### Set up .offstants and files


``` r
opt.offs(width=90)
```

```
Error in `opt.offs()`:
! could not find function "opt.offs"
```

``` r
load("../data/porometer.anal_d.v1.Rdata")
saveFile <- "../data/stomatalCondParamEstimation_porometer_d.v1.Rdata"
```

### Have a look at the final model summary and store fixed effects


``` r
summary(fit_off_final)
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

### Build a dataframe of model parameters


``` r
li600.d.param <- expand.grid(Room_Temperature = levels(li600.off.dat$Room_Temperature),
                             Device_Watering = levels(li600.off.dat$Device_Watering),
                             Line = levels(li600.off.dat$Line))
```

### Use the emmeans package to get out the expected parameter estimates, a.offg with HPD (credible intervals)


``` r
g0.emmeans <- as.data.frame(emmeans(fit_off_final,
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
li600.d.param <- left_join(li600.d.param, g0.emmeans)
```

```
Joining with `by = join_by(Room_Temperature, Device_Watering, Line)`
```

``` r
g1.emmeans <- as.data.frame(emmeans(fit_off_final,
                                    ~ Device_Watering * Line,
                                    nlpar = "g1"))
g1.emmeans <- g1.emmeans %>% rename_with(~paste0("g1_", .x), 
                                         -all_of(c("Device_Watering",
                                                   "Line")))
g1.emmeans <- within(g1.emmeans, {
  g1_tabquant <- paste0("'", round(g1_emmean, 3), "\n[", 
                        round(g1_lower.HPD, 3), ", ", 
                        round(g1_upper.HPD, 3), "]'")
  
})
li600.d.param <- left_join(li600.d.param, g1.emmeans)
```

```
Joining with `by = join_by(Device_Watering, Line)`
```

``` r
k.emmeans <- as.data.frame(emmeans(fit_off_final,
                                    ~ Room_Temperature,
                                    nlpar = "k"))
k.emmeans <- k.emmeans %>% rename_with(~paste0("k_", .x), 
                                       -all_of(c("Room_Temperature")))
k.emmeans <- within(k.emmeans, {
  k_tabquant <- paste0("'", round(k_emmean, 3), "\n[", 
                        round(k_lower.HPD, 3), ", ", 
                        round(k_upper.HPD, 3), "]'")
  
})
li600.d.param <- left_join(li600.d.param, k.emmeans)
```

```
Joining with `by = join_by(Room_Temperature)`
```

``` r
lambda.emmeans <- as.data.frame(emmeans(fit_off_final,
                                   ~ Device_Watering,
                                   nlpar = "lambda"))
lambda.emmeans <- lambda.emmeans %>% rename_with(~paste0("lambda_", .x), 
                                       -all_of(c("Device_Watering")))
lambda.emmeans <- within(lambda.emmeans, {
  lambda_tabquant <- paste0("'", round(lambda.emmean, 3), "\n[", 
                       round(lambda.lower.HPD, 3), ", ", 
                       round(lambda.upper.HPD, 3), "]'")
  
})
```

```
Error:
! object 'lambda.emmean' not found
```

``` r
li600.d.param <- left_join(li600.d.param, lambda.emmeans)
```

```
Joining with `by = join_by(Device_Watering)`
```

### Using these parameters, calculate Slmax
* Slmax looks to be reported in the units mmol m^-2 s^-2
* This involves multiplying g1 and g0 by 1000 and also k by 60
* I can calculate this using the posterior draws for each parameter and then combine them together


``` r
g0_draws <- emmeans(fit_off_final,
                    ~ Room_Temperature * Device_Watering * Line,
                    nlpar = "g0") %>%
  gather_emmeans_draws()

g1_draws <- emmeans(fit_off_final,
                    ~ Device_Watering * Line,
                    nlpar = "g1") %>%
  gather_emmeans_draws()

k_draws <- emmeans(fit_off_final,
                   ~ Room_Temperature,
                   nlpar = "k") %>%
  gather_emmeans_draws()

combined_draws <- g0_draws %>%
  rename(g0 = .value) %>%
  left_join(g1_draws %>% rename(g1 = .value),
            by = c(".draw", "Device_Watering", "Line")) %>%
  left_join(k_draws %>% rename(k = .value),
            by = c(".draw", "Room_Temperature"))

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

li600.d.param <- left_join(li600.d.param,
                            summary_slopes %>% select(Room_Temperature,
                                                      Device_Watering,
                                                      Line,
                                                      Slmax_emmean,
                                                      Slmax_lower.HPD,
                                                      Slmax_upper.HPD,
                                                      Slmax_tabquant))
```

```
Joining with `by = join_by(Room_Temperature, Device_Watering, Line)`
```

``` r
li600.d.param
```

```
   Room_Temperature Device_Watering                        Line  g0_emmean g0_lower.HPD g0_upper.HPD              g0_tabquant     g1_emmean g1_lower.HPD g1_upper.HPD
1         CER7_Heat    D1_Unlimited PBI16C009-0C-0N-0N-010N-14N 0.35504387  0.223190153    0.4875488  '0.355\n[0.223, 0.488]'  0.0006697599 -0.009754446  0.010568313
2      CER8_Ambient    D1_Unlimited PBI16C009-0C-0N-0N-010N-14N 0.53057635  0.387970458    0.6692195  '0.531\n[0.388, 0.669]'  0.0006697599 -0.009754446  0.010568313
3         CER7_Heat      D2_Limited PBI16C009-0C-0N-0N-010N-14N 0.04495835 -0.102484021    0.1874464 '0.045\n[-0.102, 0.187]'  0.0055893306 -0.005114681  0.016378624
4      CER8_Ambient      D2_Limited PBI16C009-0C-0N-0N-010N-14N 0.12736297 -0.023201437    0.2714393 '0.127\n[-0.023, 0.271]'  0.0055893306 -0.005114681  0.016378624
5         CER7_Heat    D1_Unlimited PBI17N002-0C-0N-0N-010N-16N 0.18433032  0.032223940    0.3429379  '0.184\n[0.032, 0.343]' -0.0017522679 -0.012852750  0.008808669
6      CER8_Ambient    D1_Unlimited PBI17N002-0C-0N-0N-010N-16N 0.15995677  0.007042050    0.3135489   '0.16\n[0.007, 0.314]' -0.0017522679 -0.012852750  0.008808669
7         CER7_Heat      D2_Limited PBI17N002-0C-0N-0N-010N-16N 0.09554651 -0.061785945    0.2463861 '0.096\n[-0.062, 0.246]' -0.0007683233 -0.012441540  0.010304793
8      CER8_Ambient      D2_Limited PBI17N002-0C-0N-0N-010N-16N 0.05003035 -0.100548690    0.2073912  '0.05\n[-0.101, 0.207]' -0.0007683233 -0.012441540  0.010304793
9         CER7_Heat    D1_Unlimited         PBI17N006-0N-0N-89N 0.41660011  0.264659892    0.5674929  '0.417\n[0.265, 0.567]' -0.0019451914 -0.012753209  0.008856550
10     CER8_Ambient    D1_Unlimited         PBI17N006-0N-0N-89N 0.40830001  0.251336357    0.5570386  '0.408\n[0.251, 0.557]' -0.0019451914 -0.012753209  0.008856550
11        CER7_Heat      D2_Limited         PBI17N006-0N-0N-89N 0.13764445 -0.013486463    0.2961097 '0.138\n[-0.013, 0.296]'  0.0015787944 -0.009899489  0.012649374
12     CER8_Ambient      D2_Limited         PBI17N006-0N-0N-89N 0.10136644 -0.056595466    0.2529579 '0.101\n[-0.057, 0.253]'  0.0015787944 -0.009899489  0.012649374
13        CER7_Heat    D1_Unlimited  PBI17N007-0C-0N-0N-010N-8N 0.37970102  0.229958524    0.5329856    '0.38\n[0.23, 0.533]'  0.0057784145 -0.004869020  0.016836400
14     CER8_Ambient    D1_Unlimited  PBI17N007-0C-0N-0N-010N-8N 0.26744964  0.118268863    0.4255467  '0.267\n[0.118, 0.426]'  0.0057784145 -0.004869020  0.016836400
15        CER7_Heat      D2_Limited  PBI17N007-0C-0N-0N-010N-8N 0.06696552 -0.088160484    0.2244000 '0.067\n[-0.088, 0.224]'  0.0025015957 -0.008836021  0.013963010
16     CER8_Ambient      D2_Limited  PBI17N007-0C-0N-0N-010N-8N 0.01848235 -0.137019372    0.1689897 '0.018\n[-0.137, 0.169]'  0.0025015957 -0.008836021  0.013963010
17        CER7_Heat    D1_Unlimited  PBI17N013-0C-0N-0N-010N-8N 0.38050194  0.224512634    0.5303086   '0.381\n[0.225, 0.53]' -0.0033683318 -0.013904764  0.007804901
18     CER8_Ambient    D1_Unlimited  PBI17N013-0C-0N-0N-010N-8N 0.40035052  0.249692370    0.5535988     '0.4\n[0.25, 0.554]' -0.0033683318 -0.013904764  0.007804901
19        CER7_Heat      D2_Limited  PBI17N013-0C-0N-0N-010N-8N 0.15174256 -0.018206969    0.3296940  '0.152\n[-0.018, 0.33]'  0.0024020497 -0.009443417  0.014101201
20     CER8_Ambient      D2_Limited  PBI17N013-0C-0N-0N-010N-8N 0.19524012  0.040076195    0.3512786   '0.195\n[0.04, 0.351]'  0.0024020497 -0.009443417  0.014101201
21        CER7_Heat    D1_Unlimited PBI17N016-0C-0N-0N-010N-27N 0.15945749  0.007597480    0.3085600  '0.159\n[0.008, 0.309]' -0.0018105312 -0.012611888  0.009151115
22     CER8_Ambient    D1_Unlimited PBI17N016-0C-0N-0N-010N-27N 0.19533789  0.037235064    0.3459520  '0.195\n[0.037, 0.346]' -0.0018105312 -0.012611888  0.009151115
23        CER7_Heat      D2_Limited PBI17N016-0C-0N-0N-010N-27N 0.06356841 -0.092069626    0.2213411 '0.064\n[-0.092, 0.221]' -0.0022612157 -0.013994911  0.008976878
24     CER8_Ambient      D2_Limited PBI17N016-0C-0N-0N-010N-27N 0.16779927  0.008047408    0.3199940   '0.168\n[0.008, 0.32]' -0.0022612157 -0.013994911  0.008976878
25        CER7_Heat    D1_Unlimited             PBI19N007-0N-3N 0.48483313  0.334777500    0.6366735  '0.485\n[0.335, 0.637]'  0.0031731543 -0.007668239  0.014041480
26     CER8_Ambient    D1_Unlimited             PBI19N007-0N-3N 0.54381455  0.390866194    0.6957003  '0.544\n[0.391, 0.696]'  0.0031731543 -0.007668239  0.014041480
27        CER7_Heat      D2_Limited             PBI19N007-0N-3N 0.26318417  0.110955720    0.4218045  '0.263\n[0.111, 0.422]' -0.0018189490 -0.013166687  0.009589659
28     CER8_Ambient      D2_Limited             PBI19N007-0N-3N 0.12608999 -0.029738993    0.2808076  '0.126\n[-0.03, 0.281]' -0.0018189490 -0.013166687  0.009589659
29        CER7_Heat    D1_Unlimited PBIC15034-0C-24N-010N-1N-0N 0.39512486  0.244535075    0.5536881  '0.395\n[0.245, 0.554]' -0.0040411727 -0.014711012  0.006942929
30     CER8_Ambient    D1_Unlimited PBIC15034-0C-24N-010N-1N-0N 0.41611186  0.262706998    0.5719578  '0.416\n[0.263, 0.572]' -0.0040411727 -0.014711012  0.006942929
31        CER7_Heat      D2_Limited PBIC15034-0C-24N-010N-1N-0N 0.23729447  0.078996619    0.3916621  '0.237\n[0.079, 0.392]' -0.0010129309 -0.012706840  0.009906538
32     CER8_Ambient      D2_Limited PBIC15034-0C-24N-010N-1N-0N 0.11663296 -0.041826661    0.2701649  '0.117\n[-0.042, 0.27]' -0.0010129309 -0.012706840  0.009906538
33        CER7_Heat    D1_Unlimited                    ROCKSTAR 0.40759140  0.258949683    0.5606586  '0.408\n[0.259, 0.561]'  0.0007450771 -0.009879235  0.011774663
34     CER8_Ambient    D1_Unlimited                    ROCKSTAR 0.38888940  0.235814290    0.5436577  '0.389\n[0.236, 0.544]'  0.0007450771 -0.009879235  0.011774663
35        CER7_Heat      D2_Limited                    ROCKSTAR 0.25007818  0.096043451    0.4023250   '0.25\n[0.096, 0.402]' -0.0014015058 -0.013208428  0.009635498
36     CER8_Ambient      D2_Limited                    ROCKSTAR 0.32576657  0.170348427    0.4820663   '0.326\n[0.17, 0.482]' -0.0014015058 -0.013208428  0.009635498
37        CER7_Heat    D1_Unlimited                   SUNCHASER 0.47747022  0.320942564    0.6276430  '0.477\n[0.321, 0.628]'  0.0023082835 -0.008516641  0.013218958
38     CER8_Ambient    D1_Unlimited                   SUNCHASER 0.48089341  0.326268944    0.6349031  '0.481\n[0.326, 0.635]'  0.0023082835 -0.008516641  0.013218958
39        CER7_Heat      D2_Limited                   SUNCHASER 0.19232587  0.036263572    0.3547791  '0.192\n[0.036, 0.355]' -0.0034568309 -0.014693424  0.007900996
40     CER8_Ambient      D2_Limited                   SUNCHASER 0.19279746  0.041880604    0.3541100  '0.193\n[0.042, 0.354]' -0.0034568309 -0.014693424  0.007900996
                 g1_tabquant k_emmean k_lower.HPD k_upper.HPD              k_tabquant lambda_emmean lambda_lower.HPD lambda_upper.HPD Slmax_emmean Slmax_lower.HPD
1    '0.001\n[-0.01, 0.011]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.01763447     3.016218e-08       0.07411054  -0.83072925      -1.1444408
2    '0.001\n[-0.01, 0.011]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.01763447     3.016218e-08       0.07411054  -0.70820224      -0.8977364
3   '0.006\n[-0.005, 0.016]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.04680341     1.824633e-03       0.13114073  -0.09193003      -0.4364561
4   '0.006\n[-0.005, 0.016]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.04680341     1.824633e-03       0.13114073  -0.16319288      -0.3587307
5  '-0.002\n[-0.013, 0.009]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.01763447     3.016218e-08       0.07411054  -0.43647790      -0.8020240
6  '-0.002\n[-0.013, 0.009]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.01763447     3.016218e-08       0.07411054  -0.21600726      -0.4255437
7   '-0.001\n[-0.012, 0.01]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.04680341     1.824633e-03       0.13114073  -0.22631177      -0.5862761
8   '-0.001\n[-0.012, 0.01]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.04680341     1.824633e-03       0.13114073  -0.06772304      -0.2751438
9  '-0.002\n[-0.013, 0.009]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.01763447     3.016218e-08       0.07411054  -0.98050440      -1.3439437
10 '-0.002\n[-0.013, 0.009]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.01763447     3.016218e-08       0.07411054  -0.54773826      -0.7560986
11   '0.002\n[-0.01, 0.013]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.04680341     1.824633e-03       0.13114073  -0.31913519      -0.6900884
12   '0.002\n[-0.01, 0.013]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.04680341     1.824633e-03       0.13114073  -0.13301979      -0.3418246
13  '0.006\n[-0.005, 0.017]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.01763447     3.016218e-08       0.07411054  -0.87722331      -1.2371635
14  '0.006\n[-0.005, 0.017]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.01763447     3.016218e-08       0.07411054  -0.34936268      -0.5616023
15  '0.003\n[-0.009, 0.014]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.04680341     1.824633e-03       0.13114073  -0.15137892      -0.5161741
16  '0.003\n[-0.009, 0.014]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.04680341     1.824633e-03       0.13114073  -0.02144389      -0.2256498
17 '-0.003\n[-0.014, 0.008]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.01763447     3.016218e-08       0.07411054  -0.89822103      -1.2604829
18 '-0.003\n[-0.014, 0.008]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.01763447     3.016218e-08       0.07411054  -0.53958771      -0.7432511
19  '0.002\n[-0.009, 0.014]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.04680341     1.824633e-03       0.13114073  -0.35066318      -0.7614849
20  '0.002\n[-0.009, 0.014]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.04680341     1.824633e-03       0.13114073  -0.25828971      -0.4710955
21 '-0.002\n[-0.013, 0.009]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.01763447     3.016218e-08       0.07411054  -0.37851041      -0.7338869
22 '-0.002\n[-0.013, 0.009]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.01763447     3.016218e-08       0.07411054  -0.26356249      -0.4764488
23 '-0.002\n[-0.014, 0.009]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.04680341     1.824633e-03       0.13114073  -0.15429361      -0.5314234
24 '-0.002\n[-0.014, 0.009]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.04680341     1.824633e-03       0.13114073  -0.22727272      -0.4377870
25  '0.003\n[-0.008, 0.014]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.01763447     3.016218e-08       0.07411054  -1.12886886      -1.4921657
26  '0.003\n[-0.008, 0.014]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.01763447     3.016218e-08       0.07411054  -0.72239194      -0.9255836
27  '-0.002\n[-0.013, 0.01]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.04680341     1.824633e-03       0.13114073  -0.62110761      -0.9847026
28  '-0.002\n[-0.013, 0.01]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.04680341     1.824633e-03       0.13114073  -0.17079443      -0.3834329
29 '-0.004\n[-0.015, 0.007]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.01763447     3.016218e-08       0.07411054  -0.93519016      -1.3106717
30 '-0.004\n[-0.015, 0.007]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.01763447     3.016218e-08       0.07411054  -0.56132971      -0.7701752
31  '-0.001\n[-0.013, 0.01]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.04680341     1.824633e-03       0.13114073  -0.55887118      -0.9339616
32  '-0.001\n[-0.013, 0.01]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.04680341     1.824633e-03       0.13114073  -0.15711366      -0.3646493
33   '0.001\n[-0.01, 0.012]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.01763447     3.016218e-08       0.07411054  -0.95307559      -1.3259209
34   '0.001\n[-0.01, 0.012]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.01763447     3.016218e-08       0.07411054  -0.51818680      -0.7284858
35  '-0.001\n[-0.013, 0.01]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.04680341     1.824633e-03       0.13114073  -0.58906178      -0.9430422
36  '-0.001\n[-0.013, 0.01]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.04680341     1.824633e-03       0.13114073  -0.43774545      -0.6479149
37  '0.002\n[-0.009, 0.013]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.01763447     3.016218e-08       0.07411054  -1.11304568      -1.4764472
38  '0.002\n[-0.009, 0.013]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.01763447     3.016218e-08       0.07411054  -0.63956704      -0.8480032
39 '-0.003\n[-0.015, 0.008]' 2.615539    2.538767    2.695345 '2.616\n[2.539, 2.695]'    0.04680341     1.824633e-03       0.13114073  -0.45926641      -0.8278778
40 '-0.003\n[-0.015, 0.008]' 4.587102    4.476275    4.693512 '4.587\n[4.476, 4.694]'    0.04680341     1.824633e-03       0.13114073  -0.26176705      -0.4738498
   Slmax_upper.HPD             Slmax_tabquant
1      -0.51642760 '-0.831\n[-1.144, -0.516]'
2      -0.51697670 '-0.708\n[-0.898, -0.517]'
3       0.24704882  '-0.092\n[-0.436, 0.247]'
4       0.03570141  '-0.163\n[-0.359, 0.036]'
5      -0.07051462 '-0.436\n[-0.802, -0.071]'
6      -0.01492542 '-0.216\n[-0.426, -0.015]'
7       0.14064003  '-0.226\n[-0.586, 0.141]'
8       0.13640023  '-0.068\n[-0.275, 0.136]'
9      -0.62641438 '-0.981\n[-1.344, -0.626]'
10     -0.34532991 '-0.548\n[-0.756, -0.345]'
11      0.04053865   '-0.319\n[-0.69, 0.041]'
12      0.07482503  '-0.133\n[-0.342, 0.075]'
13     -0.52064247 '-0.877\n[-1.237, -0.521]'
14     -0.14804415 '-0.349\n[-0.562, -0.148]'
15      0.22060557  '-0.151\n[-0.516, 0.221]'
16      0.18473373  '-0.021\n[-0.226, 0.185]'
17     -0.53904857  '-0.898\n[-1.26, -0.539]'
18     -0.33475711  '-0.54\n[-0.743, -0.335]'
19      0.05764602  '-0.351\n[-0.761, 0.058]'
20     -0.05263493 '-0.258\n[-0.471, -0.053]'
21     -0.02173523 '-0.379\n[-0.734, -0.022]'
22     -0.06094785 '-0.264\n[-0.476, -0.061]'
23      0.20622966  '-0.154\n[-0.531, 0.206]'
24     -0.01884647 '-0.227\n[-0.438, -0.019]'
25     -0.77833537 '-1.129\n[-1.492, -0.778]'
26     -0.51389451 '-0.722\n[-0.926, -0.514]'
27     -0.25081685 '-0.621\n[-0.985, -0.251]'
28      0.03412901  '-0.171\n[-0.383, 0.034]'
29     -0.57675825 '-0.935\n[-1.311, -0.577]'
30     -0.35330939  '-0.561\n[-0.77, -0.353]'
31     -0.19582710 '-0.559\n[-0.934, -0.196]'
32      0.05247397  '-0.157\n[-0.365, 0.052]'
33     -0.60959417  '-0.953\n[-1.326, -0.61]'
34     -0.31408474 '-0.518\n[-0.728, -0.314]'
35     -0.21970749  '-0.589\n[-0.943, -0.22]'
36     -0.22803571 '-0.438\n[-0.648, -0.228]'
37     -0.74530851 '-1.113\n[-1.476, -0.745]'
38     -0.43222093  '-0.64\n[-0.848, -0.432]'
39     -0.07518326 '-0.459\n[-0.828, -0.075]'
40     -0.05571012 '-0.262\n[-0.474, -0.056]'
```

### Format this in a nicer way for export to Excel
Keep the columns in the same order as Li6800 results and order first by parameter and then Line


``` r
li600.d.export <- li600.d.param %>% select(Room_Temperature, Device_Watering, Line, 
                                             grep("tabquant", colnames(.))) %>%
  pivot_longer(cols = 4:ncol(.),
               names_to = "params",
               values_to = "values") %>%
  pivot_wider(names_from = c("Room_Temperature", "Device_Watering"),
              values_from = "values")
li600.d.export$params <- factor(li600.d.export$params, levels = c("g0_tabquant",
                                                                  "g1_tabquant",
                                                                  "k_tabquant",
                                                                  "lambda_tabquant",
                                                                  "Slmax_tabquant"))
li600.d.export <- li600.d.export %>% arrange(params, Line) %>%
  select("params", "Line", "CER7_Heat_D2_Limited", "CER8_Ambient_D2_Limited", 
         "CER7_Heat_D1_Unlimited", "CER8_Ambient_D1_Unlimited")
li600.d.export
```

```
# A tibble: 40 × 6
   params      Line                        CER7_Heat_D2_Limited       CER8_Ambient_D2_Limited    CER7_Heat_D1_Unlimited    CER8_Ambient_D1_Unlimited
   <fct>       <fct>                       <chr>                      <chr>                      <chr>                     <chr>                    
 1 g0_tabquant PBI16C009-0C-0N-0N-010N-14N "'0.045\n[-0.102, 0.187]'" "'0.127\n[-0.023, 0.271]'" "'0.355\n[0.223, 0.488]'" "'0.531\n[0.388, 0.669]'"
 2 g0_tabquant PBI17N002-0C-0N-0N-010N-16N "'0.096\n[-0.062, 0.246]'" "'0.05\n[-0.101, 0.207]'"  "'0.184\n[0.032, 0.343]'" "'0.16\n[0.007, 0.314]'" 
 3 g0_tabquant PBI17N006-0N-0N-89N         "'0.138\n[-0.013, 0.296]'" "'0.101\n[-0.057, 0.253]'" "'0.417\n[0.265, 0.567]'" "'0.408\n[0.251, 0.557]'"
 4 g0_tabquant PBI17N007-0C-0N-0N-010N-8N  "'0.067\n[-0.088, 0.224]'" "'0.018\n[-0.137, 0.169]'" "'0.38\n[0.23, 0.533]'"   "'0.267\n[0.118, 0.426]'"
 5 g0_tabquant PBI17N013-0C-0N-0N-010N-8N  "'0.152\n[-0.018, 0.33]'"  "'0.195\n[0.04, 0.351]'"   "'0.381\n[0.225, 0.53]'"  "'0.4\n[0.25, 0.554]'"   
 6 g0_tabquant PBI17N016-0C-0N-0N-010N-27N "'0.064\n[-0.092, 0.221]'" "'0.168\n[0.008, 0.32]'"   "'0.159\n[0.008, 0.309]'" "'0.195\n[0.037, 0.346]'"
 7 g0_tabquant PBI19N007-0N-3N             "'0.263\n[0.111, 0.422]'"  "'0.126\n[-0.03, 0.281]'"  "'0.485\n[0.335, 0.637]'" "'0.544\n[0.391, 0.696]'"
 8 g0_tabquant PBIC15034-0C-24N-010N-1N-0N "'0.237\n[0.079, 0.392]'"  "'0.117\n[-0.042, 0.27]'"  "'0.395\n[0.245, 0.554]'" "'0.416\n[0.263, 0.572]'"
 9 g0_tabquant ROCKSTAR                    "'0.25\n[0.096, 0.402]'"   "'0.326\n[0.17, 0.482]'"   "'0.408\n[0.259, 0.561]'" "'0.389\n[0.236, 0.544]'"
10 g0_tabquant SUNCHASER                   "'0.192\n[0.036, 0.355]'"  "'0.193\n[0.042, 0.354]'"  "'0.477\n[0.321, 0.628]'" "'0.481\n[0.326, 0.635]'"
# ℹ 30 more rows
```

### Save data and export


``` r
save.image(saveFile)
write.csv(li600.d.export, "../data/li600.d.param.export.csv", row.names = FALSE,
          quote = TRUE)
```

