# Wheat stomatal kinetics (heat & drought)
## --- Ed Chaplin (USyd), *Apr--June 2025* ---
### Li-600 STOMATAL CONDUCTANCE MODEL PARAMETERS VS ANATOMICAL TRAIT CORRELATIONS

* Using the parameters of interest estimated from the stomatal conductance models and the anatomical trait predictions, calculate correlations on a 'environment' basis

### Load the libraries


``` r
library(knitr)
#knitr::spin("porometer.stomatalCondvsAnatomicalTraits.v1.r")
library(ggplot2)
```

```
## Keep up to date with changes at https://tidyverse.org/blog/
```

``` r
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, 
               "fig.width" = 7.5, "fig.height" = 10, "comment"=NA, fig.path = "figures/increase/")
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

``` r
library(tidyverse)
```

```
## ── Attaching core tidyverse packages ──────────────────────────────────────────────────────── tidyverse 2.0.0 ──
## ✔ forcats   1.0.1     ✔ stringr   1.6.0
## ✔ lubridate 1.9.5     ✔ tibble    3.3.1
## ✔ purrr     1.2.2     ✔ tidyr     1.3.2
## ✔ readr     2.2.0
```

```
## ── Conflicts ────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

``` r
library(lattice)
library(GGally)
```

```
## Warning: package 'GGally' was built under R version 4.6.1
```

### Set up constants and files


``` r
options(width=90)
load("../data/stomatalCondParamEstimation_porometer_i.v1.Rdata")
load("../data/stomatalCondParamEstimation_porometer_d.v1.Rdata")
load("../../../Anatomical traits/data/anat.diffs.full.rda")
saveFile <- "../data/stomatalCondvsAnatTraits.v1.Rdata"
```

### Add variables to make the dataframes compatible


``` r
Diffs.full <- within(Diffs.full, {
  Room_Temperature <- as.character(Rm_Temperature)
  Room_Temperature[Rm_Temperature == "7,Heat"] <- "CER7_Heat"
  Room_Temperature[Rm_Temperature == "8,Ambient"] <- "CER8_Ambient"
  
  Device_Watering <- as.character(Watering)
  Device_Watering[Watering == "Limited"] <- "D2_Limited"
  Device_Watering[Watering == "Unlimited"] <- "D1_Unlimited"
})

##do some work to li600.i.param and li600.d.param to make the variable names unique
li600.i.param <- li600.i.param %>% rename(g0_i = g0_emmean,
                                          g1_i = g1_emmean,
                                          k_i = k_emmean,
                                          lambda_i = lambda_emmean,
                                          Slmax_i = Slmax_emmean)

li600.d.param <- li600.d.param %>% rename(g0_d = g0_emmean,
                                          g1_d = g1_emmean,
                                          k_d = k_emmean,
                                          lambda_d = lambda_emmean,
                                          Slmax_d = Slmax_emmean)

##convert Diffs.full to wide format and then merge with li600.i.param and li600.d.param
traits.df <- Diffs.full %>% select(Room_Temperature, Device_Watering, Line, 
                                   Trait, predicted.value) %>%
  pivot_wider(names_from = c("Trait"), values_from = "predicted.value") %>%
  left_join(., li600.i.param %>% select(Room_Temperature, Device_Watering, Line,
                                        g0_i, g1_i, k_i,
                                        lambda_i, Slmax_i)) %>%
  left_join(., li600.d.param %>% select(Room_Temperature, Device_Watering, Line,
                                        g0_d, g1_d, k_d,
                                        lambda_d, Slmax_d)) %>%
  as.data.frame()
```

```
Joining with `by = join_by(Room_Temperature, Device_Watering, Line)`
Joining with `by = join_by(Room_Temperature, Device_Watering, Line)`
```

``` r
head(traits.df)
```

```
  Room_Temperature Device_Watering                        Line StomatalCount
1        CER7_Heat      D2_Limited                   SUNCHASER      24.25000
2        CER7_Heat      D2_Limited                    ROCKSTAR      24.50000
3        CER7_Heat      D2_Limited         PBI17N006-0N-0N-89N      26.00000
4        CER7_Heat      D2_Limited  PBI17N013-0C-0N-0N-010N-8N      29.00902
5        CER7_Heat      D2_Limited PBIC15034-0C-24N-010N-1N-0N      26.25000
6        CER7_Heat      D2_Limited PBI16C009-0C-0N-0N-010N-14N      29.50000
  StomatalDensity StomatalArea StomatalWidth StomatalLength PoreArea    gsmax       g0_i
1        47.14864     1572.885      31.83217       64.17884 961.0736 1.276378 0.08794347
2        47.63471     1810.276      36.37408       64.72257 978.2416 1.299678 0.02812548
3        50.55112     1643.381      33.55433       63.49741 942.4175 1.351104 0.03065694
4        56.40148     1776.963      37.51749       61.38072 880.1979 1.451460 0.06633263
5        51.03719     1520.344      30.71767       64.23318 967.4540 1.378773 0.07179241
6        57.35608     1498.850      33.45667       57.28076 777.7226 1.457145 0.04117321
        g1_i      k_i lambda_i   Slmax_i       g0_d         g1_d      k_d   lambda_d
1 0.36222523 5.952203 5.482128 0.2815461 0.19232587 -0.003456831 2.615539 0.04680341
2 0.36250003 5.120474 5.482128 0.3996154 0.25007818 -0.001401506 2.615539 0.04680341
3 0.41342218 3.146440 5.482128 0.7455489 0.13764445  0.001578794 2.615539 0.04680341
4 0.33471686 2.164836 5.482128 0.7588531 0.15174256  0.002402050 2.615539 0.04680341
5 0.26972436 7.207846 5.482128 0.1677780 0.23729447 -0.001012931 2.615539 0.04680341
6 0.09041874 2.018804 5.482128 0.1505521 0.04495835  0.005589331 2.615539 0.04680341
      Slmax_d
1 -0.45926641
2 -0.58906178
3 -0.31913519
4 -0.35066318
5 -0.55887118
6 -0.09193003
```

``` r
pairs(traits.df[,4:ncol(traits.df)])
```

![plot of chunk unnamed-chunk-3](figures/increase/unnamed-chunk-3-1.png)

### Put together proper plot


``` r
traits <- c("g0_i", "g1_i", "k_i", "Slmax_i", "g0_d", "g1_d", "Slmax_d", #"k_d"
            "StomatalDensity", "StomatalArea", "StomatalWidth", "StomatalLength",
            "gsmax")

trait_labs <- list(
  g0_i = expression(italic(r)[0[i]]),
  g1_i = expression(italic(g)[smax]),
  k_i = expression(italic(k)[i]),
  Slmax_i = expression(italic(Sl)[max[i]]),
  g0_d = expression(italic(r)[0[d]]),
  g1_d = expression(italic(g)[smin]),
 # k_d = expression(italic(k)[d]),
  Slmax_d = expression(italic(Sl)[max[d]]),
  StomatalDensity = "SD",
  StomatalArea = "SA",
  StomatalWidth = "GCW",
  StomatalLength = "GCL",
  gsmax = expression(g[smax[a]])
)


diag_lab <- function(data, mapping) {
  var <- as_label(mapping$x)
  GGally::ggally_text(
    label = trait_labs[[var]],
    xP = 0.5,
    yP = 0.5,
    size = 6
  ) + theme_void()
}

lower_scatter <- function(data, mapping) {
  
  ggplot(data = data, mapping = mapping) +
    geom_point(
      shape = 16,
      colour = "black",
      alpha = 1,
      size = 1.5
    ) +
    theme_bw()
}


p1 <- ggpairs(
  filter(traits.df, Room_Temperature == "CER7_Heat" &
           Device_Watering == "D1_Unlimited"),
  columns = traits,
  upper = list(continuous = wrap("cor", stats = TRUE, 
                                 digits = 2, size = 4,
                                 colour = "black")),
  lower = list(continuous = lower_scatter),
  diag = list(continuous = diag_lab),
  progress = FALSE) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )
p1 <- p1 + ggtitle("CER7_Heat, D1_Unlimited")
p1
```

```
Warning in geom_text(label = label, mapping = mapping, ...): Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
```

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-1.png)

``` r
ggsave("figures/CER7_Heat_D1_Unlimited_Correlation.png", width = 10, height = 10)
```

```
Warning in geom_text(label = label, mapping = mapping, ...): Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
```

``` r
ggsave("figures/CER7_Heat_D1_Unlimited_Correlation.pdf", width = 10, height = 10)
```

```
Warning in geom_text(label = label, mapping = mapping, ...): Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
```

``` r
ggsave("figures/CER7_Heat_D1_Unlimited_Correlation.tiff", width = 10, height = 10)
```

```
Warning in geom_text(label = label, mapping = mapping, ...): Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
```

``` r
p2 <- ggpairs(
  filter(traits.df, Room_Temperature == "CER7_Heat" &
           Device_Watering == "D2_Limited"),
  columns = traits,
  upper = list(continuous = wrap("cor", stats = TRUE, 
                                 digits = 2, size = 4,
                                 colour = "black")),
  lower = list(continuous = lower_scatter),
  diag = list(continuous = diag_lab),
  progress = FALSE) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )
p2 <- p2 + ggtitle("CER7_Heat, D2_Limited")
p2
```

```
Warning in geom_text(label = label, mapping = mapping, ...): Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
```

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-2.png)

``` r
ggsave("figures/CER7_Heat_D2_Limited_Correlation.png", width = 10, height = 10)
```

```
Warning in geom_text(label = label, mapping = mapping, ...): Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
```

``` r
ggsave("figures/CER7_Heat_D2_Limited_Correlation.pdf", width = 10, height = 10)
```

```
Warning in geom_text(label = label, mapping = mapping, ...): Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
```

``` r
ggsave("figures/CER7_Heat_D2_Limited_Correlation.tiff", width = 10, height = 10)
```

```
Warning in geom_text(label = label, mapping = mapping, ...): Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
```

``` r
p3 <- ggpairs(
  filter(traits.df, Room_Temperature == "CER8_Ambient" &
           Device_Watering == "D1_Unlimited"),
  columns = traits,
  upper = list(continuous = wrap("cor", stats = TRUE, 
                                 digits = 2, size = 4,
                                 colour = "black")),
  lower = list(continuous = lower_scatter),
  diag = list(continuous = diag_lab),
  progress = FALSE) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )
p3 <- p3 + ggtitle("CER8_Ambient, D1_Unlimited")
p3
```

```
Warning in geom_text(label = label, mapping = mapping, ...): Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
```

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-3.png)

``` r
ggsave("figures/CER8_Ambient_D1_Unlimited_Correlation.png", width = 10, height = 10)
```

```
Warning in geom_text(label = label, mapping = mapping, ...): Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
```

``` r
ggsave("figures/CER8_Ambient_D1_Unlimited_Correlation.pdf", width = 10, height = 10)
```

```
Warning in geom_text(label = label, mapping = mapping, ...): Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
```

``` r
ggsave("figures/CER8_Ambient_D1_Unlimited_Correlation.tiff", width = 10, height = 10)
```

```
Warning in geom_text(label = label, mapping = mapping, ...): Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
```

``` r
p4 <- ggpairs(
  filter(traits.df, Room_Temperature == "CER8_Ambient" &
           Device_Watering == "D2_Limited"),
  columns = traits,
  upper = list(continuous = wrap("cor", stats = TRUE, 
                                 digits = 2, size = 4,
                                 colour = "black")),
  lower = list(continuous = lower_scatter),
  diag = list(continuous = diag_lab),
  progress = FALSE) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )
p4 <- p4 + ggtitle("CER8_Ambient, D2_Limited")
p4
```

```
Warning in geom_text(label = label, mapping = mapping, ...): Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
```

![plot of chunk unnamed-chunk-4](figures/increase/unnamed-chunk-4-4.png)

``` r
ggsave("figures/CER8_Ambient_D2_Limited_Correlation.png", width = 10, height = 10)
```

```
Warning in geom_text(label = label, mapping = mapping, ...): Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
```

``` r
ggsave("figures/CER8_Ambient_D2_Limited_Correlation.pdf", width = 10, height = 10)
```

```
Warning in geom_text(label = label, mapping = mapping, ...): Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
```

``` r
ggsave("figures/CER8_Ambient_D2_Limited_Correlation.tiff", width = 10, height = 10)
```

```
Warning in geom_text(label = label, mapping = mapping, ...): Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
Ignoring unknown parameters: `stats`
```

### Save data and export


``` r
save.image(saveFile)
```

