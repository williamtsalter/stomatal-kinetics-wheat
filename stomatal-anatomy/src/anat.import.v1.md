---
output:
  pdf_document: default
  html_document: default
---
# Wheat stomatal kinetics (heat & drought)
## --- Ed Chaplin (USyd), *Apr--June 2025* ---
### ANATOMICAL DATA - IMPORT  

* Importing anatomical data
* Some data is collected on the individual stomata level, while other traits are aggregated to a pot level

### Load the libraries


``` r
library(knitr)
#knitr::spin("anat.import.v1.r")
library(tidyverse)
library(dae)
```

### Set up constants and files


``` r
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, fig.path = "figures/", 
               "fig.width" = 7.5, "fig.height" = 10, "comment"=NA, warning = FALSE)
options(width = 90)
load("../../data/indv.v2.Rdata")
importFile <- "../data/Stomatal Anatomy Data.xlsx"
saveFile <- "../data/anat.dat.v1.Rdata"
```

### Import indv.dat from WU analysis, as this will have all of the design terms I need
### Anatomical traits

* First import the traits that were measured at the pot level
* This sheet contains the treatment and rep info
* Use that sheet to then import the traits measured at the stomata level



``` r
responses.anat <- c("StomatalCount", "StomatalDensity", "StomatalArea", 
                    "StomatalWidth", "StomatalLength", "PoreArea", "gsmax")
respunits.anat <- c("", "per mm^2", "microm^2", "microm","microm", "microm", "mol m^-2 s^-1")
```

#### Import pot level data


``` r
anat.pot.dat <- readxl::read_xlsx(importFile, "Processed Stomatal anatomy Data") %>% as.data.frame()
names(anat.pot.dat)[c(7, 8, 9, 11, 13, 15, 17)] <- responses.anat
summary(anat.pot.dat)
```

```
      Filename        Line           Temperature       Watering        Rep       
 Length   :159   Min.   : 1.000   Length   :159   Length   :159   Min.   :1.000  
 N.unique :159   1st Qu.: 3.000   N.unique :  2   N.unique :  2   1st Qu.:1.500  
 N.blank  :  0   Median : 6.000   N.blank  :  0   N.blank  :  0   Median :3.000  
 Min.nchar: 42   Mean   : 5.509   Min.nchar:  4   Min.nchar:  7   Mean   :2.503  
 Max.nchar: 65   3rd Qu.: 8.000   Max.nchar:  7   Max.nchar:  9   3rd Qu.:3.500  
                 Max.   :10.000                                   Max.   :4.000  
 Line_Temperature_Watering_Rep StomatalCount   StomatalDensity  StomatalArea 
 Length   :159                 Min.   :16.00   Min.   :31.11   Min.   :1045  
 N.unique :159                 1st Qu.:23.00   1st Qu.:44.72   1st Qu.:1517  
 N.blank  :  0                 Median :27.00   Median :52.50   Median :1624  
 Min.nchar: 16                 Mean   :26.35   Mean   :51.22   Mean   :1640  
 Max.nchar: 22                 3rd Qu.:29.00   3rd Qu.:56.38   3rd Qu.:1751  
                               Max.   :41.00   Max.   :79.72   Max.   :2726  
 Stdev Stomatal Area (µm^2) StomatalWidth   Stdev Stomatal Width (µm) StomatalLength 
 Min.   :116.2              Min.   :27.94   Min.   :1.697             Min.   :48.11  
 1st Qu.:185.2              1st Qu.:31.67   1st Qu.:3.090             1st Qu.:59.95  
 Median :235.9              Median :34.09   Median :3.850             Median :62.91  
 Mean   :250.8              Mean   :34.20   Mean   :3.977             Mean   :62.36  
 3rd Qu.:300.6              3rd Qu.:35.79   3rd Qu.:4.637             3rd Qu.:65.83  
 Max.   :882.9              Max.   :50.81   Max.   :9.158             Max.   :72.90  
 Stdev Stomatal Length (µm)    PoreArea      Stdev pore area (µm^2)     gsmax       
 Min.   : 2.297             Min.   : 542.4   Min.   :  70.17        Min.   :0.8403  
 1st Qu.: 4.528             1st Qu.: 844.9   1st Qu.: 127.19        1st Qu.:1.2128  
 Median : 5.263             Median : 921.7   Median : 152.02        Median :1.3498  
 Mean   : 5.696             Mean   : 913.6   Mean   : 171.93        Mean   :1.3419  
 3rd Qu.: 6.401             3rd Qu.:1010.4   3rd Qu.: 188.63        3rd Qu.:1.4710  
 Max.   :22.882             Max.   :1240.4   Max.   :1290.05        Max.   :2.0044  
```

#### Import stomata level data


``` r
anat.stoma.dat <- readxl::read_xlsx(importFile, "Raw All Output Anatomy Data") %>% as.data.frame()
names(anat.stoma.dat)[c(2:5, 8)] <- c("Stomata", responses.anat[3:6])
summary(anat.stoma.dat)
```

```
      Filename       Stomata       StomatalArea  StomatalWidth   StomatalLength  
 Length   :4189   Min.   : 1.00   Min.   : 428   Min.   :14.85   Min.   : 31.85  
 N.unique : 159   1st Qu.: 7.00   1st Qu.:1411   1st Qu.:30.59   1st Qu.: 57.27  
 N.blank  :   0   Median :14.00   Median :1598   Median :33.41   Median : 62.29  
 Min.nchar:  42   Mean   :14.06   Mean   :1625   Mean   :34.08   Mean   : 61.95  
 Max.nchar:  65   3rd Qu.:20.00   3rd Qu.:1796   3rd Qu.:36.75   3rd Qu.: 66.48  
                  Max.   :41.00   Max.   :6075   Max.   :77.95   Max.   :193.08  
 Pore length (µm) Pore width (µm)     PoreArea      Pore depth (µm) 
 Min.   : 17.29   Min.   : 8.647   Min.   : 234.9   Min.   : 4.324  
 1st Qu.: 31.10   1st Qu.:15.548   1st Qu.: 759.4   1st Qu.: 7.774  
 Median : 33.82   Median :16.911   Median : 898.5   Median : 8.456  
 Mean   : 33.64   Mean   :16.819   Mean   : 902.4   Mean   : 8.409  
 3rd Qu.: 36.10   3rd Qu.:18.049   3rd Qu.:1023.5   3rd Qu.: 9.025  
 Max.   :104.84   Max.   :52.420   Max.   :8632.6   Max.   :26.210  
```

``` r
##incorporate treatment info into anat.stoma.dat
anat.stoma.dat <- left_join(anat.stoma.dat, 
                            anat.pot.dat %>% select(Filename, Line, Temperature,
                                                    Watering, Rep, 
                                                    Line_Temperature_Watering_Rep))
```

```
Joining with `by = join_by(Filename)`
```

#### Define variables to allow for appending design info from indv.dat


``` r
anat.pot.dat$Replicate <- factor(anat.pot.dat$Rep)
anat.pot.dat$Line.ID <- sprintf("L%02d", anat.pot.dat$Line)
anat.pot.dat$Line.ID <- factor(anat.pot.dat$Line.ID)
anat.pot.dat$Rm_Temperature <- anat.pot.dat$Temperature
anat.pot.dat$Rm_Temperature[anat.pot.dat$Temperature == "Heat"] <- "7,Heat"
anat.pot.dat$Rm_Temperature[anat.pot.dat$Temperature == "Ambient"] <- "8,Ambient"
anat.pot.dat$Rm_Temperature <- factor(anat.pot.dat$Rm_Temperature)
anat.pot.dat$Watering <- factor(anat.pot.dat$Watering)

anat.pot.dat <- left_join(anat.pot.dat %>% select(-Line), indv.dat,
                          by = c("Watering", "Replicate", "Rm_Temperature", "Line.ID"))

##do the same for the anat.stoma.dat dataframe
anat.stoma.dat$Replicate <- factor(anat.stoma.dat$Rep)
anat.stoma.dat$Line.ID <- sprintf("L%02d", anat.stoma.dat$Line)
anat.stoma.dat$Line.ID <- factor(anat.stoma.dat$Line.ID)
anat.stoma.dat$Rm_Temperature <- anat.stoma.dat$Temperature
anat.stoma.dat$Rm_Temperature[anat.stoma.dat$Temperature == "Heat"] <- "7,Heat"
anat.stoma.dat$Rm_Temperature[anat.stoma.dat$Temperature == "Ambient"] <- "8,Ambient"
anat.stoma.dat$Rm_Temperature <- factor(anat.stoma.dat$Rm_Temperature)
anat.stoma.dat$Watering <- factor(anat.stoma.dat$Watering)

anat.stoma.dat <- left_join(anat.stoma.dat %>% select(-Line), indv.dat,
                          by = c("Watering", "Replicate", "Rm_Temperature", "Line.ID"))
```

### Finalise data


``` r
anat.pot.dat <- anat.pot.dat[order(anat.pot.dat$Cell.ID),]
summary(anat.pot.dat)
```

```
      Filename      Temperature       Watering       Rep       
 Length   :159   Length   :159   Limited  :79   Min.   :1.000  
 N.unique :159   N.unique :  2   Unlimited:80   1st Qu.:1.500  
 N.blank  :  0   N.blank  :  0                  Median :3.000  
 Min.nchar: 42   Min.nchar:  4                  Mean   :2.503  
 Max.nchar: 65   Max.nchar:  7                  3rd Qu.:3.500  
                                                Max.   :4.000  
                                                               
 Line_Temperature_Watering_Rep StomatalCount   StomatalDensity  StomatalArea 
 Length   :159                 Min.   :16.00   Min.   :31.11   Min.   :1045  
 N.unique :159                 1st Qu.:23.00   1st Qu.:44.72   1st Qu.:1517  
 N.blank  :  0                 Median :27.00   Median :52.50   Median :1624  
 Min.nchar: 16                 Mean   :26.35   Mean   :51.22   Mean   :1640  
 Max.nchar: 22                 3rd Qu.:29.00   3rd Qu.:56.38   3rd Qu.:1751  
                               Max.   :41.00   Max.   :79.72   Max.   :2726  
                                                                             
 Stdev Stomatal Area (µm^2) StomatalWidth   Stdev Stomatal Width (µm) StomatalLength 
 Min.   :116.2              Min.   :27.94   Min.   :1.697             Min.   :48.11  
 1st Qu.:185.2              1st Qu.:31.67   1st Qu.:3.090             1st Qu.:59.95  
 Median :235.9              Median :34.09   Median :3.850             Median :62.91  
 Mean   :250.8              Mean   :34.20   Mean   :3.977             Mean   :62.36  
 3rd Qu.:300.6              3rd Qu.:35.79   3rd Qu.:4.637             3rd Qu.:65.83  
 Max.   :882.9              Max.   :50.81   Max.   :9.158             Max.   :72.90  
                                                                                     
 Stdev Stomatal Length (µm)    PoreArea      Stdev pore area (µm^2)     gsmax       
 Min.   : 2.297             Min.   : 542.4   Min.   :  70.17        Min.   :0.8403  
 1st Qu.: 4.528             1st Qu.: 844.9   1st Qu.: 127.19        1st Qu.:1.2128  
 Median : 5.263             Median : 921.7   Median : 152.02        Median :1.3498  
 Mean   : 5.696             Mean   : 913.6   Mean   : 171.93        Mean   :1.3419  
 3rd Qu.: 6.401             3rd Qu.:1010.4   3rd Qu.: 188.63        3rd Qu.:1.4710  
 Max.   :22.882             Max.   :1240.4   Max.   :1290.05        Max.   :2.0044  
                                                                                    
 Replicate    Line.ID     Rm_Temperature      Cell.ID          Pot.ID      Room    Lane  
 1:40      L01    :16   7,Heat   :79     Length   :159   Length   :159   CER7:79   A:28  
 2:39      L02    :16   8,Ambient:80     N.unique :159   N.unique :159   CER8:80   B:24  
 3:40      L03    :16                    N.blank  :  0   N.blank  :  0             C:27  
 4:40      L05    :16                    Min.nchar:  7   Min.nchar:  6             D:28  
           L06    :16                    Max.nchar:  7   Max.nchar:  6             E:24  
           L07    :16                                                              F:28  
           (Other):63                                                                    
    Position      PCoord    Side    Location  Block  SubUnit      Cell        Treatment  
 2      :12   1      :  8   N:79   Inner:80   1:80   1:80    1      :16   Length   :159  
 3      :12   2      :  8   S:80   Outer:79   2:79   2:79    2      :16   N.unique : 40  
 4      :12   3      :  8                                    3      :16   N.blank  :  0  
 5      :12   4      :  8                                    4      :16   Min.nchar:  3  
 6      :12   5      :  8                                    5      :16   Max.nchar:  3  
 9      :12   6      :  8                                    6      :16                  
 (Other):87   (Other):111                                    (Other):63                  
                          Line        AWatering      HWatering 
 SUNCHASER                  :16   rest     :79   rest     :80  
 ROCKSTAR                   :16   Limited  :40   Limited  :39  
 PBI17N006-0N-0N-89N        :16   Unlimited:40   Unlimited:40  
 PBIC15034-0C-24N-010N-1N-0N:16                                
 PBI16C009-0C-0N-0N-010N-14N:16                                
 PBI19N007-0N-3N            :16                                
 (Other)                    :63                                
                         ALine                            HLine                    TW    
 rest                       :79   rest                       :80   Ambient-Limited  :40  
 SUNCHASER                  : 8   SUNCHASER                  : 8   Ambient-Unlimited:40  
 ROCKSTAR                   : 8   ROCKSTAR                   : 8   Heat-Limited     :39  
 PBI17N006-0N-0N-89N        : 8   PBI17N006-0N-0N-89N        : 8   Heat-Unlimited   :40  
 PBI17N013-0C-0N-0N-010N-8N : 8   PBIC15034-0C-24N-010N-1N-0N: 8                         
 PBIC15034-0C-24N-010N-1N-0N: 8   PBI16C009-0C-0N-0N-010N-14N: 8                         
 (Other)                    :40   (Other)                    :39                         
                 WT     ASide  ABlock    WU.Total     
 Limited-Ambient  :40   N:79   1:80   Min.   :0.4881  
 Limited-Heat     :39   S:80   2:79   1st Qu.:1.5072  
 Unlimited-Ambient:40                 Median :2.3507  
 Unlimited-Heat   :40                 Mean   :2.7552  
                                      3rd Qu.:3.9597  
                                      Max.   :6.2256  
                                      NAs    :2       
```

``` r
anat.stoma.dat <- anat.stoma.dat[order(anat.stoma.dat$Cell.ID),]
summary(anat.stoma.dat)
```

```
      Filename       Stomata       StomatalArea  StomatalWidth   StomatalLength  
 Length   :4189   Min.   : 1.00   Min.   : 428   Min.   :14.85   Min.   : 31.85  
 N.unique : 159   1st Qu.: 7.00   1st Qu.:1411   1st Qu.:30.59   1st Qu.: 57.27  
 N.blank  :   0   Median :14.00   Median :1598   Median :33.41   Median : 62.29  
 Min.nchar:  42   Mean   :14.06   Mean   :1625   Mean   :34.08   Mean   : 61.95  
 Max.nchar:  65   3rd Qu.:20.00   3rd Qu.:1796   3rd Qu.:36.75   3rd Qu.: 66.48  
                  Max.   :41.00   Max.   :6075   Max.   :77.95   Max.   :193.08  
                                                                                 
 Pore length (µm) Pore width (µm)     PoreArea      Pore depth (µm)     Temperature  
 Min.   : 17.29   Min.   : 8.647   Min.   : 234.9   Min.   : 4.324   Length   :4189  
 1st Qu.: 31.10   1st Qu.:15.548   1st Qu.: 759.4   1st Qu.: 7.774   N.unique :   2  
 Median : 33.82   Median :16.911   Median : 898.5   Median : 8.456   N.blank  :   0  
 Mean   : 33.64   Mean   :16.819   Mean   : 902.4   Mean   : 8.409   Min.nchar:   4  
 3rd Qu.: 36.10   3rd Qu.:18.049   3rd Qu.:1023.5   3rd Qu.: 9.025   Max.nchar:   7  
 Max.   :104.84   Max.   :52.420   Max.   :8632.6   Max.   :26.210                   
                                                                                     
      Watering         Rep       Line_Temperature_Watering_Rep Replicate    Line.ID    
 Limited  :2115   Min.   :1.00   Length   :4189                1: 985    L06    : 488  
 Unlimited:2074   1st Qu.:2.00   N.unique : 159                2:1059    L07    : 458  
                  Median :3.00   N.blank  :   0                3:1085    L10    : 434  
                  Mean   :2.53   Min.nchar:  16                4:1060    L08    : 424  
                  3rd Qu.:4.00   Max.nchar:  22                          L04    : 417  
                  Max.   :4.00                                           L09    : 412  
                                                                         (Other):1556  
   Rm_Temperature      Cell.ID           Pot.ID       Room      Lane       Position   
 7,Heat   :2088   Length   :4189   Length   :4189   CER7:2088   A:736   11     : 340  
 8,Ambient:2101   N.unique : 159   N.unique : 159   CER8:2101   B:600   5      : 338  
                  N.blank  :   0   N.blank  :   0               C:751   10     : 325  
                  Min.nchar:   7   Min.nchar:   6               D:713   3      : 323  
                  Max.nchar:   7   Max.nchar:   6               E:660   2      : 318  
                                                                F:729   9      : 317  
                                                                        (Other):2228  
     PCoord     Side      Location    Block    SubUnit       Cell          Treatment   
 11     : 232   N:2087   Inner:2082   1:2060   1:2127   7      : 438   Length   :4189  
 9      : 228   S:2102   Outer:2107   2:2129   2:2062   3      : 436   N.unique :  40  
 5      : 225                                           5      : 435   N.blank  :   0  
 3      : 224                                           2      : 429   Min.nchar:   3  
 10     : 218                                           6      : 426   Max.nchar:   3  
 7      : 214                                           8      : 420                   
 (Other):2848                                           (Other):1605                   
                          Line          AWatering        HWatering   
 PBI16C009-0C-0N-0N-010N-14N: 488   rest     :2088   rest     :2101  
 PBI19N007-0N-3N            : 458   Limited  :1071   Limited  :1044  
 PBI17N002-0C-0N-0N-010N-16N: 434   Unlimited:1030   Unlimited:1044  
 PBI17N007-0C-0N-0N-010N-8N : 424                                    
 PBI17N013-0C-0N-0N-010N-8N : 417                                    
 PBI17N016-0C-0N-0N-010N-27N: 412                                    
 (Other)                    :1556                                    
                         ALine                              HLine     
 rest                       :2088   rest                       :2101  
 PBI16C009-0C-0N-0N-010N-14N: 231   PBI16C009-0C-0N-0N-010N-14N: 257  
 PBI17N013-0C-0N-0N-010N-8N : 218   PBI19N007-0N-3N            : 240  
 PBI19N007-0N-3N            : 218   PBI17N002-0C-0N-0N-010N-16N: 220  
 PBI17N016-0C-0N-0N-010N-27N: 215   PBI17N007-0C-0N-0N-010N-8N : 211  
 PBI17N002-0C-0N-0N-010N-16N: 214   PBI17N013-0C-0N-0N-010N-8N : 199  
 (Other)                    :1005   (Other)                    : 961  
                 TW                       WT       ASide    ABlock      WU.Total     
 Ambient-Limited  :1071   Limited-Ambient  :1071   N:2044   1:2070   Min.   :0.4881  
 Ambient-Unlimited:1030   Limited-Heat     :1044   S:2145   2:2119   1st Qu.:1.5072  
 Heat-Limited     :1044   Unlimited-Ambient:1030                     Median :2.2327  
 Heat-Unlimited   :1044   Unlimited-Heat   :1044                     Mean   :2.7279  
                                                                     3rd Qu.:3.9576  
                                                                     Max.   :6.2256  
                                                                     NAs    :49      
```

### Save data


``` r
save(anat.pot.dat, file="../data/anat.pot.dat.v1.rda")
save(anat.stoma.dat, file="../data/anat.stoma.dat.v1.rda")
save.image(saveFile)
```

