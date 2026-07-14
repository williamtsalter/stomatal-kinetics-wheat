#'# Wheat stomatal kinetics (heat & drought)
#'## --- Ed Chaplin (USyd), *Apr--June 2025* ---
#'### ANATOMICAL DATA - IMPORT  
#'
#' * Importing anatomical data
#' * Some data is collected on the individual stomata level, while other traits are aggregated to a pot level
#'
#'### Load the libraries
library(knitr)
#knitr::spin("anat.import.v1.r")
library(tidyverse)
library(dae)

#'### Set up constants and files
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, fig.path = "figures/", 
               "fig.width" = 7.5, "fig.height" = 10, "comment"=NA, warning = FALSE)
options(width = 90)

#'### Load design terms (Cell.ID, spatial/design factors, WU.Total) needed for the analysis
#' Previously taken from the external water-use object `../../data/indv.v2.Rdata`.
#' To make this repository self-contained, `indv.dat` is now rebuilt from the bundled
#' `../data/indv.design.csv` via `load.indv.design.R` (identical contents, types and
#' factor level orders). To revert to the original object, replace the source() call with:
#'   load("../../data/indv.v2.Rdata")
source("load.indv.design.R")


importFile <- "../data/Stomatal Anatomy Data.xlsx"
saveFile <- "../data/anat.dat.v1.Rdata"

#'### Anatomical traits
#'
#' * First import the traits that were measured at the pot level
#' * This sheet contains the treatment and rep info
#' * Use that sheet to then import the traits measured at the stomata level
#'

responses.anat <- c("StomatalCount", "StomatalDensity", "StomatalArea", 
                    "StomatalWidth", "StomatalLength", "PoreArea", "gsmax")
respunits.anat <- c("", "per mm^2", "microm^2", "microm","microm", "microm", "mol m^-2 s^-1")

#'#### Import pot level data
anat.pot.dat <- readxl::read_xlsx(importFile, "Processed Stomatal anatomy Data") %>% as.data.frame()
names(anat.pot.dat)[c(7, 8, 9, 11, 13, 15, 17)] <- responses.anat
summary(anat.pot.dat)

#'#### Import stomata level data
anat.stoma.dat <- readxl::read_xlsx(importFile, "Raw All Output Anatomy Data") %>% as.data.frame()
names(anat.stoma.dat)[c(2:5, 8)] <- c("Stomata", responses.anat[3:6])
summary(anat.stoma.dat)

##incorporate treatment info into anat.stoma.dat
anat.stoma.dat <- left_join(anat.stoma.dat, 
                            anat.pot.dat %>% select(Filename, Line, Temperature,
                                                    Watering, Rep, 
                                                    Line_Temperature_Watering_Rep))


#'#### Define variables to allow for appending design info from indv.dat
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
anat.stoma.dat$Stomata <- factor(anat.stoma.dat$Stomata)

anat.stoma.dat <- left_join(anat.stoma.dat %>% select(-Line), indv.dat,
                          by = c("Watering", "Replicate", "Rm_Temperature", "Line.ID"))

#'### Finalise data
anat.pot.dat <- anat.pot.dat[order(anat.pot.dat$Cell.ID),]
summary(anat.pot.dat)

anat.stoma.dat <- anat.stoma.dat[order(anat.stoma.dat$Cell.ID),]
summary(anat.stoma.dat)

#'### Save data
save(anat.pot.dat, file="../data/anat.pot.dat.v1.rda")
save(anat.stoma.dat, file="../data/anat.stoma.dat.v1.rda")
save.image(saveFile)
