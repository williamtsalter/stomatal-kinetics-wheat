#'# Wheat stomatal kinetics (heat & drought)
#'## --- Ed Chaplin (USyd), *Apr--June 2025* ---
#'### Li-6800 DATA
#'
#' * Importing stomatal conductance data obtained from Li-6800
#' * There are individual spreadsheets for each plant sampled with the Li-6800
#' * Reps 1, 2 and 3 were used for measurements. Rep 4 was not used for Li-6800 measurements
#' * There are a total of 12 individual Li-6800 spreadsheets 
#' * Each spreadsheet contains a 'Raw Data' tab and a 'Sorted Data' tab. The 'Sorted Data' tab is the one to use for analyses.
#' * The key data parameters to assess are 'gsw' (which we can refer to as 'gs') (Column P) and 'A' (net photosynthetic carbon assimilation) (Column K)
#' * The light level (to determine if it is phase 1 dark, phase 2 light or phase 3 dark) is denoted in column DG
#'
#'### Load the libraries
library(knitr)
#knitr::spin("li6800.import.v1.r")
library(ggplot2)
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, 
               "fig.width" = 7.5, "fig.height" = 10, "comment"=NA, fig.path = "figures/")
library(dae)
library(dplyr)
library(growthPheno) # calcLagged
packageVersion("growthPheno")

#'### Set up constants and files
options(width=90)
inputFiles <- list.files("../data")
inputFiles <- inputFiles[grepl("L01_", inputFiles)]
writeFile  <- "../data/li6800.v1.Rdata"

#'### Function for reading in each of the Li-6800 files
readli6800files <- function(file)
{
  tbl <- suppressMessages(readxl::read_xlsx(paste0("../data/", file),
                                            sheet = "Sorted Data",
                                            col_names = FALSE)[,c(1:6, 11, 16, 21, (26*4+7))]) ##column DG = 26*4+7)
  
  
  ##tidy up
  colnames(tbl) <- c("obs", "time", "elapsed","date","hhmmss","time.elapsed", "A", 
                     "gsw", "TleafEB", "light")
  tbl <- tbl %>% slice(-c(1:3))
  
  ##convert variables
  tbl <- tbl %>% mutate(across(c("A", "gsw", "TleafEB", "light"), as.numeric)) %>% as.data.frame()
  tbl$Timestamp <- as.POSIXct(tbl$date, format = "%Y%m%d %H:%M:%S")
  
  ##calculate time difference between each observation
  tbl$time.diff <- c(0, diff(tbl$Timestamp))
  tbl$xElapsed <- cumsum(tbl$time.diff)
  tbl$xElapsed.min <- tbl$xElapsed/60
  
  ##append treatment information using filename
  trt.info <- strsplit(file, "_")
  tbl$Line <- trt.info[[1]][1]
  tbl$Temperature <- trt.info[[1]][2]
  tbl$Watering <- trt.info[[1]][3]
  tbl$Rep <- trt.info[[1]][4]
  
  return(tbl)
}

li6800.dat <- do.call(rbind, lapply(inputFiles, readli6800files))
head(li6800.dat)
nrow(li6800.dat) #4392

li6800.dat$Rep <- factor(li6800.dat$Rep)
li6800.dat$Rep <- dae::fac.recast(li6800.dat$Rep,
                                  newlevels = c(1:3))

##put together a unique identifier for each dataset (day)
li6800.dat$TWR <- paste(li6800.dat$Temperature, li6800.dat$Watering, li6800.dat$Rep, sep ="_")

#'### How many measurements on each day (for each TWR)?
table(li6800.dat$TWR)

#'### Check if the lights were off for the first measurement? If so, set A and gs to NA
li6800.dat %>% filter(xElapsed < 60 & light > 1) %>% select("TWR", "time.diff", "light")
##this confirms that the light was only on for the first measurement of some treatments

li6800.dat[(li6800.dat$xElapsed == 0 & li6800.dat$light > 1 & !is.na(li6800.dat$A)), 
           c("A", "gsw")] <- NA

#'### Add location of the plant that was sampled
#' This information orginally came from the Porometer data spreadsheet, but I have created my own little spreadsheet to import here
plant.sampled.df <- read.csv("../data/PlantPositionsLookupTable.csv")
##remove those plants that have data used for analysis as "No", as this was not supplied to us
plant.sampled.df <- plant.sampled.df[plant.sampled.df$Data.Used.for.Analyses. == "Yes",]

plant.sampled.df <- within(plant.sampled.df,{
  Temperature <- Room
  Temperature <- factor(Temperature, levels = c("Ambient", "Heat"),
                        labels = c("Amb", "Heat"))
  Watering <- factor(Watering, levels = c("Limited", "Unlimited"),
                        labels = c("Lim", "Unlim"))
  Rep <- factor(Rep)
})

li6800.dat <- left_join(li6800.dat, plant.sampled.df[, c("Rep", "Watering", "Temperature", "Location")])
nrow(li6800.dat)

#'## Define the light and dark regions
hist(li6800.dat$light)
li6800.dat$Light <- "Off"
li6800.dat$Light[li6800.dat$light > 1] <- "On"
table(li6800.dat$TWR, li6800.dat$Light)

#'## Find anomalies in gsw
#' set the value to NA
li6800.dat[li6800.dat$gsw < -0.01 & !is.na(li6800.dat$gsw),]
li6800.dat$gsw[li6800.dat$gsw < -0.01] <- NA

#'## Visualise the A and gsw data
plotProfiles(li6800.dat, response = "A",
             individuals = "Location", times = "xElapsed.min",
             facet.x = "Temperature", facet.y = "Watering",  colour.column = "TWR",
             y.title = "A",
             x.title = "Seconds",
             title = "A",
             breaks.spacing.x = 30,
             alpha = 0.75)

plotProfiles(li6800.dat, response = "gsw",
             individuals = "Location", times = "xElapsed.min",
             facet.x = "Temperature", facet.y = "Watering",  colour.column = "TWR",
             y.title = "gsw",
             x.title = "Seconds",
             title = "gsw",
             breaks.spacing.x = 30,
             alpha = 0.75)

save.image(writeFile)
