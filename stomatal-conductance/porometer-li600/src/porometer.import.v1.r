#'# Wheat stomatal kinetics (heat & drought)
#'## --- Ed Chaplin (USyd), *Apr--June 2025* ---
#'### Porometer DATA
#'
#' * Importing porometer data obtained from handheld Li-600 device
#' * There are individual tabs for each treatment sampled with the Li-600, accompanied by a All Data Compiled tab
#' * All reps were used for Li-600 measurements
#' * The key data parameter to assess is 'gsw' (which we can refer to as 'gs') (Column AS)
#' * The light level (off and on) is denoted in column AQ
#' * Measurements in each day were taken by two Operators (Operators confounded with Watering trts)
#' * Lines within Chambers and Watering trts were measured repeatedly in a systematic way by doing laps on the pots
#' * I've defined a time.elapsed variable to measure the overall time taken by each Operator during the sampling window
#' * The times that the lights were on and off for each Operator do vary slightly
#'
#'### Load the libraries
library(knitr)
#knitr::spin("porometer.import.v1.r")
library(ggplot2)
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, 
               "fig.width" = 7.5, "fig.height" = 10, "comment"=NA, fig.path = "figures/")
library(dplyr)

#'### Set up constants and files
options(width=90)
inputFile <- "../data/Porometer Data - 20250701.xlsx"
writeFile  <- "../data/porometer.v1.Rdata"

#'### Read in the porometer data from the All Data Compiled tab
li600.dat <- suppressMessages(readxl::read_xlsx(inputFile,
                         sheet = "All Data Compiled",
                         col_names = FALSE)[,c(1:28, 35:38, 43:45)])

##tidy up
colnames(li600.dat) <- li600.dat[1,]
colnames(li600.dat)[c(1:2, 29:35)] <- c("Sample.ID", "Barcode.ID",
                               "ObsNum", "Time", "MinsFromStart","Date_Sys",
                               "Light", "Lap", "gsw")
li600.dat <- li600.dat %>% slice(-c(1:3))

##convert variables to appropriate format
li600.dat <- li600.dat %>% mutate(across(c("gsw", "Order", "Lap", "MinsFromStart", "ObsNum",
                                           "CEBPosition", "CEBBlock","CEBSubUnit", "CEBCell",
                                           "CEBPCoord", "Replicate", "Entry"),
                                         as.numeric)) %>% as.data.frame()

##Date variable is incorrect - use Data_Sys instead, but this also contains errors
##Email from Ed confirms that changes need to be made to the Date_Sys column
##From Ed: 1.	Heat_Limited_1 should be 2025-06-13 NOT 2025-06-12
##         2.	Ambient_Limited_4 should be 2025-06-16 NOT 2025-06-15 
li600.dat$Date_Sys[li600.dat$Temperature == "Heat" &
                     li600.dat$Watering == "Limited" &
                     li600.dat$Replicate == 1] <- "2025-06-13"
li600.dat$Date_Sys[li600.dat$Temperature == "Ambient" &
                     li600.dat$Watering == "Limited" &
                     li600.dat$Replicate == 4] <- "2025-06-16"

li600.dat$Date_Sys<- lubridate::ymd(li600.dat$Date_Sys)
li600.dat$Timestamp <- paste(li600.dat$Date_Sys, li600.dat$Time)
li600.dat$Timestamp <- lubridate::ymd_hms(li600.dat$Timestamp)
str(li600.dat)
## Correct Operator variable
li600.dat$Operator[li600.dat$Operator == "Fiona Nourish"] <- "Fiona"

##calculate time difference between each observation
li600.dat <- li600.dat[order(li600.dat$Date_Sys, li600.dat$Operator, li600.dat$Timestamp),]
li600.dat <- li600.dat %>% group_by(Date_Sys, Operator) %>% 
  mutate(time.diff = as.numeric(Timestamp - lag(Timestamp, 
                                                default = first(Timestamp)))/60,
         time.elapsed = cumsum(time.diff)) %>% as.data.frame()

head(li600.dat)
nrow(li600.dat) #5003

#'### Understand the structure in the measurments
#' * Ed says, One complete replicate from one growth chamber (heat or ambient) was sampled per day
#' * Within each phase, stomatal conductance was recorded sequentially by moving from plant to plant, continuously within the phase, collecting approximately 10-15 measurements per plant per phase.
#' * To increase the temporal resolution of measurements, two operators with two porometers sampled plants concurrently.
#' * One operator measured the well watered treatment (10 plants) while the other measured the drought treatment (10 plants) within the same replicate. 
table(li600.dat$Date_Sys)
table(li600.dat$Treatment, li600.dat$Date_Sys)
table(li600.dat$Treatment, li600.dat$Operator)
##this demonstrates confounding between Operator and Treatment - need to keep this in mind for models
table(li600.dat$Date_Sys, li600.dat$CEBRoom)
table(li600.dat$Date_Sys, li600.dat$Temperature)
table(li600.dat$Date_Sys, li600.dat$Watering)
table(li600.dat$Date_Sys, li600.dat$Line)
table(li600.dat$Date_Sys, li600.dat$Replicate)
##confounding between Dates and Room_Temperature, along with Replicate

##Define some new variables
li600.dat$Device <- 1
li600.dat$Device[li600.dat$Operator != "Ed"] <- 2

li600.dat <- within(li600.dat, {
  Room_Temperature <- paste(CEBRoom, Temperature, sep="_")
  Operator_Watering <- paste(Operator, Watering, sep="_")
  Device_Watering <- paste(paste0("D", Device), Watering, sep="_")
})

#'### Format variables
li600.dat <- within(li600.dat, {
  CEBRoom <- factor(CEBRoom)
  CEBLane <- factor(CEBLane)
  CEBPosition <- factor(CEBPosition)
  CEBSide <- factor(CEBSide)
  CEBLocation <- factor(CEBLocation)
  CEBBlock <- factor(CEBBlock)
  CEBSubUnit <- factor(CEBSubUnit)
  CEBCell <- factor(CEBCell)
  Treatment <- factor(Treatment)
  Replicate <- factor(Replicate)
  Room_Temperature <- factor(Room_Temperature)
  Operator_Watering <- factor(Operator_Watering)
  Device_Watering <- factor(Device_Watering)
  Line <- factor(Line)
})


#'## Visualise the gsw data
#'In this case, we'll see that some Lines will be exposed to light on/off conditions for much longer than other Lines before being measured
##building manually, as want to present points
##when did the lights go off and on?
suppressMessages(li600.dat %>% filter(Light == "On") %>%
  group_by(Date_Sys, Operator) %>%
  summarise(min = min(time.elapsed)))

suppressMessages(li600.dat %>% filter(Light == "Off" &
                       time.elapsed > 45.6) %>%
  group_by(Date_Sys, Operator) %>%
  summarise(min = min(time.elapsed)))
##Inconsistent across operators???

ggplot(li600.dat %>% filter(Date_Sys == "2025-06-12"),
       aes(x = time.elapsed, y= gsw, group = Sample.ID)) +
  geom_point(aes(colour = Line)) +
  geom_point(data = li600.dat %>% filter(Date_Sys == "2025-06-12" & Light == "On"),
             shape = 1, colour = "black") +
  geom_line(aes(colour = Line)) +
  xlab("Minutes") + ylab("gsw") +
  ggtitle("gsw - CE8_Ambient, Rep 2 - 2025-06-12") +
  facet_grid(Operator_Watering ~ Line) +
  theme_bw() +
  theme(legend.position = "top")

ggplot(li600.dat %>% filter(Date_Sys == "2025-06-12" & Light == "On"),
       aes(x = time.elapsed, y= gsw, group = Sample.ID)) +
  geom_point(aes(colour = Line)) +
  geom_line(aes(colour = Line)) +
  xlab("Minutes") + ylab("gsw") +
  ggtitle("gsw - Light ON - CE8_Ambient, Rep 2 - 2025-06-12") +
  facet_grid(Operator_Watering ~ Line) +
  theme_bw() +
  theme(legend.position = "top")

ggplot(li600.dat %>% filter(Temperature == "Ambient"),
       aes(x = time.elapsed, y= gsw, group = Sample.ID)) +
  geom_point(aes(colour = Replicate)) +
  geom_point(data = li600.dat %>% filter(Temperature == "Ambient" & Light == "On"),
             shape = 1, colour = "black") +
  geom_line(aes(colour = Replicate)) +
  xlab("Minutes") + ylab("gsw") +
  ggtitle("gsw - CE8_Ambient") +
  facet_grid(Device_Watering ~ Line) +
  theme_bw() +
  theme(legend.position = "top")

ggplot(li600.dat %>% filter(Temperature == "Ambient" &
                              Light == "On"),
       aes(x = time.elapsed, y= gsw, group = Sample.ID)) +
  geom_point(aes(colour = Replicate)) +
  geom_line(aes(colour = Replicate)) +
  xlab("Minutes") + ylab("gsw") +
  ggtitle("gsw - CE8_Ambient - Light ON") +
  facet_grid(Device_Watering ~ Line) +
  theme_bw() +
  theme(legend.position = "top")

ggplot(li600.dat %>% filter(Temperature == "Heat"),
       aes(x = time.elapsed, y= gsw, group = Sample.ID)) +
  geom_point(aes(colour = Replicate)) +
  geom_point(data = li600.dat %>% filter(Temperature == "Heat" & Light == "On"),
             shape = 1, colour = "black") +
  geom_line(aes(colour = Replicate)) +
  xlab("Minutes") + ylab("gsw") +
  ggtitle("gsw - CE7_Heat") +
  facet_grid(Device_Watering ~ Line) +
  theme_bw() +
  theme(legend.position = "top")

ggplot(li600.dat %>% filter(Temperature == "Heat" &
                              Light == "On"),
       aes(x = time.elapsed, y= gsw, group = Sample.ID)) +
  geom_point(aes(colour = Replicate)) +
  geom_line(aes(colour = Replicate)) +
  xlab("Minutes") + ylab("gsw") +
  ggtitle("gsw - CE7_Heat - Light ON") +
  facet_grid(Device_Watering ~ Line) +
  theme_bw() +
  theme(legend.position = "top")

save.image(writeFile)
