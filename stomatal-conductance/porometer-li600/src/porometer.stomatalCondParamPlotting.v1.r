#'# Wheat stomatal kinetics (heat & drought)
#'## --- Ed Chaplin (USyd), *Apr--June 2025* ---
#'### Li-600 STOMATAL CONDUCTANCE MODEL PLOTTING PARAMETERS
#'
#' * Plotting the resulting parameter estimates
#' * 
#'
#'### Load the libraries
library(knitr)
#knitr::spin("porometer.stomatalCondParamPlotting.v1.r")
library(ggplot2)
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, 
               "fig.width" = 7.5, "fig.height" = 10, "comment"=NA, fig.path = "figures/")
library(dae)
library(dplyr)
library(brms) ##package used to fit bayesian regression model using Stan
library(tidyr)
library(emmeans)
library(tidybayes)
library(pracma)
library(patchwork)

#'### Set up constants and files
options(width=90)
load("../data/stomatalCondParamEstimation_porometer_i.v1.Rdata")
load("../data/stomatalCondParamEstimation_porometer_d.v1.Rdata")
saveFile <- "../data/stomatalCondParamPlotting_porometer.v1.Rdata"

#'### Bring together parameters from across different analyses
head(li600.i.param)
head(li600.d.param)

##bring together the parameters into one dataframe, but first make sure names are unique
li600.i.param <- li600.i.param %>% rename_with(~paste0("i_", .x),
                                         -all_of(c("Room_Temperature",
                                                   "Device_Watering",
                                                   "Line")))
li600.d.param <- li600.d.param %>% rename_with(~paste0("d_", .x),
                                                 -all_of(c("Room_Temperature",
                                                           "Device_Watering",
                                                           "Line")))

li600.all.params <- left_join(li600.i.param,
                              li600.d.param)
li600.all.params

#'### Put together plots of each parameter separately to begin with and then we'll stitch them together

##increasing params
param.i <- li600.all.params %>% select(Room_Temperature, Device_Watering, Line,
                                       grep("i_", names(.))) %>%
  select(-c(grep("tabquant", names(.)))) %>% 
  pivot_longer(col = 4:ncol(.), values_to = "vals", names_to = "param") %>%
  separate(param, into = c("i", "param", "quantity"), sep="_") %>%
  pivot_wider(id_cols = c("Room_Temperature", "Device_Watering", "Line", "i", "param"),
              names_from = "quantity", values_from = "vals")
param.i <- within(param.i,{
  paramLabel <- factor(param, levels = c("g0", "g1", "k", "lambda", "Slmax"),
                  labels = c(expression(italic(r)[0[i]]),
                             expression(italic(g)[smax]),
                             expression(italic(k)[i]),
                             expression(lambda[i]),
                             expression(italic(Sl)[max[i]])))
})

g0.g1.plot <- ggplot(param.i %>% filter(param %in% c("g0", "g1")), 
                     aes(x = Line, y = emmean, 
                         group = Line)) +
  geom_point(aes(colour = paramLabel), 
             stat = "identity",
             size = 2,
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower.HPD, 
                    ymax = upper.HPD,
                    colour = paramLabel),
                width = 0.2,
                position = position_dodge(width = 0.9)) +
  facet_grid(Device_Watering ~ Room_Temperature) +
  theme_bw() +
  #scale_x_discrete(labels = scales::label_parse()) +
  scale_colour_discrete(labels = scales::label_parse(),
                        name = "Parameter")+
  xlab("Genotype") + 
  ylab(expression(Parameter~estimate~(mol~m^-2~s^-1))) +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("(a)")
g0.g1.plot

ki.plot <- ggplot(param.i %>% filter(param %in% c("k")), 
                     aes(x = Line, y = emmean, 
                         group = Line)) +
  geom_point(colour = "grey50", 
             stat = "identity",
             size = 2,
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower.HPD, 
                    ymax = upper.HPD),
                width = 0.2,
                position = position_dodge(width = 0.9),
                colour = "grey50",) +
  facet_grid(Device_Watering ~ Room_Temperature) +
  theme_bw() +
  #scale_x_discrete(labels = scales::label_parse()) +
  scale_colour_discrete(labels = scales::label_parse(),
                        name = "Parameter")+
  xlab("Genotype") + 
  ylab(expression(italic(k)[i]~(min))) +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("(b)")
ki.plot

Slmax.i.plot <- ggplot(param.i %>% filter(param %in% c("Slmax")), 
                  aes(x = Line, y = emmean, 
                      group = Line)) +
  geom_point(colour = "grey50", 
             stat = "identity",
             size = 2,
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower.HPD, 
                    ymax = upper.HPD),
                width = 0.2,
                position = position_dodge(width = 0.9),
                colour = "grey50",) +
  facet_grid(Device_Watering ~ Room_Temperature) +
  theme_bw() +
  #scale_x_discrete(labels = scales::label_parse()) +
  scale_colour_discrete(labels = scales::label_parse(),
                        name = "Parameter")+
  xlab("Genotype") + 
  ylab(expression(italic(Sl)[max[i]]~(mmol~m^-2~s^-2))) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, hjust=1, vjust =0.5)) +
  ggtitle("(c)")
Slmax.i.plot

(g0.g1.plot) / (ki.plot) / (Slmax.i.plot)
ggsave("figures/li600_i_params.png", height = 16, width = 12)
ggsave("figures/li600_i_params.pdf", height = 16, width = 12)
ggsave("figures/li600_i_params.tiff", height = 16, width = 12)

##decreasing params
param.d <- li600.all.params %>% select(Room_Temperature, Device_Watering, Line,
                                       grep("d_", names(.))) %>%
  select(-c(grep("tabquant", names(.)))) %>% 
  pivot_longer(col = 4:ncol(.), values_to = "vals", names_to = "param") %>%
  separate(param, into = c("d", "param", "quantity"), sep="_") %>%
  pivot_wider(id_cols = c("Room_Temperature", "Device_Watering", "Line", "d", "param"),
              names_from = "quantity", values_from = "vals")
param.d <- within(param.d,{
  paramLabel <- factor(param, levels = c("g0", "g1", "k", "lambda", "Slmax"),
                       labels = c(expression(italic(r)[0[d]]),
                                  expression(italic(g)[smin]),
                                  expression(italic(k)[d]),
                                  expression(lambda[d]),
                                  expression(italic(Sl)[max[d]])))
})

g0.g1.d.plot <- ggplot(param.d %>% filter(param %in% c("g0", "g1")), 
                     aes(x = Line, y = emmean, 
                         group = Line)) +
  geom_point(aes(colour = paramLabel), 
             stat = "identity",
             size = 2,
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower.HPD, 
                    ymax = upper.HPD,
                    colour = paramLabel),
                width = 0.2,
                position = position_dodge(width = 0.9)) +
  facet_grid(Device_Watering ~ Room_Temperature) +
  theme_bw() +
  #scale_x_discrete(labels = scales::label_parse()) +
  scale_colour_discrete(labels = scales::label_parse(),
                        name = "Parameter")+
  xlab("Genotype") + 
  ylab(expression(Parameter~estimate~(mol~m^-2~s^-1))) +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("(a)")
g0.g1.d.plot

kd.plot <- ggplot(param.d %>% filter(param %in% c("k")), 
                  aes(x = Line, y = emmean, 
                      group = Line)) +
  geom_point(colour = "grey50", 
             stat = "identity",
             size = 2,
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower.HPD, 
                    ymax = upper.HPD),
                width = 0.2,
                position = position_dodge(width = 0.9),
                colour = "grey50",) +
  facet_grid(Device_Watering ~ Room_Temperature) +
  theme_bw() +
  #scale_x_discrete(labels = scales::label_parse()) +
  scale_colour_discrete(labels = scales::label_parse(),
                        name = "Parameter")+
  xlab("Genotype") + 
  ylab(expression(italic(k)[d]~(min))) +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("(b)")
kd.plot

Slmax.d.plot <- ggplot(param.d %>% filter(param %in% c("Slmax")), 
                       aes(x = Line, y = emmean, 
                           group = Line)) +
  geom_point(colour = "grey50", 
             stat = "identity",
             size = 2,
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower.HPD, 
                    ymax = upper.HPD),
                width = 0.2,
                position = position_dodge(width = 0.9),
                colour = "grey50",) +
  facet_grid(Device_Watering ~ Room_Temperature) +
  theme_bw() +
  #scale_x_discrete(labels = scales::label_parse()) +
  scale_colour_discrete(labels = scales::label_parse(),
                        name = "Parameter")+
  xlab("Genotype") + 
  ylab(expression(italic(Sl)[max[d]]~(mmol~m^-2~s^-2))) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, hjust=1, vjust =0.5)) +
  ggtitle("(c)")
Slmax.d.plot

(g0.g1.d.plot) / (kd.plot) / (Slmax.d.plot)
ggsave("figures/li600_d_params.png", height = 16, width = 12)
ggsave("figures/li600_d_params.pdf", height = 16, width = 12)
ggsave("figures/li600_d_params.tiff", height = 16, width = 12)

#'### Put together a plot of the fits from the model for each trt combo
li600.on.gs.plot <- li600.on.plot
li600.on.gs.plot$fitted.fixed <- fitted(fit_on_final, re_formula = NA)[,1]

ggplot(li600.on.gs.plot, aes(time.adj, gsw)) + 
  geom_point(alpha = 0.5, colour = "gray50") +
  geom_line(aes(y=fitted.fixed), colour = "black", linewidth = 1)+
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  theme_bw() +
  xlab("Time (min)")+
  ylab(expression(italic(g)[s]~(mol~m^-2~s^-1))) +
  theme(strip.text.x = element_text(size = 8))
ggsave("figures/li600_gs_on_FittedResponse.png", width=20, height = 12)
ggsave("figures/li600_gs_on_FittedResponse.pdf", width=20, height = 12)
ggsave("figures/li600_gs_on_FittedResponse.tiff", width=20, height = 12)

li600.off.gs.plot <- li600.off.plot
li600.off.gs.plot$fitted.fixed <- fitted(fit_off_final, re_formula = NA)[,1]

ggplot(li600.off.gs.plot, aes(time.adj, gsw)) + 
  geom_point(alpha = 0.5, colour = "gray50") +
  geom_line(aes(y=fitted.fixed), colour = "black", linewidth = 1)+
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  theme_bw() +
  xlab("Time (min)")+
  ylab(expression(italic(g)[s]~(mol~m^-2~s^-1))) +
  theme(strip.text.x = element_text(size = 8))
ggsave("figures/li600_gs_off_FittedResponse.png", width=20, height = 12)
ggsave("figures/li600_gs_off_FittedResponse.pdf", width=20, height = 12)
ggsave("figures/li600_gs_off_FittedResponse.tiff", width=20, height = 12)

#'### Put together a combined plot
li600.gs.all.plot <- bind_rows(li600.on.gs.plot, li600.off.gs.plot)
li600.gs.all.plot$Light <- factor(li600.gs.all.plot$Light,
                                   levels = c("On", "Off"))
li600.gs.all.plot$new.time <- li600.gs.all.plot$time.adj + 
  min(li600.gs.all.plot$time.elapsed)
li600.gs.all.plot$new.time[li600.gs.all.plot$Light == "Off"] <- 
  li600.gs.all.plot$new.time[li600.gs.all.plot$Light == "Off"] +
  (min(li600.gs.all.plot$time.elapsed[li600.gs.all.plot$Light == "Off"]) -
     min(li600.gs.all.plot$time.elapsed))

ggplot(li600.gs.all.plot, aes(new.time, gsw)) + 
  geom_point(alpha = 0.5, colour = "gray50", shape = 17) +
  geom_line(aes(y=fitted.fixed, group = Light, colour = Light), 
          linewidth = 1) +
  geom_vline(aes(xintercept = min(time.elapsed[Light == "Off"]))) +
  scale_colour_manual(values = c("orange", "black")) +
  facet_grid(Room_Temperature + Device_Watering ~ Line) +
  theme_bw() +
  xlab("Time (min)")+
  ylab(expression(italic(g)[s]~(mol~m^-2~s^-1))) +
  theme(legend.position = "none",
        strip.text.x = element_text(size = 8))
ggsave("figures/li600.gs.fitted.response.png", width=20, height = 12)
ggsave("figures/li600.gs.fitted.response.pdf", width=20, height = 12)
ggsave("figures/li600.gs.fitted.response.tiff", width=20, height = 12)


#'### Save data and export
save.image(saveFile)

