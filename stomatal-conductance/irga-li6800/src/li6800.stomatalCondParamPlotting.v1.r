#'# Wheat stomatal kinetics (heat & drought)
#'## --- Ed Chaplin (USyd), *Apr--June 2025* ---
#'### Li-6800 STOMATAL CONDUCTANCE MODEL PLOTTING PARAMETERS
#'
#' * Plotting the resulting parameter estimates
#' * 
#'
#'### Load the libraries
library(knitr)
#knitr::spin("li6800.stomatalCondParamPlotting.v1.r")
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
load("../data/stomatalCondParamEstimation_i.v1.Rdata")
load("../data/stomatalCondParamEstimation_d.v1.Rdata")
load("../data/stomatalCondParamEstimation_A_i.v1.Rdata")
saveFile <- "../data/stomatalCondParamPlotting.v1.Rdata"

#'### Bring together parameters from across different analyses
head(li6800.A.i.param)
head(li6800.i.param)
head(li6800.d.param)

##bring together the parameters into one dataframe, but first make sure names are unique
li6800.i.param <- li6800.i.param %>% rename_with(~paste0("i_", .x),
                                         -all_of(c("Room_Temperature",
                                                   "Watering")))
li6800.d.param <- li6800.d.param %>% rename_with(~paste0("d_", .x),
                                                 -all_of(c("Room_Temperature",
                                                           "Watering")))



li6800.all.params <- left_join(li6800.A.i.param,
                               li6800.i.param) %>% left_join(li6800.d.param)
li6800.all.params

##recode the factors
li6800.all.params <- within(li6800.all.params, {
  Room_Temperature <- fac.recode(Room_Temperature,
                                 newlevels = c("CER7_Heat",
                                               "CER8_Ambient"))
  Watering <- fac.recode(Watering,
                         newlevels = c("Limited",
                                       "Unlimited"))
  
})

#'### Put together plots of each parameter separately to begin with and then we'll stitch them together

##increasing params
param.i <- li6800.all.params %>% select(Room_Temperature, Watering, grep("i_", 
                                                                         names(.))) %>%
  # bind_cols(li6800.all.params %>% select(grep("A_", names(.)),
  #                                             grep("t95", names(.)),
  #                                             grep("A95", names(.)),
  #                                             grep("Limit", names(.)))) %>%
  select(-c(grep("tabquant", names(.)))) %>% 
  pivot_longer(col = 3:ncol(.), values_to = "vals", names_to = "param") %>%
  separate(param, into = c("i", "param", "quantity"), sep="_") %>%
  pivot_wider(id_cols = c("Room_Temperature", "Watering", "i", "param"),
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
                     aes(x = Watering, y = emmean, 
                         group = Watering)) +
  geom_point(aes(colour = paramLabel), 
             stat = "identity",
             size = 2,
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower.HPD, 
                    ymax = upper.HPD,
                    colour = paramLabel),
                width = 0.1,
                position = position_dodge(width = 0.9)) +
  facet_grid( ~ Room_Temperature) +
  theme_bw() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_colour_discrete(labels = scales::label_parse(),
                        name = "Parameter")+
  xlab("Watering treatment") + 
  ylab(expression(Parameter~estimate~(mol~m^-2~s^-1))) +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("(a)")
g0.g1.plot

ki.plot <- ggplot(param.i %>% filter(param %in% c("k")), 
                     aes(x = Watering, y = emmean, 
                         group = Watering)) +
  geom_point(colour = "grey50", 
             stat = "identity",
             size = 2,
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower.HPD, 
                    ymax = upper.HPD),
                width = 0.1,
                position = position_dodge(width = 0.9),
                colour = "grey50") +
  facet_grid( ~ Room_Temperature) +
  theme_bw() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_colour_discrete(labels = scales::label_parse(),
                        name = "Parameter")+
  xlab("Watering treatment") + 
  ylab(expression(italic(k)[i]~(min))) +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("(b)")
ki.plot

Slmax.i.plot <- ggplot(param.i %>% filter(param %in% c("Slmax")), 
                  aes(x = Watering, y = emmean, 
                      group = Watering)) +
  geom_point(colour = "grey50", 
             stat = "identity",
             size = 2,
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower.HPD, 
                    ymax = upper.HPD),
                width = 0.1,
                position = position_dodge(width = 0.9),
                colour = "grey50") +
  facet_grid( ~ Room_Temperature) +
  theme_bw() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_colour_discrete(labels = scales::label_parse(),
                        name = "Parameter")+
  xlab("Watering treatment") + 
  ylab(expression(italic(Sl)[max[i]]~(mmol~m^-2~s^-2))) +
  theme(legend.position = "top") +
  ggtitle("(c)")
Slmax.i.plot

##increasing params A
param.i.A <- li6800.all.params %>% select(Room_Temperature, Watering, grep("A_", 
                                                                         names(.))) %>%
  bind_cols(li6800.all.params %>% select(grep("t95", names(.)),
                                              grep("A95", names(.)),
                                              grep("Limit", names(.)))) %>%
  select(-c(grep("tabquant", names(.)))) %>% 
  pivot_longer(col = 3:ncol(.), values_to = "vals", names_to = "param") %>%
  separate(param, into = c("param", "A", "quantity"), sep="_")
param.i.A$quantity[param.i.A$param %in% c("t95", 
                                          "A95", 
                                          "Limit")] <- param.i.A$A[param.i.A$param %in% c("t95", "A95", "Limit")]
param.i.A <- param.i.A %>% select(-A) %>%
  pivot_wider(id_cols = c("Room_Temperature", "Watering", "param"),
              names_from = "quantity", values_from = "vals")
param.i.A <- within(param.i.A,{
  paramLabel <- factor(param, levels = c("g0", "A95", "g1", "k","t95", "Limit"), 
                       labels = c(expression(italic(A)[0]),
                                  expression(italic(A)[95]),
                                  expression(italic(A)[max]),
                                  expression(italic(k)[italic(A)]),
                                  expression(italic("T")[95]*italic(A)),
                                  expression(italic(g)[slimit]*italic(A))))
})

g0.g1.A.plot <- ggplot(param.i.A %>% filter(param %in% c("g0", "g1", "A95")), 
                     aes(x = Watering, y = emmean, 
                         group = Watering)) +
  geom_point(aes(colour = paramLabel), 
             stat = "identity",
             size = 2,
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower.HPD, 
                    ymax = upper.HPD,
                    colour = paramLabel),
                width = 0.1,
                position = position_dodge(width = 0.9)) +
  facet_grid( ~ Room_Temperature) +
  theme_bw() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_colour_discrete(labels = scales::label_parse(),
                        name = "Parameter")+
  xlab("Watering treatment") + 
  ylab(expression(Parameter~estimate~(mu*mol~m^-2~s^-1))) +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("(d)")
g0.g1.A.plot

t95.A.plot <- ggplot(param.i.A %>% filter(param %in% c("t95")), 
                   aes(x = Watering, y = emmean, 
                       group = Watering)) +
  geom_point(colour = "grey50",
             stat = "identity",
             size = 2,
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower.HPD, 
                    ymax = upper.HPD),
                width = 0.1,
                position = position_dodge(width = 0.9),
                colour = "grey50") +
  facet_grid( ~ Room_Temperature) +
  theme_bw() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_colour_discrete(labels = scales::label_parse(),
                        name = "Parameter")+
  xlab("Watering treatment") + 
  ylab(expression(italic("T")[95]*italic(A)~(min))) +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("(e)")
t95.A.plot

limit.A.plot <- ggplot(param.i.A %>% filter(param %in% c("Limit")), 
                     aes(x = Watering, y = emmean*100, 
                         group = Watering)) +
  geom_point(colour = "grey50",
             stat = "identity",
             size = 2,
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower.HPD*100, 
                    ymax = upper.HPD*100),
                width = 0.1,
                position = position_dodge(width = 0.9),
                colour = "grey50") +
  facet_grid( ~ Room_Temperature) +
  theme_bw() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_colour_discrete(labels = scales::label_parse(),
                        name = "Parameter")+
  xlab("Watering treatment") + 
  ylab(expression(italic(g)[slimit]*italic(A)~("%"))) +
  theme(legend.position = "top") +
  ggtitle("(f)")
limit.A.plot

(g0.g1.plot + g0.g1.A.plot) / (ki.plot + t95.A.plot) / (Slmax.i.plot + limit.A.plot)
ggsave("figures/li6800_i_params.png", height = 12, width = 12)
ggsave("figures/li6800_i_params.pdf", height = 12, width = 12)
ggsave("figures/li6800_i_params.tiff", height = 12, width = 12)

##decreasing params
param.d <- li6800.all.params %>% select(Room_Temperature, Watering, grep("d_", 
                                                                         names(.))) %>%
  # bind_cols(li6800.all.params %>% select(grep("A_", names(.)),
  #                                             grep("t95", names(.)),
  #                                             grep("A95", names(.)),
  #                                             grep("Limit", names(.)))) %>%
  select(-c(grep("tabquant", names(.)))) %>% 
  pivot_longer(col = 3:ncol(.), values_to = "vals", names_to = "param") %>%
  separate(param, into = c("d", "param", "quantity"), sep="_") %>%
  pivot_wider(id_cols = c("Room_Temperature", "Watering", "d", "param"),
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
                     aes(x = Watering, y = emmean, 
                         group = Watering)) +
  geom_point(aes(colour = paramLabel), 
             stat = "identity",
             size = 2,
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower.HPD, 
                    ymax = upper.HPD,
                    colour = paramLabel),
                width = 0.1,
                position = position_dodge(width = 0.9)) +
  facet_grid( ~ Room_Temperature) +
  theme_bw() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_colour_discrete(labels = scales::label_parse(),
                        name = "Parameter")+
  xlab("Watering treatment") + 
  ylab(expression(Parameter~estimate~(mol~m^-2~s^-1))) +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("(a)")
g0.g1.d.plot

kd.plot <- ggplot(param.d %>% filter(param %in% c("k")), 
                  aes(x = Watering, y = emmean, 
                      group = Watering)) +
  geom_point(colour = "grey50", 
             stat = "identity",
             size = 2,
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower.HPD, 
                    ymax = upper.HPD),
                width = 0.1,
                position = position_dodge(width = 0.9),
                colour = "grey50") +
  facet_grid( ~ Room_Temperature) +
  theme_bw() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_colour_discrete(labels = scales::label_parse(),
                        name = "Parameter")+
  xlab("Watering treatment") + 
  ylab(expression(italic(k)[d]~(min))) +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("(b)")
kd.plot

Slmax.d.plot <- ggplot(param.d %>% filter(param %in% c("Slmax")), 
                       aes(x = Watering, y = emmean, 
                           group = Watering)) +
  geom_point(colour = "grey50", 
             stat = "identity",
             size = 2,
             position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower.HPD, 
                    ymax = upper.HPD),
                width = 0.1,
                position = position_dodge(width = 0.9),
                colour = "grey50") +
  facet_grid( ~ Room_Temperature) +
  theme_bw() +
  scale_x_discrete(labels = scales::label_parse()) +
  scale_colour_discrete(labels = scales::label_parse(),
                        name = "Parameter")+
  xlab("Watering treatment") + 
  ylab(expression(italic(Sl)[max[d]]~(mmol~m^-2~s^-2))) +
  theme(legend.position = "top") +
  ggtitle("(c)")
Slmax.d.plot

(g0.g1.d.plot) / (kd.plot) / (Slmax.d.plot)
ggsave("figures/li6800_d_params.png", height = 10, width = 8)
ggsave("figures/li6800_d_params.pdf", height = 10, width = 8)
ggsave("figures/li6800_d_params.tiff", height = 10, width = 8)

#'### Put together a plot of the fits from the model for each trt combo
li6800.on.A.plot <- li6800.on.plot
li6800.on.A.plot <- within(li6800.on.A.plot, {
  Room_Temperature <- fac.recode(Room_Temperature,
                                 newlevels = c("CER7_Heat",
                                               "CER8_Ambient"))
  Watering <- fac.recode(Watering,
                         newlevels = c("Limited",
                                       "Unlimited"))
  
})

ggplot(li6800.on.A.plot, aes(tmin, A)) + 
  geom_point(alpha = 0.5, colour = "gray50") +
  geom_line(aes(y=fitted.fixed), colour = "black", linewidth = 1)+
  facet_grid(Room_Temperature ~ Watering) +
  theme_bw() +
  xlab("Time (min)")+
  ylab(expression(italic(A)~(mu*mol~m^-2~s^-1)))
ggsave("figures/Li6800_A_FittedResponse.png", width=8, height = 8)
ggsave("figures/Li6800_A_FittedResponse.pdf", width=8, height = 8)
ggsave("figures/Li6800_A_FittedResponse.tiff", width=8, height = 8)

li6800.on.gs.plot <- li6800.on.plot
li6800.on.gs.plot <- within(li6800.on.gs.plot, {
  Room_Temperature <- fac.recode(Room_Temperature,
                                 newlevels = c("CER7_Heat",
                                               "CER8_Ambient"))
  Watering <- fac.recode(Watering,
                         newlevels = c("Limited",
                                       "Unlimited"))
  
})
li6800.on.gs.plot$fitted.fixed <- fitted(fit_on_final, re_formula = NA)[,1]

ggplot(li6800.on.gs.plot, aes(tmin, gsw)) + 
  geom_point(alpha = 0.5, colour = "gray50") +
  geom_line(aes(y=fitted.fixed), colour = "black", linewidth = 1)+
  facet_grid(Room_Temperature ~ Watering) +
  theme_bw() +
  xlab("Time (min)")+
  ylab(expression(italic(g)[s]~(mol~m^-2~s^-1)))
ggsave("figures/Li6800_gs_on_FittedResponse.png", width=8, height = 8)
ggsave("figures/Li6800_gs_on_FittedResponse.pdf", width=8, height = 8)
ggsave("figures/Li6800_gs_on_FittedResponse.tiff", width=8, height = 8)

li6800.off.gs.plot <- li6800.off.plot
li6800.off.gs.plot <- within(li6800.off.gs.plot, {
  Room_Temperature <- fac.recode(Room_Temperature,
                                 newlevels = c("CER7_Heat",
                                               "CER8_Ambient"))
  Watering <- fac.recode(Watering,
                         newlevels = c("Limited",
                                       "Unlimited"))
  
})
li6800.off.gs.plot$fitted.fixed <- fitted(fit_off_final, re_formula = NA)[,1]

ggplot(li6800.off.gs.plot, aes(tmin, gsw)) + 
  geom_point(alpha = 0.5, colour = "gray50") +
  geom_line(aes(y=fitted.fixed), colour = "black", linewidth = 1)+
  facet_grid(Room_Temperature ~ Watering) +
  theme_bw() +
  xlab("Time (min)")+
  ylab(expression(italic(g)[s]~(mol~m^-2~s^-1)))
ggsave("figures/Li6800_gs_off_FittedResponse.png", width=8, height = 8)
ggsave("figures/Li6800_gs_off_FittedResponse.pdf", width=8, height = 8)
ggsave("figures/Li6800_gs_off_FittedResponse.tiff", width=8, height = 8)

#'### Put together a combined plot
li6800.gs.all.plot <- bind_rows(li6800.on.gs.plot, li6800.off.gs.plot)
li6800.gs.all.plot$Light <- factor(li6800.gs.all.plot$Light,
                                   levels = c("On", "Off"))
li6800.gs.all.plot$fitted.A <- NA
li6800.gs.all.plot$fitted.A[li6800.gs.all.plot$Light == "On"] <- li6800.on.A.plot$fitted.fixed

gs.plot <- ggplot(li6800.gs.all.plot, aes(xElapsed.min, gsw)) + 
  geom_point(alpha = 0.5, colour = "gray50", shape = 17) +
  geom_line(aes(y=fitted.fixed, group = Light, colour = Light), 
            linewidth = 1)+
  geom_vline(aes(xintercept = min(xElapsed.min[Light == "Off"])),
                 colour = "black") +
  geom_vline(data = param.i.A %>% filter(param == "t95"),
             aes(xintercept = emmean+60.5),
             linetype = "dashed") +
  scale_colour_manual(values = c("orange", "black")) +
  facet_grid(Room_Temperature ~ Watering) +
  theme_bw() +
  xlab("Time (min)")+
  ylab(expression(italic(g)[s]~(mol~m^-2~s^-1))) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("(a)")

A.plot <- ggplot(li6800.gs.all.plot, aes(xElapsed.min, A)) + 
  geom_point(alpha = 0.5, colour = "gray50", shape = 16) +
  geom_line(aes(y=fitted.A, group = Light, colour = Light), 
            linewidth = 1)+
  geom_vline(aes(xintercept = min(xElapsed.min[Light == "Off"])),
             colour = "black") +
  geom_vline(data = param.i.A %>% filter(param == "t95"),
             aes(xintercept = emmean+60.5),
             linetype = "dashed") +
  scale_colour_manual(values = c("red", "black")) +
  facet_grid(Room_Temperature ~ Watering) +
  theme_bw() +
  xlab("Time (min)")+
  ylab(expression(italic(A)~(mu*mol~m^-2~s^-1))) +
  theme(legend.position = "none") +
  ggtitle("(b)")

gs.plot / A.plot
ggsave("figures/Li6800_gs.A.fitted.response.png", width=8, height = 8)
ggsave("figures/Li6800_gs.A.fitted.response.pdf", width=8, height = 8)
ggsave("figures/Li6800_gs.A.fitted.response.tiff", width=8, height = 8)

#'### Put together a theoretical plot of the gs and A responses
li6800.theo.plot <- data.frame(Time = seq(from = min(li6800.gs.all.plot$xElapsed.min, na.rm = TRUE),
                                           to = max(li6800.gs.all.plot$xElapsed.min, na.rm = TRUE),
                                           length.out = 500))
li6800.theo.plot$time.adj <- li6800.theo.plot$Time - min(li6800.theo.plot$Time)

##what was the time where the light was switched off?
min(li6800.gs.all.plot$xElapsed.min[li6800.gs.all.plot$Light == "Off"], na.rm = TRUE)
li6800.theo.plot$Light <- "On"
li6800.theo.plot$Light[li6800.theo.plot$Time >= 120.5] <- "Off"

##add A response
li6800.theo.plot$A.resp <-  1.192237 + (24.46462 -  1.192237)*exp(-exp(((-li6800.theo.plot$time.adj/5.957095)+1)))
##what was the mean observed value of A when the light was off?
summary(li6800.gs.all.plot$A[li6800.gs.all.plot$Light == "Off"])
li6800.theo.plot$A.resp[li6800.theo.plot$Light == "Off"] <- 0.02

##add gs responses
li6800.theo.plot$gs.i.resp <-  -0.004365413 + (0.3152114 -  -0.004365413)*exp(-exp((((2.92598-li6800.theo.plot$time.adj)/4.936533)+1)))
li6800.theo.plot$gs.i.resp[li6800.theo.plot$Light == "Off"] <- NA

li6800.theo.plot$gs.d.resp <-  0.33 + (0.02427425 -  0.33)*exp(-exp((((0.6671626-li6800.theo.plot$time.adj)/4.818909)+1)))
li6800.theo.plot$gs.d.resp[li6800.theo.plot$Light == "Off"] <- li6800.theo.plot$gs.d.resp[li6800.theo.plot$Light == "On"]
li6800.theo.plot$gs.d.resp[li6800.theo.plot$Light == "On"] <- NA

##solve the above equation for t95 and then for A95
-5.957095*(log(-log((0.95*24.46462 - 1.192237)/(24.46462-1.192237)))-1) ##23.34524
1.192237 + (24.46462 -  1.192237)*exp(-exp(((-23.34524/5.957095)+1))) ##23.24139

##figure out what Slmax is and where it occurs
(0.3152114 - -0.004365413) / (exp(1) * 4.936533) ##0.02381545
(0.02427425 - 0.33) / (exp(1) * 4.818909) ##-0.02333935
##great, so the maximum rate of change occurs at t = lambda + k
2.92598+4.936533 ##7.862513
0.6671626+4.818909 ##5.486072
##what's the value of the function at this time?
-0.004365413 + (0.3152114 -  -0.004365413)*exp(-exp((((2.92598-7.862513)/4.936533)+1))) ##0.1132003
0.33 + (0.02427425 -  0.33)*exp(-exp((((0.6671626-5.486072)/4.818909)+1))) ## 0.2175298

gs.theo.plot <- ggplot(li6800.theo.plot, aes(Time, gs.i.resp)) + 
  geom_point(colour = "orange", shape = 17) +
  geom_line(colour = "orange") +
  geom_point(aes(y=gs.d.resp), colour = "black", shape = 17) +
  geom_line(aes(y=gs.d.resp), colour = "black") +
  geom_vline(aes(xintercept = min(Time[Light == "Off"])),
             colour = "black") +
  geom_vline(aes(xintercept = 23.34524+60.5), colour = "black", linetype = "dashed") +
  geom_abline(aes(intercept = 0.1132003 - 0.02381545*(7.862513+60.5), slope = 0.02381545),
              linetype = "dotdash", linewidth = 0.75, colour = "brown") +
  geom_abline(aes(intercept = 0.2175298 - -0.02333935*(5.486072+120.6176), slope =-0.02333935),
              linetype = "dotdash", linewidth = 0.75, colour = "grey75") +
  theme_bw() +
  xlab("Time (min)")+
  ylab(expression(italic(g)[s]~(mol~m^-2~s^-1))) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("(a)")

A.theo.plot <- ggplot(li6800.theo.plot, aes(Time, A.resp)) + 
  geom_point(colour = "red", shape = 16) +
  geom_line(colour = "red") +
  geom_vline(aes(xintercept = min(Time[Light == "Off"])),
             colour = "black") +
  geom_vline(aes(xintercept = 23.34524+60.5), colour = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 23.24139), colour = "black", linetype = "dashed") +
  theme_bw() +
  xlab("Time (min)")+
  ylab(expression(italic(A)~(mu*mol~m^-2~s^-1))) +
  theme(legend.position = "none") +
  ggtitle("(b)")

gs.theo.plot / A.theo.plot
ggsave("figures/SchematicResponse.png", width=8, height = 8)
ggsave("figures/SchematicResponse.pdf", width=8, height = 8)
ggsave("figures/SchematicResponse.tiff", width=8, height = 8)

#'### Save data and export
save.image(saveFile)

