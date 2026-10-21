#'# EXP 801 DS Wheat Heat Drought
#'## --- Ed Chaplin (USyd), *Apr--June 2025* ---
#'### Li-600 STOMATAL CONDUCTANCE MODEL PARAMETERS VS ANATOMICAL TRAIT CORRELATIONS
#'
#' * Using the parameters of interest estimated from the stomatal conductance models and the anatomical trait predictions, calculate correlations on a 'environment' basis
#' 
#' * Changes to v2: Adding multiple comparison adjustment to p-values in correlation matrix
#'
#'### Load the libraries
library(knitr)
#knitr::spin("Exp801.porometer.stomatalCondvsAnatomicalTraits.v2.r")
library(ggplot2)
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, 
               "fig.width" = 7.5, "fig.height" = 10, "comment"=NA, fig.path = "figures/increase/")
library(dplyr)
library(tidyverse)
library(lattice)
library(GGally)

#'### Set up constants and files
options(width=90)
load("../data/stomatalCondParamEstimation_porometer_i.v1.Rdata")
load("../data/stomatalCondParamEstimation_porometer_d.v1.Rdata")
load("../../../Anatomical traits/data/anat.diffs.full.rda")
saveFile <- "../data/stomatalCondvsAnatTraits.v2.Rdata"

#'### Add variables to make the dataframes compatible
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
head(traits.df)
pairs(traits.df[,4:ncol(traits.df)])

#'### Put together proper plot
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

## Create a custom function to calculate corrected p-values for the plot matrix
custom_cor_with_adjustment <- function(data, mapping, method = "bonferroni") {
  
  # Identify the current x and y variables from the mapping
  x_var <- rlang::as_label(mapping$x)
  y_var <- rlang::as_label(mapping$y)
  
  # Select only the numeric columns that are likely being used in the pairs plot
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  
  # Step A: Compute the raw p-value matrix for all numeric columns
  p_matrix <- matrix(NA, nrow = length(numeric_cols), ncol = length(numeric_cols),
                     dimnames = list(numeric_cols, numeric_cols))
  
  for (i in numeric_cols) {
    for (j in numeric_cols) {
      if (i != j) {
        # Extract pair complete observations to match ggpairs default behavior
        valid_data <- na.omit(data[, c(i, j)])
        if (nrow(valid_data) > 2) {
          p_matrix[i, j] <- suppressWarnings(cor.test(valid_data[[i]], valid_data[[j]])$p.value)
        }
      }
    }
  }
  
  # Step B: Flatten the matrix to adjust all p-values simultaneously
  raw_p_values <- as.vector(p_matrix)
  adjusted_p_values <- p.adjust(raw_p_values, method = method) # e.g., "BH", "bonferroni"
  
  # Reconstruct back into the matrix format
  adjusted_p_matrix <- matrix(adjusted_p_values, nrow = nrow(p_matrix), 
                              dimnames = dimnames(p_matrix))
  
  # Step C: Extract the correlation coefficient and the specific adjusted p-value
  current_cor <- suppressWarnings(cor(data[[x_var]], data[[y_var]], use = "complete.obs"))
  current_adj_p <- adjusted_p_matrix[x_var, y_var]
  ##build dataframe to add significance stars
  sig_adj_p <- data.frame("adj.p" = current_adj_p)
  sig_adj_p$sig <- ""
  sig_adj_p$sig[sig_adj_p$adj.p < 0.10] <- "."
  sig_adj_p$sig[sig_adj_p$adj.p < 0.05] <- "*"
  sig_adj_p$sig[sig_adj_p$adj.p < 0.01] <- "**"
  sig_adj_p$sig[sig_adj_p$adj.p < 0.001] <- "***"
  
  # Step D: Construct the text display (similar to ggpairs default style)
  label_text <- paste0(
    "Corr: ", "\n", format(round(current_cor, 2), nsmall = 2), sig_adj_p$sig
  )
  
  # Create a blank plot centered with our custom text label
  ggplot(data = data, mapping = mapping) +
    annotate("text", x = mean(range(data[[x_var]], na.rm = TRUE)), 
             y = mean(range(data[[y_var]], na.rm = TRUE)), 
             label = label_text, size = 4) +
    theme_void() +
    theme(panel.background = element_rect(fill = "white", color = "grey90"))
}


p1 <- ggpairs(
  filter(traits.df, Room_Temperature == "CER7_Heat" &
           Device_Watering == "D1_Unlimited"),
  columns = traits,
  # upper = list(continuous = wrap("cor", stats = TRUE,
  #                                digits = 2, size = 4,
  #                                colour = "black")),
  upper = list(continuous = wrap(custom_cor_with_adjustment, method = "bonferroni")),
  lower = list(continuous = lower_scatter),
  diag = list(continuous = diag_lab),
  progress = FALSE) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )
p1 <- p1 + ggtitle("CER7_Heat, D1_Unlimited")
p1
ggsave("figures/CER7_Heat_D1_Unlimited_Correlation.png", width = 10, height = 10)
ggsave("figures/CER7_Heat_D1_Unlimited_Correlation.pdf", width = 10, height = 10)
ggsave("figures/CER7_Heat_D1_Unlimited_Correlation.tiff", width = 10, height = 10)

p2 <- ggpairs(
  filter(traits.df, Room_Temperature == "CER7_Heat" &
           Device_Watering == "D2_Limited"),
  columns = traits,
  # upper = list(continuous = wrap("cor", stats = TRUE, 
  #                                digits = 2, size = 4,
  #                                colour = "black")),
  upper = list(continuous = wrap(custom_cor_with_adjustment, method = "bonferroni")),
  lower = list(continuous = lower_scatter),
  diag = list(continuous = diag_lab),
  progress = FALSE) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )
p2 <- p2 + ggtitle("CER7_Heat, D2_Limited")
suppressWarnings(p2)
ggsave("figures/CER7_Heat_D2_Limited_Correlation.png", width = 10, height = 10)
ggsave("figures/CER7_Heat_D2_Limited_Correlation.pdf", width = 10, height = 10)
ggsave("figures/CER7_Heat_D2_Limited_Correlation.tiff", width = 10, height = 10)

p3 <- ggpairs(
  filter(traits.df, Room_Temperature == "CER8_Ambient" &
           Device_Watering == "D1_Unlimited"),
  columns = traits,
  # upper = list(continuous = wrap("cor", stats = TRUE, 
  #                                digits = 2, size = 4,
  #                                colour = "black")),
  upper = list(continuous = wrap(custom_cor_with_adjustment, method = "bonferroni")),
  lower = list(continuous = lower_scatter),
  diag = list(continuous = diag_lab),
  progress = FALSE) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )
p3 <- p3 + ggtitle("CER8_Ambient, D1_Unlimited")
suppressWarnings(p3)
ggsave("figures/CER8_Ambient_D1_Unlimited_Correlation.png", width = 10, height = 10)
ggsave("figures/CER8_Ambient_D1_Unlimited_Correlation.pdf", width = 10, height = 10)
ggsave("figures/CER8_Ambient_D1_Unlimited_Correlation.tiff", width = 10, height = 10)

p4 <- ggpairs(
  filter(traits.df, Room_Temperature == "CER8_Ambient" &
           Device_Watering == "D2_Limited"),
  columns = traits,
  # upper = list(continuous = wrap("cor", stats = TRUE, 
  #                                digits = 2, size = 4,
  #                                colour = "black")),
  upper = list(continuous = wrap(custom_cor_with_adjustment, method = "bonferroni")),
  lower = list(continuous = lower_scatter),
  diag = list(continuous = diag_lab),
  progress = FALSE) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )
p4 <- p4 + ggtitle("CER8_Ambient, D2_Limited")
suppressWarnings(p4)
ggsave("figures/CER8_Ambient_D2_Limited_Correlation.png", width = 10, height = 10)
ggsave("figures/CER8_Ambient_D2_Limited_Correlation.pdf", width = 10, height = 10)
ggsave("figures/CER8_Ambient_D2_Limited_Correlation.tiff", width = 10, height = 10)

#'### Save data and export
save.image(saveFile)

