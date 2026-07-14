library(brms)

setwd("C:\\Users\\viale002\\OneDrive - Wageningen University & Research\\gs_fit")

# Stomatal conductance model
nlform <- bf(gs ~ g1 + (g0 - g1) * exp(-((th / tau)^theta)),
             g0 ~ 1 + (1 | repl),
             g1 ~ 1 + (1 | repl),
             tau ~ 1,
             theta ~ 1,
             nl = TRUE)

priors <- c(prior(normal(0.05, 0.02), nlpar = "g0"),
            prior(normal(0.5, 0.2), nlpar = "g1"),
            prior(normal(0, 0.05), class = "sd", nlpar = "g0"),
            prior(normal(0, 0.1), class = "sd", nlpar = "g1"),
            prior(normal(1, 1), nlpar = "tau", lb = 0),
            prior(normal(2, 2), nlpar = "theta", lb = 0),
            prior(student_t(3, 0, 0.01), class = "sigma"))

process <- function(d1, fit_template) {
    dd <- subset(d1, PAR == 1 & !is.na(gs))

    # t=0 at the light step 
    dd$th <- dd$th - dd$th[1]

    fit <- update(
        fit_template,
        newdata = dd,
        recompile = FALSE,
        refresh = 0
    )

    # Check diagnostics, Important!
    print(summary(fit))

    #plot(fit)
    #pp_check(fit)

    plot(gs~th, dd, pch=16)
    points(fitted(fit)[,1]~th, dd, pch=16, col="red")

    return(fit)
}

d <- read.csv("McAusland_compiled.csv")
d$th <- d$time / 3600

# pick a small but valid dataset
dd_init <- subset(d, PAR == 1 & plant == d$plant[1])
dd_init$th <- dd_init$th - dd_init$th[1]

# Configure and compile the model
fit_init <- brm(
  formula = nlform,
  data = dd_init,
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 2000,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 13,
    metric = "diag_e"
  )
)

# Main loop that select each plant per species
fits <- list()

for (pl in unique(d$plant)) {
  d1 <- subset(d, plant == pl)
  cat("Plant:", pl, "\n")
  fits[[pl]] <- process(d1, fit_init)
}