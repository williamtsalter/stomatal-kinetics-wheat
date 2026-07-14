#'# Wheat stomatal kinetics (heat & drought)
#'## --- Ed Chaplin (USyd), *Apr--June 2025* ---
#'
#'### Analysis of anatomical traits
#' 

#'### Load packages
library(knitr)
#knitr::spin("anat.anout.v1.r")
library(ggplot2)
library(asreml)
library(asremlPlus)

#'### Set up constants and files
options(width=95)
load("../data/anat.Analyses.v1.RData")
source("anat.anal.v1.stub.r")

#'## Show the results
#+ "show-v1"
for (resp in responses.anat[kvec]) {
  showIndvModelResults(models[[resp]], response=resp)
}

#'## MODEL SUMMARIES
#'### Wald tables
for (resp in responses.anat[kvec]) {
  tbl <- models[[resp]]$wald.tab[-1,]
  cat("\n\n", resp)
  print(kable(tbl, digits = 3))
}

#'### Chosen terms
tableChosenModels(models, responses.anat)

#'### Table the residual variances
tableVarModels(models, responses.anat)

#'## Extremal residuals
#+ "outliers-v1"
typelabls <- c("Anat")
resptlabs <- rep(typelabls,
                 times = length(responses.anat))
print(cbind(responses.anat, resptlabs), quote=F)
res.lim.label <- 5
Outliers.dat <- data.frame()
tbl <- data.frame()
for (resp in responses.anat[kvec]) 
{
  if(resp %in% responses.pot){
    tmp.dat <- models[[resp]]$asreml.obj$call$data
    res.dat <- within(tmp.dat, res.std <- residuals(models[[resp]]$asreml.obj, type="stdCond"))
    res.dat <- res.dat[!is.na(res.dat$res.std), ]
    outliers.dat <- subset(res.dat, abs(res.std) > res.lim.label)
    outliers.dat <- outliers.dat[c("Cell.ID", "res.std", "Rm_Temperature", 
                                   "Watering", "Line")]
    if (nrow(outliers.dat) > 0)
    { 
      Outliers.dat <- rbind(Outliers.dat,
                            data.frame(Response=resp, outliers.dat))
      tbl <- rbind(tbl,
                   data.frame(Response = resp, 
                              n        = nrow(outliers.dat),
                              rmin     = min(res.dat$res.std, na.rm=T),
                              rmax     = max(res.dat$res.std, na.rm=T)))
    }
  } else{
    tmp.dat <- models[[resp]]$asreml.obj$call$data
    res.dat <- within(tmp.dat, res.std <- residuals(models[[resp]]$asreml.obj, type="stdCond"))
    res.dat <- res.dat[!is.na(res.dat$res.std), ]
    outliers.dat <- subset(res.dat, abs(res.std) > res.lim.label)
    outliers.dat <- outliers.dat[c("Cell.ID", "Stomata", "res.std", "Rm_Temperature", 
                                   "Watering", "Line")]
    if (nrow(outliers.dat) > 0)
    { 
      Outliers.dat <- rbind(Outliers.dat,
                            data.frame(Response=resp, outliers.dat))
      tbl <- rbind(tbl,
                   data.frame(Response = resp, 
                              n        = nrow(outliers.dat),
                              rmin     = min(res.dat$res.std, na.rm=T),
                              rmax     = max(res.dat$res.std, na.rm=T)))
    }
  }
}

if (nrow(Outliers.dat)==0) {
  
  cat("\n\n #### No outliers found\n\n")
  
} else {
  
  Outliers.dat <- with(Outliers.dat, Outliers.dat[order(Cell.ID, Stomata, Response),])
  Outliers.dat <- within(Outliers.dat, res <- round(res.std,2))
  tmp.dat <- Outliers.dat[c("Cell.ID","Stomata", "res","Response", "Rm_Temperature", 
                            "Watering", "Line")]
  names(tmp.dat)[1] <- c("Plant")
  outliers.dat <- unique(tmp.dat[,-(3:4)])
  tbl <- within(tbl, {
    rmin <- ifelse(rmin < 0, round(rmin,2), NA)
    rmax <- ifelse(rmax > 0, round(rmax,2), NA)
  })
  
  msg1 <- paste("Standardised residuals outside the range (-",
                res.lim.label, ", ",res.lim.label,")",
                sep="")
  msg2 <- "Unique plants in outlier list"
  
  cat("\n\n #### Number of outliers per response\n\n"); print(tbl, right=F)
  cat("\n\n",msg1,"\n\n"); print(tmp.dat, right=F)
  cat("\n\n",msg2,"\n\n"); print(outliers.dat, right=F)
  
  tmp2.dat <- within(tmp.dat, {
    Magnitude <- abs(res)
    Direction <- factor(sign(res), levels = c(1,-1), labels = c("pos","neg"))
    Response <- factor(Response, levels = responses.anat)
    Resptype <- dae::fac.recode(Response, newlevels = resptlabs, levels = rev(typelabls))
  })
  
  lims <- c(res.lim.label,
            max(4, 0.5*ceiling(max(tmp2.dat$Magnitude)/0.5)))  #round to nearest half
  plt <- ggplot(tmp2.dat, aes(Plant, as.character(Response),
                              shape=Direction, colour=Magnitude, fill=Magnitude)) +
    facet_grid(rows = vars(Resptype), space="free", scales="free") +
    geom_point(size=2) +
    scale_colour_gradient(low="green", high="red", limits = lims) +
    scale_fill_gradient(low="green", high="red", limits = lims) +
    scale_shape_manual(values = c(24,25), drop=F) +
    labs(y=NULL) +
    scale_y_discrete(drop=F) +
    ggtitle("Outliers by Plant and Response",
            subtitle = "Coloured by magnitude of standardized residual") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
          strip.text = element_text(size=10, face="bold"))
  
  print(plt)
}

