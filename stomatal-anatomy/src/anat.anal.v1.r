#'# Wheat stomatal kinetics (heat & drought)
#'## --- Ed Chaplin (USyd), *Apr--June 2025* ---
#'
#'### Analysis of anatomical traits
#' 

#'### Load packages
library(knitr)
#knitr::spin("anat.anal.v1.r")
library(ggplot2)
library(asreml)
library(asremlPlus)
library(parallel) #for parallel processing
library(foreach)
library(doParallel)
packageVersion("asreml")
packageVersion("asremlPlus")

#'### Set up constants and files
options(width=95)
source("anat.anal.v1.stub.r")
load("../data/anat.dat.v1.Rdata")
saveFile <- "../data/anat.Analyses.v1.RData"

#+ "anal-v1"
#'## Initialize parallel processing
ncores <- detectCores()
cl <- makeCluster(ncores)
registerDoParallel(cl)

responses.pot <- responses.anat[c(1, 2, 7)]
responses.stoma <- responses.anat[!(responses.anat %in% responses.pot)]

# Fit the models to the pot level first
#+ "anal-pot-v1"
kvec <- 1:length(responses.anat)
models <- foreach(response = responses.anat[kvec], 
                  .errorhandling = "pass", 
                  .inorder=TRUE,
                  .packages = c("asreml", "asremlPlus", "dae"))  %dopar%
  {
    if(response %in% responses.pot){
      mod <- fitIndvModel(response=response, selection = "backward", 
                          fix.mod = ~ ASide/ABlock + Rm_Temperature/(AWatering*ALine + HWatering*HLine),
                          ran.mod = ~  Room:Side/Location/SubUnit,
                          term.marginality = marginality(mod.canon)$trts,
                          classify = c("Rm_Temperature:AWatering:ALine:HWatering:HLine"),
                          add.facs2preds = c("Watering", "Line", "Line.ID", "WT"),
                          hetero.facs =  c("Watering", "Rm_Temperature"),
                          dat = anat.pot.dat)
    } else{
      mod <- fitIndvModel(response=response, selection = "backward", 
                          fix.mod = ~ ASide/ABlock + Rm_Temperature/(AWatering*ALine + HWatering*HLine),
                          ran.mod = ~  Room:Side/Location/SubUnit/Cell,
                          term.marginality = marginality(mod.canon)$trts,
                          classify = c("Rm_Temperature:AWatering:ALine:HWatering:HLine"),
                          add.facs2preds = c("Watering", "Line", "Line.ID", "WT"),
                          hetero.facs =  c("Watering", "Rm_Temperature"),
                          dat = anat.stoma.dat)
    }
  }
names(models) <- responses.anat[kvec]

msgs <- (do.call(rbind, lapply(models, 
                               function(mod) if ("message" %in% names(mod)) return(mod$message))))
if (!is.null(msgs))
{
  colnames(msgs) <- "messages"
  cat("\n#### Output any error messages\n\n")
  print(msgs)
}

#'## Hand back clusters
stopCluster(cl)

#'## Save workspace - the compression takes longer, but produces a substantially smaller file
save.image(saveFile, compress = "xz")

