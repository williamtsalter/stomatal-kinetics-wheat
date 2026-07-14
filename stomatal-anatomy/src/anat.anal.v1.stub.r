#'# Wheat stomatal kinetics (heat & drought)
#'## --- Ed Chaplin (USyd), *Apr--June 2025* ---
#'## ANALYSIS of individual traits
#'
#'### Changes in versionn 2:
#'
#' * Change Temperature to Rm_Temperature.
#' 

#'### Load the libraries
library(knitr)
#knitr::spin("indv.anal.v2.stub.r")
library(ggplot2)
library(asreml)
library(asremlPlus)
library(qqplotr)
library(gridExtra)
library(dae)
source("anal.stub.funcs.r")

#'### Set up constants
options(width=95)
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, fig.path = "figures/",
               "fig.width" = 7.5, "fig.height" = 10)


#'## Function to show the results of the analysis of the `response` or $k$th trait
showIndvModelResults <- function(model, k=0, response=NULL)
{
  if (k > nresp) {
    warning(paste0("Unrecognised trait 'responses[", k,"]'.\n"))
    k <- 0
  }
  if (is.null(response))
    response <- ifelse(k > 0,  responses[k], "TRAIT")
  else
    k <- match(response,responses)
  
  dat  <- model$asreml.obj$call$data
  pwr  <- model$pwr
  offs <- model$offs
  transformed <- pwr!=1 | offs != 0

  response0 <- response # pre-transformed trait
  
  if (transformed) # transformed trait
    response  <- names(dat)[length(names(dat))] # must have been appended to DAT
  
  #Output results from final analysis
  cat("\n\n#### Results for", response, "\n\n")
  
  #Check for an error
  if ("message" %in% names(model)) 
  {    
    mess <- model$message
    names(mess) <- "messages"
    print(mess)
    print(model$call)
  } else
  {
    current.asrt <- as.asrtests(model$asreml.obj, model$wald.tab, model$tests)
    if (current.asrt$asreml.obj$converge) {
      cat("\n  Final model converged \n\n")
    } else {
      cat("\n  Final model failed to converge \n\n")
    }
    printFormulae(current.asrt$asreml.obj)
    print(current.asrt)
    
    cat("\n  Information criteria for spatial models \n\n")
    print(model$spatial.IC)
    
    cat("\nChosen terms\n")
    print(model$chosen.terms)  
    
    cat("\nLSDs\n")
    print(model$diffs$LSD)
    
    preds.dat <- model$preds
    LSIs <- "lower.halfLeastSignificant.limit" %in% names(preds.dat)
    if (LSIs) {
      titl <- "Half-LSIs"
      errorBar <- geom_errorbar(aes(ymin=lower.halfLeastSignificant.limit,
                                    ymax=upper.halfLeastSignificant.limit),
                                width=0.5, alpha=0.75)
    } else {
      titl <- "95% CIs"
      errorBar <- geom_errorbar(aes(ymin=lower.Confidence.limit,
                                    ymax=upper.Confidence.limit),
                                width=0.5, alpha=0.75)
    }
    pred.plt <- ggplot(preds.dat, aes(x=Line.ID, y=predicted.value, fill=WT)) +
      facet_grid(cols = vars(Rm_Temperature), rows = vars(Watering)) +
      geom_bar(stat="identity", alpha=0.75) +
      errorBar +
      labs(x=NULL, y = "Predicted value") +
      ggtitle(paste(response0,titl,sep=" -- ")) +
      theme_bw() +
      theme(legend.position = "none",
            strip.text = element_text(face = "bold", size=11),
            axis.text.x = element_text(size=9, hjust=1, vjust=0.5, angle=90))

    # Diagnostic plots of residuals
    res.dat <- getDiagnosticEffFit(model$asreml.obj,
                                   classify.facs = c("Rm_Temperature", "Watering", "Line.ID"),
                                   preds = preds.dat)
    res.dat <- res.dat[!is.na(res.dat$res.std), ]
    res.dat <- within(res.dat, { out.lab <- substring(Cell.ID,4)})

    plts <- plotResiduals(res.dat, 
                          outlier.label = "out.lab",
                          response = response,
                          colour.col = "WT",  xdodge = TRUE,
                          ggplotFuncs = list(facet_grid(rows = vars(Watering)),
                                             theme(legend.position = "bottom")),
                          var.mod = model$var.mod)


    suppressWarnings(grid.arrange(plts$resraw.plt, plts$qqplt, ncol=1))
    print(plts$resstd.plt)
    print(pred.plt)
   }
  
  invisible(model)
}

