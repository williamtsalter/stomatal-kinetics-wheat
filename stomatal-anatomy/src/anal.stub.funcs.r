#'# Reusable functions for use with the indv.anal/anout scripts.
#'
#'* The defaults for these functions must not be changed.
#'* Except for showIndvModelResults, changes must not be made to these functions without approval of CB.
#'
#'* The following functions are currently included:
#'
#'    + fitIndvModel, showIndvModelResults;
#'    + choose1or2WayHetVarModel, convertRownames2Column, fixBoundVCs, getDiagnosticEffFit, 
#'      getVarModel, plotResiduals, addR2adj2model;
#'    + tableAICR2ASRemlobjs, tableAICR2SpatialICs, tableChosenModels, tableEMMstats, 
#'      tableResVars, tableSEstats, tableVarModels.  
#'    

#'## Example heading and introduction for an indv.anal script
#'# EXP 589b Barley Drought
#'## --- Poznan, *Winter 2021* ---
#'## ANALYSIS of indv data from a split-plot design in a single Smarthouse
#'
#'### Changes in Version v2:
#'
#'* Updated to use the functions in `anal.stub.funcs.r`
#'

#'### Load the libraries
library(knitr)
#knitr::spin("anal.stub.funcs.r")
library(ggplot2)
library(asreml)
library(asremlPlus)
library(dae)
library(dplyr) # residual plot: summary statistics
library(stringr) #for splitting terms into factors
library(qqplotr)
library(gridExtra)
packageVersion("asreml")
packageVersion("asremlPlus")
packageVersion("dae")
#source("anal.stub.funcs.r")

#'### Set up constants
options(width=95)
opts_chunk$set("tidy" = FALSE, "dev" = "png", "dpi" = 300, fig.path = "figures/",
               "fig.width" = 7.5, "fig.height" = 10)

#'# The functions 

#'## Analyse `response` or $k$th trait
#'
#'### This function requires the following functions: choose1or2WayHetVarModel, fixBoundVCs 
#'   and getVarModel.
#'   
#'### It assumes that responses are available
#'
fitIndvModel <- function(
    k=1, response=NULL, resp.vec = responses, 
    fix.mod = ~ NULL,   #supply a RHS formula only by replacing NULL
    ran.mod = ~ NULL,   #supply a formula by replacing NULL
    res.mod = ~ NULL,   #supply a formula by replacing NULL
    pwr=1, offs=0, 
    hetero.facs = NULL, #a character or a list of a characters, 
    #                    each character having 1 or 2 factor names
    selection = "backward",
    fix.bound.val = 1e-04,   #value for fixed bound variance components
    term.marginality = NULL, #use marginality function to extract from anatomy;
    #                        if NULL, chosen.terms <- NULL
    classify,                #must supply a classify
    present = NULL,          #a character vector giving the present set of variables for 
    #                         which predictions are to be based on cambinations with data. 
    add.facs2preds = NULL,   # a vector of factors in data that are to be added to 
    #                        the predictions.frames (used with nested models)
    LSDby=NULL, useLSI=TRUE, 
    spatial.models = "none",      #change to "all" to choose a spatial model
    sections = NULL,  #a factor for which separate separate spatial models are to be fitted 
    spatial.row.fac = "Lane",     #for row.covar, uses or makes paste0("c", spatial.row.fac)
    spatial.col.fac = "Position", #for col.covar, uses or makes paste0("c", spatial.col.fac) 
    dat = indv.dat, 
    asremlArgs = NULL)     # a list of the form list(arg1 = value1, arg2 = value2, ..)
{ 
  if (is.null(response))  response <- resp.vec[k]  else  k <- match(response,resp.vec)
  class.facs <- stringr::str_split_1(classify, ":") #used for factors in present argument
  
  transformed <- pwr!=1 | offs!=0
  
  if (transformed) {
    #stop if any negative or, for log transformation, zeroes in the data, after any offset 
    #has been applied
    tmp  <- dat[[response]]
    bad  <- !is.na(tmp) & (tmp < -offs | (pwr==0 & tmp == -offs)) 
    nBad <- length(bad[bad])
    if (nBad > 0) {
      stop(paste("Found", nBad,
                 ifelse(pwr==0, "non-positive", "negative"),
                 "observation(s)."))
    }
    names(response) <- response # enable automated naming of transformed variables
    trans.obj <- powerTransform(var.name = response, power=pwr, offset = offs, 
                                titles = response, data = dat)
    dat <- trans.obj$data
    response <- trans.obj$tvar.name
  }
  
  #Add response to fixed model and fit the inital model
  rlang::f_lhs(fix.mod)<-as.symbol(response)
  cat("\n\n#### fixed model\n")
  print(fix.mod)
  asreml.options(keep.order = TRUE)
  current.asr <- do.call(asreml, 
                         args=c(list(fixed = fix.mod,
                                     random = ran.mod,
                                     residual = res.mod,
                                     data = dat, maxiter=50),
                                asremlArgs))
  current.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                              label = "Initial model")
  
  # Test for unequal variance
  if (is.list(hetero.facs) && length(hetero.facs) == 1)
    hetero.facs <- hetero.facs[[1]]
  if (is.list(hetero.facs) && length(hetero.facs) > 1)
  { 
    for (het.facs in hetero.facs)
    {  
      current.asrt <- choose1or2WayHetVarModel(current.asrt, selection = selection,
                                               hetero.term = "Residual", hetero.facs = het.facs)
    }
    hetero.facs <- unlist(hetero.facs)
  } else
    current.asrt <- choose1or2WayHetVarModel(current.asrt, selection = selection,
                                             hetero.term = "Residual", hetero.facs = hetero.facs)
  var.model <- getVarModel(current.asrt$asreml.obj, hetero.facs = hetero.facs)
  
  # remove boundaries, if any
  current.asrt <- rmboundary(current.asrt)
  
  #Fix any bound variance components to fix.bound.val 
  current.asrt <- fixBoundVCs(current.asrt, fix.bound.val = fix.bound.val)
  
  if (!("none" %in% spatial.models))
  { 
    #Get nonspatial diffs for comparison with spatial diffs
    if (is.null(present)) present <- class.facs
    nonspatial.diffs <- predictPlus(current.asrt$asreml.obj,
                                    classify = classify, present = present, 
                                    transform.power = pwr, offset = offs, 
                                    wald.tab = current.asrt$wald.tab, 
                                    tables = "none", Vmatrix = TRUE)
    if (!is.null(add.facs2preds))
      nonspatial.diffs$predictions <- dplyr::left_join(nonspatial.diffs$predictions, 
                                                       unique(dat[c(class.facs, add.facs2preds)]))
    
    #Get nonspatial denDF for the treatment terms  
    nspatdDF <- current.asrt$wald.tab$denDF
    names(nspatdDF) <- rownames( current.asrt$wald.tab)
    
    #Set up numeric covariates
    spatial.facs <- c(spatial.row.fac, spatial.col.fac) 
    if (!all(spatial.facs %in% names(dat)))
      stop(paste0("The following spatial.facs are not in the data.frame supplied via dat: ",
           spatial.facs[!(spatial.facs %in% names(dat))]))
    spatial.covs <- paste0("c", spatial.facs)
    if (!all(spatial.covs %in% names(dat)))
    {
      dat[spatial.covs[1]] <- as.numfac(unlist(dat[spatial.facs[1]]), center = TRUE)
      dat[spatial.covs[2]] <- as.numfac(unlist(dat[spatial.facs[2]]), center = TRUE)
      current.asrt$asreml.obj$call$data <- dat; if (is.symbol(dat)) dat <- eval(dat)
    }
    
    #Try the spatial models
    spatial.fit <- chooseSpatialModelOnIC(current.asrt, 
                                          sections = sections,
                                          row.covar = spatial.covs[1], col.covar = spatial.covs[2], 
                                          row.factor = spatial.facs[1], col.factor = spatial.facs[2], 
                                          trySpatial = spatial.models, 
                                          dropFixed = paste(spatial.covs, collapse = " + "), 
                                          dropRandom = paste(sapply(spatial.covs, 
                                                              function(cov) 
                                                                paste0("spl(", cov, ") + dev(", 
                                                                       cov, ")")),
                                                             collapse = " + "),
                                          allow.fixedcorrelation = FALSE, 
                                          allow.unconverged = FALSE,
                                          asreml.option = "grp", 
                                          rotateX = TRUE, ngridangles = NULL,
                                          return.asrts = "best")
    current.asrt <- spatial.fit$asrts[[1]]
    
    #Use Residual DF from the nonspatial model for denDF that are NA for the treatment terms  
    dDF <- current.asrt$wald.tab$denDF
    knames <- rownames(current.asrt$wald.tab)
    names(dDF) <- knames
    knames <- names(nspatdDF)[names(nspatdDF) %in% knames]
    dDF[knames] <- nspatdDF[knames]
    dDFfault <- "supp"
    current.asrt$wald.tab <- recalcWaldTab(current.asrt, 
                                           dDF.fault = dDFfault, dDF.values = dDF)
  } else
  {  
    nonspatial.diffs <- NULL
    dDFfault <- "none"
    dDF <- NULL
  }
  
  # test treatment terms
  if (!is.null(term.marginality))
  {
    choice <- chooseModel(object = current.asrt, 
                          terms.marginality = term.marginality, 
                          denDF="numeric", ssType = "conditional", 
                          dDF.fault = dDFfault, dDF.values = dDF)
    current.asrt <- choice$asrtests.obj
    chosen.terms <- unlist(choice$sig.terms)
  } else
    chosen.terms <- NULL
  
  if (is.null(present)) present <- class.facs
  if (useLSI) {
    diffs <- predictPlus(current.asrt$asreml.obj,
                         classify = classify, present = present, 
                         transform.power = pwr, offset = offs, 
                         wald.tab = current.asrt$wald.tab, 
                         error.intervals = "half", 
                         LSDtype=ifelse(is.null(LSDby), "overall", "factor"), 
                         LSDby=LSDby,
                         avsed.tolerance = NA,
                         tables = "none", Vmatrix = TRUE)
    cat("\n")
  } else {
    diffs <- predictPlus(current.asrt$asreml.obj,
                         classify = classify, present = present, 
                         transform.power = pwr, offset = offs, 
                         wald.tab = current.asrt$wald.tab, 
                         tables = "none", Vmatrix = TRUE)
  }
  
  if (!is.null(add.facs2preds))
    diffs$predictions <- dplyr::left_join(diffs$predictions, 
                                          unique(dat[c(class.facs, add.facs2preds)]))
  
  
  if (transformed) {
    preds.dat <- diffs$backtransforms
    preds.dat <- dplyr::left_join(preds.dat, unique(dat[c(class.facs, add.facs2preds)]))
    names(preds.dat)[match("backtransformed.predictions", names(preds.dat))] <- "predicted.value"
  } else
    preds.dat <- diffs$predictions
  
  #Add standardized residuals
  asreml.obj <- current.asrt$asreml.obj
  asreml.obj <- update(asreml.obj, aom = TRUE)
  # models testing results and store summaries
  model <- list()
  model[["asreml.obj"]] <- asreml.obj
  model[["wald.tab"]] <- current.asrt$wald.tab
  model[["tests"]] <- current.asrt$test.summary
  model[["chosen.terms"]] <- chosen.terms
  model[["diffs"]] <- diffs
  model[["preds"]] <- preds.dat
  model[["var.mod"]] <- var.model
  model[["pwr"]] <- pwr
  model[["offs"]] <- offs
  if (!("none" %in% spatial.models))
  { 
    model[["spatial.IC"]] <- spatial.fit$spatial.IC
    model[["best.spatial.mod"]] <- spatial.fit$best.spatial.mod
    model[["best.spatial.AIC"]] <- spatial.fit$best.spatial.IC
    model[["nonspatial.preds"]] <- nonspatial.diffs$predictions
  }
  
  invisible(model)
}


#'## Function to show the results of the analysis of the `response` or $k$th trait
showIndvModelResults <- function(model, k=0, response=NULL, resp.vec = responses)
{
  if (k > nresp) {
    warning(paste0("Unrecognised trait 'resp.vec[", k,"]'.\n"))
    k <- 0
  }
  if (is.null(response))
    response <- ifelse(k > 0,  resp.vec[k], "TRAIT")
  else
    k <- match(response,resp.vec)
  
  dat  <- model$asreml.obj$call$data; if (is.symbol(dat)) dat <- eval(dat)
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
    
    sIC <- model$spatial.IC
    if (!is.null(sIC)) {
      cat("\n  Information criteria for spatial models \n\n")
      print(sIC)
    }
    
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
    pred.plt <- ggplot(preds.dat,
                       aes(x=Soil, y=predicted.value, fill=Soil)) +
      facet_grid(cols = vars(Watering)) +
      geom_bar(stat="identity", alpha=0.75) +
      errorBar +
      labs(x=NULL, y = "Predicted value") +
      ggtitle(paste(response0,titl,sep=" -- ")) +
      theme_bw() +
      theme(legend.position = "none",
            strip.text = element_text(face = "bold", size=11),
            axis.text.x = element_text(size=9, hjust=1, vjust=0.5, angle=90))
    
    # Get the residuals/effects and fitted/EMM values for diagnostic plots
    if (!is.null(model$best.spatial.mod) && model$best.spatial.mod == "corr")
      res.dat <- getDiagnosticEffFit(model$asreml.obj, terms = "Lane:Position", use = "G.aom", 
                                     classify.facs = c("Salinity","Line"), 
                                     preds = model$diffs$predictions)
    else    
      res.dat <- getDiagnosticEffFit(model$asreml.obj, 
                                     classify.facs = c("Salinity","Line"), 
                                     preds = model$diffs$predictions)

    soil.levs <- levels(res.dat$Soil)
    soil.levs <- gsub(" ", "", soil.levs); soil.levs <- gsub("Lucerne", "Luc", soil.levs)
    soil.levs <- gsub("Aquasil", "Asil", soil.levs)
    res.dat <- within(res.dat, Soil <-  fac.recast(Soil, newlevels = soil.levs, 
                                                   levels.order = soil.levs))
    #Do individual residual plots
    plts <- plotResiduals(res.dat, x = "fit", x.facet = "Watering", 
                          outlier.label = "Pot.ID", 
                          response = response, colour.col = "Soil",  
                          var.mod = model$var.mod)

    suppressWarnings(grid.arrange(pred.plt, qplt, ncol=1))
    
    grid.arrange(resraw.plot, resfit.plot, ncol=1)
  }
  
  invisible(model)
}

#'### Function to test hetero.term for one-way heterogeneity for a single factor (fac1) or two-way heterogeneity for two factors (fac1 and fac2) in a list component of hetero.facs
#'
#' * hetero.facs can be a character or a list of a characters, with each character having 
#'   1 or 2 factor names. If a single character vector then heterogeneity for the factor(s) 
#'   will be investigated. If the list has more than one component, then heterogeneity will 
#'   be investigated for each component sequentially. 
#'   Thus, hetero.facs = list(c("Line", "Salinity")) is a single-component list and the two-way 
#'   heterogeneity will be explored for the two factors. On the other hand, 
#'   hetero.facs = list("Line", "Salinity") is a two-component list that specifies sequential. 
#'   testing of the models. In this case, the heterogeneity of Line will be investigated first 
#'   and then the heterogeneity of Salinity will be investigated. That is, for a two-component 
#'   list, the only possible outcomes are one-way or additive heterogeneity. It is likely that 
#'   `selection = "forward"` should be used so that heterogeneity for the second factor is 
#'   compared with the heterogeneity for the first factor. Othewrise, heterogeneity for the 
#'   second factor will be compared with the homogeneous model. 
#' * If hetero.facs, and hence fac1, is NULL, no heterogeneous variance models are fitted. 
#' * If hetero.facs is not NULL, it will be used to set fac1 and fac2, ignoring 
#'   whether fac1 or fac2 has been set.
#' * The starting model in the supplied asrt should be a homogeneous variance model.
#' 
#' * There are two selections methods: forward or backward.
#'   + Backward selection has only been implemented for the the Residual term. 
#'   + Backward selection begins by fitting the maximal heterogeneous variance model 
#'     as a starting model and then investigates whether the model can be simplified.
#'   + If hetero.facs is a multicomponent list, then it is likely that forward selection 
#'     is used so that each variance model is compared with the previously fitted variance 
#'     model. 
#'
#' * term specifies a term or the Residual, that is to be tested for heterogeneity.
#' * If the term is Residual or contains Seed, then it is designated as a Seed term and is assumed to 
#'    be the Residual term in the analysis.
#'   + Residual terms are fitted using the dsum function and there needs to be factors in the data 
#'    for (i) fac1, (ii) fac2, and (iii) paste0(fac1,fac2).
#' * If term is neither Residual nor contains Seed, it will assumed to be an Indv term for the individuals 
#'    in the experiment, unless it includes Spike as a factor so that it will be designated as a Spike term.
#' 
#' * For any term that is specified using a set of factors, the factors in the term need to include 
#'   fac1 followed by fac2 and separated by a colon. 
#'   + also, the data must be ordered for fac1 then fac2.
#'   
#' * If any variances are singular, then they are fixed at 0.0001 and a warning issued.
#'    
choose1or2WayHetVarModel <- function(asrtests.obj, hetero.term, selection = "forward", 
                                     fac1=NULL, fac2 = NULL, hetero.facs = NULL)
{
  #If hetero.facs is not NULL, use it to set fac1 and fac2
  if (all(!is.null(hetero.facs)))
  {
    fac1 <- hetero.facs[1]
    if (length(hetero.facs) > 1)
      fac2 <- hetero.facs[2]
    else
      fac2 <- NULL
  }
  #pass through the function if fac1 is NULL 
  if (!is.null(fac1))
  { 
    if (!(selection) %in% c("forward", "backward"))
      stop("selection must be set to either forward or backward")
    
    #Check that required factors are in the asreml.obj data.frame
    dat <- asrtests.obj$asreml.obj$call$data; if (is.symbol(dat)) dat <- eval(dat)
    dat.cols <- names(dat)
    het.facs <- c(fac1,fac2)
    if (!all(het.facs %in% dat.cols))
      stop("Some factors that are to be investigaed for variance heterogeneity model are not in data: ", 
           paste0(het.facs[!(het.facs %in% dat.cols)], collapse =", "))
    
    #Begin fitting the models
    asreml.options(keep.order=TRUE, fail = "soft", ai.sing = TRUE, extra = 5)
    term.unit <- "Indv"
    if (grepl("Seed", hetero.term) || hetero.term == "Residual") 
      term.unit <- "Residual"
    else
    { 
      if (grepl("Spike", hetero.term)) term.unit<- "Spike"
    }
    
    test.name <- paste("Try heterogeneity for",fac1, "for", term.unit)
    if (term.unit != "Residual")
    {
      if (selection != "forward")
        stop("backward selection has not been implement for ", term.unit)
      old.term <- hetero.term
      new.term  <- gsub(fac1, paste0("idh(", fac1,")"), old.term)
      asrtests.obj <- changeModelOnIC(asrtests.obj,
                                      dropRandom = old.term, 
                                      addRandom = new.term, 
                                      label=test.name, checkboundaryonly = TRUE,
                                      denDF = "numeric", ssType = "conditional", 
                                      update = FALSE)
      if (!is.null(fac2)) #two-way
      {
        if (getTestEntry(asrtests.obj, label = test.name)$action =="Swapped") #both fac1 and fac2
        {
          old.term <- new.term
          test.name <- paste("Try heterogeneity in", fac1, "&", fac2, "for", term.unit)
          inter <- paste(fac1, fac2, sep = ":")
          new.term <- gsub(inter, paste0("idh(", gsub(":", "", inter), ")"), hetero.term)
        } else #fac2 only
        {
          test.name <- paste("Try heterogeneity in", fac2, "for", term.unit)
          new.term <- gsub(fac2, paste0("idh(", fac2, ")"), old.term)
        }
        asrtests.obj <- changeModelOnIC(asrtests.obj,
                                        dropRandom = old.term, 
                                        addRandom = new.term, 
                                        label=test.name, checkboundaryonly = TRUE,
                                        denDF = "numeric", ssType = "conditional", 
                                        update = FALSE)
      }
    } else #A Residual hetero.term
    {
      #Check for ancillary factors needed to fit heterogeneity to the Residual and add if not there 
      
      #Add any required Wfacs not in dat
      if (length(hetero.facs)  == 1)
      { 
        if (!(hetero.facs[1] %in% names(dat)))
          stop(paste0("Cannot find the factor ", hetero.facs[1]))
        dat <- dat[order(unlist(dat[hetero.facs])), ]
      } else
      { 
        if (length(hetero.facs) == 2)
        { 
          inter <- paste0(hetero.facs, collapse = "")
          if (!(inter %in% names(dat.cols)))
            dat[inter] <- fac.combine(dat[hetero.facs], combine.levels=T, sep="")
          dat <- dat[order(unlist(dat[inter])), ]
        } else
          stop("hetero.facs is not of length 1 or 2")
      }
      
      if (asrtests.obj$asreml.obj$formulae$residual == ~NULL)
      {
        asr <- asrtests.obj$asreml.obj
        asr <- newfit(asr, residual. = ~ units)
        asrtests.obj$asreml.obj <- asr
      }
      #Fit and compare heterogeneous models for a Residual
      if (selection == "forward")
      {    
        test.name <- paste("Try heterogeneity for",fac1, "for", term.unit)
        new.term  <- paste0("dsum(~ units | ", fac1, ")")
        asrtests.obj <- changeModelOnIC(asrtests.obj,
                                        newResidual = new.term, 
                                        label=test.name, checkboundaryonly = TRUE,
                                        denDF = "numeric", ssType = "conditional", 
                                        update = FALSE, 
                                        data = dat) #pass dat via data, not by overwriting in the call 
        
        if (!is.null(fac2)) #two-way
        {
          if (getTestEntry(asrtests.obj, label = test.name)$action =="Swapped") 
          {
            test.name <- paste("Try heterogeneity in", fac1, "&", fac2, "for", term.unit)
            inter <- paste0(fac1, fac2)
            new.term  <- paste0("dsum(~ units | ", inter, ")")
          } else 
          {
            # #Order data for fac1 and replace data in the call
            dat <- dat[order(unlist(dat[fac2])), ]
            #Set label and new residual
            test.name <- paste("Try heterogeneity in", fac2, "for", term.unit)
            new.term  <- paste0("dsum(~ units | ", fac2, ")")
          }
          asrtests.obj <- changeModelOnIC(asrtests.obj,
                                          newResidual = new.term, label=test.name, 
                                          checkboundaryonly = TRUE,
                                          denDF = "numeric", ssType = "conditional", 
                                          update = TRUE, 
                                          data = dat) #pass dat via data, not by overwriting in the call
        }
      } else #backward elimination
      {
        if (is.null(fac2)) #one-way heterogeneity
        {
          #Fit one-way heterogeneity
          lab <- paste("Start with residual heterogeneity for",fac1, 
                       "for", term.unit)
          new.term  <- paste0("dsum(~ units | ", fac1, ")")
          asrtests.obj <- changeTerms(asrtests.obj,
                                      newResidual = new.term, 
                                      label=lab, checkboundaryonly = TRUE,
                                      denDF = "numeric", ssType = "conditional", 
                                      update = FALSE, 
                                      data = dat) #pass dat via data, not by overwriting in the call
          #Test for one-way heterogeneity
          test.name <- paste("Try homogeneity for",fac1, "for", term.unit)
          new.term  <- "- (.)"
          asrtests.obj <- changeModelOnIC(asrtests.obj,
                                          newResidual = new.term, 
                                          label=test.name, checkboundaryonly = TRUE,
                                          denDF = "numeric", ssType = "conditional", 
                                          update = FALSE, 
                                          data = dat) #pass dat via data, not by overwriting in the call
        } else #two-way heterogeneity
        {
          #Fit two-way heterogeneity
          lab <- paste("Start with residual heterogeneity in", 
                       fac1, "&", fac2, "for", term.unit)
          new.term  <- paste0("dsum(~ units | ", inter, ")")
          asrtests.obj <- changeTerms(asrtests.obj,
                                      newResidual = new.term, 
                                      label=lab, checkboundaryonly = TRUE,
                                      denDF = "numeric", ssType = "conditional", 
                                      update = FALSE, 
                                      data = dat) #pass dat via data, not by overwriting in the call
          
          #Test for simplification of the two-way model
          #Order data for fac1 and replace data in the call
          tmp <- dat[order(unlist(dat[fac1])), ]
          #Test new Residual hetero.term
          test.name <- paste("Try homogeneity in",fac2, "for", term.unit)
          new.term  <- paste0("dsum(~ units | ", fac1, ")")
          asrtests.obj <- changeModelOnIC(asrtests.obj,
                                          newResidual = new.term, 
                                          label=test.name, checkboundaryonly = TRUE,
                                          denDF = "numeric", ssType = "conditional", 
                                          update = TRUE, 
                                          data = tmp) #pass tmp via data, not by overwriting in the call
          if (getTestEntry(asrtests.obj, label = test.name)$action !="Swapped") 
          {
            #Order data for fac2 and replace data in the call
            tmp <- dat[order(unlist(dat[fac2])), ]
            test.name <- paste("Try homogeneity in", fac1, "for", term.unit)
            new.term  <- paste0("dsum(~ units | ", fac2, ")")
          } else 
          {
            test.name <- paste("Try homogeneity in", fac1, "&", fac2, "for", term.unit)
            new.term  <- "- (.)" 
            tmp <- dat
          }
          asrtests.obj <- changeModelOnIC(asrtests.obj,
                                          newResidual = new.term, label=test.name, 
                                          checkboundaryonly = TRUE,
                                          denDF = "numeric", ssType = "conditional", 
                                          update = TRUE, 
                                          data = tmp) #pass dat via data, not by overwriting in the call
        }
      }
    }
    vpar.S <- asrtests.obj$asreml.obj$vparameters.con[asrtests.obj$asreml.obj$vparameters.con == "S"]
    if (length(vpar.S) > 0)
    {
      vpt <- current.asrt$asreml.obj$vparameters.type[names(vpar.S)]
      vpar.s <- vpar.S[vpt %in% c("V","G")]
      if (length(vpar.S) > 0)
      {
        asr <- setvarianceterms(current.asrt$asreml.obj$call, 
                                terms = names(vpar.s), initial.values = fix.bound.val,
                                bounds = "F", ignore.suffices = FALSE)
        asrtests.obj$asreml.obj <- asr
        warning("One or more singular variances have been fixed at 0.0001")
      }
    }
  }
  
  return(asrtests.obj)  
}

#'## Function to determine the fitted var.model in an asreml.obj
getVarModel <- function(asreml.obj, hetero.facs, idh = FALSE)
{  
  res.comp <- names(asreml.obj$vparameters)
  if (idh)
  {
    res.facs <- unlist(lapply(hetero.facs, 
                              function(het.fac, comp) 
                              {
                                if (!any(grepl(paste(het.fac, het.fac, sep = "!"), res.comp))) 
                                  het.fac <- NULL
                                return(het.fac)
                              }, res.comp))
  } else
  {
    res.comp <- res.comp[grepl("(!R$)", res.comp)]
    res.facs <- unlist(lapply(hetero.facs, 
                              function(het.fac, comp) 
                              {
                                if (!all(grepl(het.fac, res.comp))) het.fac <- NULL
                                return(het.fac)
                              }, res.comp))
  } 
  var.model <- ifelse(is.null(res.facs),
                      "Homogeneous variance",
                      paste("Variance differs by", paste(res.facs, collapse = " and ")))
  
  return(var.model)
}

#'### Function to fix any variance components that are bound to `fix.bound.val`
fixBoundVCs <- function(asrtests.obj, fix.bound.val = 1e-04)
{
  asr <- asrtests.obj$asreml.obj
  asr <- newfit(asr)
  ICb4 <- infoCriteria(asr)
  vcomp <- (summary(asr)$varcomp)
  vcomp <- vcomp[vcomp$bound == "B",]
  iter <- 0
  repeat
  {
    #repeat until 10 iteration or no bound components
    if (iter == 10 || nrow(vcomp) == 0) break
    #Fix all currently known bound components
    terms <- rownames(vcomp)
    call <- asr$call
    asr <- setvarianceterms(call, terms = terms, 
                            ignore.suffices = FALSE, 
                            initial.values = fix.bound.val, bounds = "F")
    #Check no bound terms
    vcomp <- (summary(asr)$varcomp)
    vcomp <- vcomp[vcomp$bound == "B",]
    iter <- iter + 1
  }
  #Recreate asrtests object and add a line to the test summary
  ICnow <- infoCriteria(asr)
  asrtests.obj <- with(asrtests.obj, as.asrtests(asr, wald.tab = wald.tab, 
                                                 test.summary = test.summary, 
                                                 label = "Fix bound variance components"))
  kline <- nrow(asrtests.obj$test.summary)
  asrtests.obj$test.summary <- within(asrtests.obj$test.summary, 
                                      {
                                        DF[kline] <- 0
                                        denDF[kline] <- ICb4$varDF - ICnow$varDF
                                        AIC[kline] <- ICb4$AIC - ICnow$AIC
                                        BIC[kline] <- ICb4$BIC - ICnow$BIC
                                        action[kline] <- "Swapped"
                                      })
  return(asrtests.obj)
}

#'## Function to extract the random effects for terms or residuals and some fitted values for diagnostic checking 
#'
#'* The arguments `asreml.obj`, `terms`, `use` and `sep` are as described for 
#'  convEffectNames2DataFrame.asreml in the manual.
#'* Note that if `terms` is NULL, then the residuals are used to form res.raw and res.dat; 
#'  otherwise, they will be the effects derived from the random effects for `terms` and 
#'  should occur as names of the elements of the `vparameters` component of the `asreml.obj`; 
#'  multiple terms based on the same factors can be provided.
#'* The argument `classify.facs` is used in manipulating the EMMs to acts as fitted values; when term is 
#'  NULL, the fitted values will be used.
#'* The argument `preds` is a predictions.frame that contains the EMMs to be used as 
#'  fitted values and is only used if `classify.facs` is not NULL.
#'* The data data.frame will be extracted from the `call` of the supplied `asreml.obj` and will be 
#'  used to form a data.frame to which the residuals and fitted/EMM values are to be added: 
#'    + It must contain the factors and variables in terms and `classify.facs` and, for the 
#'      factors and variables in `terms`, it must have at least the combinations of them that 
#'     occur in the data.
#'* The `tol` argument deletes the rows of the res.dat frame whose res.raw values are lessn than tol. 
getDiagnosticEffFit <- function(asreml.obj, terms = NULL, use = "G.aom", classify.facs = NULL, 
                                preds, tol = .Machine$double.eps ^ 0.5, sep =  ":")
{
  dat  <- asreml.obj$call$data; if (is.symbol(dat)) dat <- eval(dat)
  if (all(is.null(terms)))
    terms.facs <- NULL
  else
    terms.facs <-  stringr::str_split_1(terms[1], ":")
  if (!all(is.null(classify.facs)) || !all(is.null(terms)))
    meta.dat <- unique(dat[union(terms.facs, classify.facs)])
  else 
    meta.dat <- NULL 
  
  if (all(is.null(terms)))
  {
    res.dat <- within(dat, 
                      {
                        res.raw <- residuals(asreml.obj, type="response")
                        res.std <- residuals(asreml.obj, type="stdCond")
                      })
  } else
  { 
    #Get the standardized terms effects
    G.dat <- lapply(terms, 
                    function(term, asreml.obj, use, sep)
                      t <- convEffectNames2DataFrame.asreml(asreml.obj, term = term, 
                                                            use = use, sep = sep),
                    asreml.obj = asreml.obj, use = use, sep = sep)
    G.dat <- do.call(rbind, G.dat)
    if (is.null(G.dat))
    {
      warning("Could not extract random effects for ", terms, "; using residuals instead")
      res.dat <- within(dat, 
                        {
                          res.raw <- residuals(asreml.obj, type="response")
                          res.std <- residuals(asreml.obj, type="stdCond")
                        })
    } else
    {
      effnam <- attr(G.dat, which = "effect.names")
      G.std <- asreml.obj$aom$G[effnam, "stdCondBlup"]
      G.dat$res.std <- G.std
      
      #Get the raw G effects
      ran.LP <- asreml.obj$coefficients$random[effnam, "effect"]
      G.dat$res.raw <- ran.LP
      
      #Form res.dat
      res.dat <- dplyr::left_join(G.dat, meta.dat)
    }
  } 
  
  #Add the predicted values for classify.facs or fitted values to res.dat
  if (is.null(classify.facs))
  { 
    res.dat$fit <- fitted(asreml.obj)
    res.dat <- res.dat[!is.na(res.dat$res.raw), ]
  } else
  {
    fit.dat <- preds[c(classify.facs,"predicted.value")]
    fit.dat <- dplyr::left_join(meta.dat, fit.dat)
    
    names(fit.dat)[match("predicted.value", names(fit.dat))] <- "fit"
    res.dat <- dplyr::left_join(res.dat, fit.dat)
    res.dat <- res.dat[!is.na(res.dat$res.raw), ]
  } 
  res.dat <- res.dat[abs(res.dat$res.raw) > tol,]
  return(res.dat)
}

#'### Function to plot a set of residuals
plotResiduals <- function(res.dat, x = "fit", 
                          x.facet = NULL,    #a character naming one or more columns of res.dat
                          xdodge = FALSE, dodge.width = 0.5, 
                          xtit = "Predicted values", resentity = "individual", 
                          response, colour.col = NULL, #each a character
                          var.mod = NULL,              #a character
                          yquantile = NULL,            #a numeric between 0 and 1
                          outlier.label = NULL,        #a character naming a column of res.dat
                          which.plots = c("qqplot", "raw.residuals", "std.residuals"), 
                          ggplotFuncs = NULL)
{
  #Check which.plots
  if (!any(which.plots %in% c("qqplot", "raw.residuals", "std.residuals")))
    stop("which.plots is not one or more of qqplot, raw.residuals and std.residuals")
  
  plts <- list()
  if (any(c("raw.residuals", "std.residuals") %in% which.plots))
  {  
    #Set common theme elements
    plotFuncs <- list(
      xlab(xtit), #no units as incorrect for transformed response
      geom_hline(yintercept=0, linetype = "solid", alpha=0.3),
      theme_bw(),
      theme(legend.position = "none"))#,
    #Add the facet grid
    if (!is.null(x.facet)) 
      plotFuncs <- c(plotFuncs, list(facet_grid(cols = vars(.data[[x.facet]]), 
                                                labeller = label_both)))
    
    if ("raw.residuals" %in% which.plots)
    {   
      #Create raw residuals plot
      yvar <- "res.raw"
      ytit <- paste("Raw", resentity, "residuals")
      plts[["resraw.plt"]] <- 
        ggplot(data=res.dat, aes(y=res.raw, x=.data[[x]], colour=.data[[colour.col]])) +
        geom_point(na.rm = TRUE, alpha=0.75, size=1.0) + 
        ylab(ytit) +
        ggtitle(paste(response,ytit,sep=" - ")) +
        geom_vline(xintercept=0, linetype = "dashed", alpha=0.4) +
        geom_abline(slope=-1, linetype = "dashed", alpha=0.3) +
        plotFuncs
    }
    
    if ("std.residuals" %in% which.plots)
    {
      #Set up for standardized residuals plot
      yvar <- "res.std"
      ytit <- paste("Standardized", resentity, "residuals")
      res.lim <- c(3,-3)
      if (!is.null(outlier.label))
        outliers.dat <- subset(res.dat, abs(res.std) > 3)
      
      #Start the res.std plot with the colour.col specified, if required 
      if (is.null(colour.col))
        resstd.plot <- ggplot(data=res.dat, 
                              aes(y=res.std, x=.data[[x]])) 
      else
        resstd.plot <- ggplot(data=res.dat, 
                              aes(y=res.std, x=.data[[x]], colour=.data[[colour.col]]))
      
      #Add the point plot, with dodging if required
      if (xdodge)
        resstd.plot <- resstd.plot + 
        geom_point(na.rm = TRUE, alpha=0.75, size=1, position = position_dodge2(width = dodge.width))
      else
        resstd.plot <- resstd.plot + 
        geom_point(na.rm = TRUE, alpha=0.75, size=1)
      
      #Add quantile when plotting other than the predicted valus
      if (!is.null(yquantile) && x != "Predicted.value")
      {
        resquant.dat <- summarise(group_by(res.dat, .data[[x]], .data[[x.facet]]),
                                  res.quant.hi = quantile(res.std, 1 - yquantile, na.rm = TRUE),
                                  res.quant.lo = quantile(res.std, yquantile, na.rm = TRUE))
        resstd.plot <- resstd.plot +
          geom_point(data = resquant.dat, aes(y = res.quant.hi), colour = "black", size = 2) +
          geom_point(data = resquant.dat, aes(y = res.quant.lo), colour = "black", size = 2)
      }
      
      resstd.plot <- resstd.plot + 
        ylab(ytit) +
        ggtitle(paste(response, "- Std Residuals vs",xtit),
                subtitle = var.mod) +
        geom_hline(yintercept=res.lim, colour = "maroon") +
        (if (!is.null(outlier.label))         
          geom_text(data = outliers.dat,
                    aes(label = .data[[outlier.label]],
                        y = res.std - 0.2*sign(res.std)),
                    colour="black", size=3)
        ) +
        scale_y_continuous(breaks = seq(min(-3, floor(min(res.dat$res.std, na.rm = TRUE))),
                                        max(3, ceiling(max(res.dat$res.std, na.rm = TRUE))),
                                        by = 2)) +
        plotFuncs
      
      if (!is.null(ggplotFuncs))
      {
        for (f in ggplotFuncs)
          resstd.plot <- resstd.plot + f
      }
      plts[["resstd.plt"]] <- resstd.plot
    }
  }
  
  #Do the qqplot if required
  if ("qqplot" %in% which.plots)
  { 
    plts[["qqplt"]] <- 
      ggplot(data=res.dat, mapping = aes(sample = res.std)) +
      stat_qq_point(alpha=1, size=1) +
      stat_qq_line(alpha=0.4, colour="blue") + stat_qq_band(bandType = "ts") +
      labs(x = "Theoretical Quantiles", y = "Sample Quantiles",
           title = paste(response,"- Normal probability plot")) +
      geom_hline(yintercept=0, linetype = "dashed", alpha=0.3) +
      geom_vline(xintercept=0, linetype = "dashed", alpha=0.3) +
      theme_bw()
  }
  return(plts) #at list haviing one or more of resraw.plt, resstd.plt and qqplt
}

#'## Function to make a table of the fitted Residual variances
#'
#' * fac.hetero = NULL results in only a list of the `var.mod` values to be produced
#' * res.idh indicates the function used to specify heterogeneous variances
#' 
convertRownames2Column <- function(x, name = "Trait")
{
  x <- cbind(rownames(x),as.data.frame(x))
  names(x)[1] <- name
  rownames(x) <- NULL
  return(x)
}

tableResVars <- function(models, dat = indv.dat, 
                         fac.hetero = NULL,      #factor that has heterogeneous variances
                         Control.lev = "Control",#the level to which other levels are compared 
                         res.idh = TRUE,         #whether idh or dsum was used to obtain hetero variances
                         trait.types = NULL,     #stubs for each of the different imaging trait types.
                         trait.digits = 3,       #the number of significant digits for each trait type.
                         printTable = TRUE)       
{
  if (is.null(fac.hetero))
  {
    var.summ <- do.call(rbind, lapply(models, function(mod) mod$var.mod))
    var.summ <- gsub("Variance differs by", "Heterogeneous", var.summ)
    colnames(var.summ) <- "Variance model"
    var.summ <- convertRownames2Column(var.summ, name = "Trait")
    print(kable(var.summ, row.names = FALSE, align = "l"))
  } else
  {
    levs <- levels(dat[[fac.hetero]])
    nlevs <- length(levs)
    var.summ <- lapply(models, function(mod, fac.hetero, levs, res.idh) 
    {
      var.mod <- mod[["var.mod"]]
      var.tbl <- summary(mod$asreml.obj)$varcomp
      nlevs <- length(levs)
      if (grepl("Homogeneous", var.mod))
      {
        var.tbl <- var.tbl[grepl("\\!R", rownames(var.tbl)), ]
        if (nrow(var.tbl) > 1)
          var.tbl <- var.tbl[grepl("!R", rownames(var.tbl)), ]
        var.dat <- rep(var.tbl$component, nlevs)
      } else
      {
        if (!res.idh) #dsum was used
        {
          var.dat <- sapply(levs, 
                            function(lev, var.tbl) 
                              var.tbl[grepl(paste0(lev,"\\!R"), rownames(var.tbl)), "component"],
                            var.tbl = var.tbl)
        }
        else #idh was used
          var.dat <- sapply(levs, 
                            function(lev, var.tbl) 
                              var.tbl[grepl(paste0(lev,"$"), rownames(var.tbl)), "component"],
                            var.tbl = var.tbl)
        var.mod <- paste("Heterogeneous", fac.hetero)
      }
      var.dat <- as.data.frame(cbind(var.mod, as.data.frame(rbind(var.dat))))
      names(var.dat) <- c("Variance model", levs)
      rownames(var.dat) <- NULL
      return(var.dat)
    }, fac.hetero = fac.hetero, levs = levs, res.idh = res.idh)
    var.summ <- do.call(rbind, var.summ)
    var.summ <- var.summ[order(rownames(var.summ)), ]
    
    var.summ[paste(substr(levs[1],1,2), substr(levs[2:nlevs],1,2), sep = ":")] <- 
      as.data.frame(lapply(var.summ[levs[2:nlevs]], 
                           function(col, col1) {col <- col1/col}, col1 = var.summ[levs[1]]))
    var.summ[Control.lev] <- signif(var.summ[Control.lev],6)
    var.summ <- convertRownames2Column(var.summ)
    
    #Output the table
    if (printTable)
    { 
      if (all(is.null(trait.types)))
      {
        print(kable(var.summ, row.names = FALSE, 
                    align = paste0(c("l", "l", rep("r", nlevs*2 - 1)), collapse = ""), 
                    digits = c(0,0,rep(trait.digits,nlevs),2), 
                    caption = "Residual variances"))
      } else
        for (trait in trait.types)
        { 
          names(trait.digits) <- trait.types
          t.summ <- var.summ[grepl(paste0(trait,".[0-9]"), var.summ$Trait), ]
          print(kable(t.summ, row.names = FALSE, 
                      align = paste0(c("l", "l", rep("r", nlevs*2 - 1)), collapse = ""), 
                      digits = c(0,0,rep(trait.digits[trait],nlevs),2), 
                      caption = paste("Residual variances for", trait)))
        }
    }
}
  invisible(var.summ)
}

#Function to tabulate chosen models
tableChosenModels <- function(models, responses, respnames = responses, printTable = TRUE)
{
  tbl <- do.call(rbind, lapply(responses,
                               function(kresp, models) 
                               {  
                                 mod <- paste(unlist(models[[kresp]]$chosen.terms), 
                                              collapse=" + ")
                                 if (mod == "") mod <- "Overall mean"
                                 return (mod)
                               }, models))
  tbl <- data.frame(Trait=respnames, `Chosen model`=tbl)
  if (printTable)
    print(tbl, right = FALSE)
  invisible(tbl)
}

#Function to table Variance Models
tableVarModels <- function(models, responses, respnames = responses, printTable = TRUE)
{ 
  tbl <- do.call(rbind, lapply(models, function(mod) mod$var.mod)) 
  tbl <- data.frame(Trait=respnames, `Variance model`=tbl)
  rownames(tbl) <- NULL
  if (printTable)
    print(tbl, right = FALSE)
  invisible(tbl)
}

#'## Function to compute the R2.adj value for a model and to add it as a composnent of the model list
addR2adj2model <- function(model, which.fixed = ~ ., which.random = ~ .)
{
  set.daeTolerance(1e-04, 1e-04)
  model[["R2.adj"]] <- R2adj(model$asreml.obj, 
                             include.which.fixed = which.fixed, 
                             include.which.random = which.random,
                             orthogonalize = "eigen")
  return(model)
}

#'## Function to make a table of AICs and R2(adj) values for the model fits stored in a models object
#'
#'* It uses the asreml.obj for a response to calculate the AIC for the fitted model and, if 
#'  R2.adj is not stored in the models, the R2(adj) value for the fixed terms specified by 
#'  `which.fixed` and all the random terms fitted.
tableAICR2ASRemlobjs <- function(models, which.fixed = ~ ., respnames = NULL)
{
  set.daeTolerance(1e-04, 1e-04)
  AICs <- lapply(models, 
                 function(mod) 
                 {
                   AIC <- infoCriteria(mod$asreml.obj, IClikelihood = "full")$AIC
                   R2.adj <- mod$R2.adj
                   if (is.null(R2.adj))
                     R2.adj =  R2adj(mod$asreml.obj, 
                                     include.which.fixed = which.fixed, 
                                     include.which.random = ~ .,
                                     orthogonalize = "eigen")
                   AIC <- data.frame(AIC = AIC, 
                                     R2.adj = R2.adj)
                   return(AIC)
                 })
  AICs <- do.call(rbind, AICs)
  if (!is.null(respnames)) rownames(AICs) <- respnames
  ncols <- ncol(AICs)
  AICs <- convertRownames2Column(AICs)
  print(kable(AICs, row.names = FALSE, align = c("l", rep("r",ncols)), digits = 1))
  invisible(AICs)
}

#'## Function to make a table of AICs used in comparing several local spatial models and R2(adj) values for the selected models
#'
#'* It uses the spatial.IC components formed in fitting several local spatial models 
#'  fitted to each of the responses to calculate the AICs and an asreml.obj to calculate 
#'  the R2(adj) values for large-scale and local spatial terms.
tableAICR2SpatialICs <- function(models, which.fixed = ~ .)
{
  set.daeTolerance(1e-04, 1e-04)
  AICs <- lapply(models, 
                 function(mod) 
                 {
                   AIC <- mod$spatial.IC[,"AIC"]
                   names(AIC) <- rownames(mod$spatial.IC)
                   min.mod <- mod$best.spatial.mod
                   R2.adj =  R2adj(mod$asreml.obj, 
                                   include.which.fixed = which.fixed, 
                                   orthogonalize = "eigen", 
                                   include.which.random = ~ .)
                   AIC <- data.frame(AIC.nonspatial = AIC[1], 
                                     AIC.fitted = mod$best.spatial.AIC, 
                                     AIC.diff = AIC[min.mod]-AIC[1], 
                                     model.fitted = min.mod,
                                     R2adj.fitted = R2.adj,
                                     AIC.TPPSC2 = AIC["TPPSC2"])
                   return(AIC)
                 })
  AICs <- do.call(rbind, AICs)
  ncols <- ncol(AICs)
  AICs <- convertRownames2Column(AICs)
  print(kable(AICs, row.names = FALSE, align = c("l", rep("r",ncols)), digits = 1))
  invisible(AICs)
}

#'## Function to make a table of EMM statistics
tableEMMstats <- function(models, classify.facs)
{  
  EMM.medns <- lapply(models, 
                      function(mod)
                      {
                        EMMs.fit <- mod$diffs$predictions
                        names(EMMs.fit)[match("predicted.value", names(EMMs.fit))] <- "EMMs.fit"
                        EMMs.ns <- mod$nonspatial.preds
                        names(EMMs.ns)[match("predicted.value", names(EMMs.ns))] <- "EMMs.ns"
                        EMMs.dat <- suppressMessages(
                          dplyr::inner_join(EMMs.fit[c(classify.facs,"EMMs.fit")],
                                            EMMs.ns[c(classify.facs,"EMMs.ns")]))
                        medians <- data.frame(medn.ns = median(EMMs.dat$EMMs.ns, na.rm = TRUE), 
                                              medn.fit = median(EMMs.dat$EMMs.fit, na.rm = TRUE), 
                                              medn.diffs = median(EMMs.dat$EMMs.fit - 
                                                                    EMMs.dat$EMMs.ns, 
                                                                  na.rm = TRUE),
                                              r2.EMMs = cor(EMMs.dat$EMMs.fit, 
                                                            EMMs.dat$EMMs.ns)^2,
                                              MAD = median(abs(EMMs.dat$EMMs.fit - 
                                                                 EMMs.dat$EMMs.ns), 
                                                           na.rm = TRUE))
                        return(medians)
                      })
  EMM.medns <- do.call(rbind, EMM.medns)
  EMM.medns <- within(EMM.medns, 
                      {
                        fit.model <- AICs$model.fitted
                        MAD100 = MAD/medn.ns*100
                      })
  ncols <- ncol(EMM.medns)
  EMM.medns <- convertRownames2Column(EMM.medns)
  print(kable(EMM.medns, row.names = FALSE, align = c("l", rep("r",ncols)), 
              digits = 3, caption = "Medians calculated from EMM values"))
  
  invisible(EMM.medns)
}

#'## Function to make a table of SE statistics
tableSEstats <- function(models, classify.facs)
{  
  se.medns <- lapply(models, 
                     function(mod)
                     {
                       se.fit <- mod$diffs$predictions
                       names(se.fit)[match("standard.error", names(se.fit))] <- "se.fit"
                       se.ns <- mod$nonspatial.preds
                       names(se.ns)[match("standard.error", names(se.ns))] <- "se.ns"
                       se.dat <- suppressMessages(
                         dplyr::inner_join(se.fit[c(classify.facs,"se.fit")],
                                           se.ns[c(classify.facs,"se.ns")]))
                       
                       medians <- data.frame(medn.ns = median(se.dat$se.ns, na.rm = TRUE), 
                                             medn.fit = median(se.dat$se.fit, na.rm = TRUE), 
                                             medn.diffs = median(se.dat$se.fit - 
                                                                   se.dat$se.ns, 
                                                                 na.rm = TRUE))
                       return(medians)
                     })
  se.medns <- do.call(rbind, se.medns)
  se.medns <- within(se.medns, 
                     {
                       fit.model <- AICs$model.fitted
                       medn.diffs100 = medn.diffs/medn.ns*100
                     })
  ncols <- ncol(se.medns)
  se.medns <- convertRownames2Column(se.medns)
  print(kable(se.medns, row.names = FALSE, align = c("l", rep("r",ncols)), 
              digits = 3, caption = "Medians calculated from standard errors"))
  
  invisible(se.medns)
}

