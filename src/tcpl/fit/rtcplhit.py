import pandas as pd
import rpy2

tcplhit2_in_r = """
library(tcplfit2)

tcplhit2_core_is <- function(params, conc, resp, cutoff, onesd, bmed = 0, conthits = T, 
                           aicc = F, identifiers = NULL,bmr_magic=1.349) {
  # initialize parameters to NA
  a <- b <- tp <- p <- q <- ga <- la <- er <- top <- ac50 <- ac50_loss <- ac5 <- ac10 <- ac20 <- acc <- ac1sd <- bmd <- NA_real_
  bmdl <- bmdu <- caikwt <- mll <- NA_real_
  # get aics and degrees of freedom
  aics <- sapply(params$modelnames, function(x) {
    params[[x]][["aic"]]
  })
  dfs <- sapply(params$modelnames, function(x) {
    if(x == "cnst") 1L
    else length(params[[x]][gsub("_sd","",names(params[[x]])[grepl("_sd",names(params[[x]]))])])
  })
  aics <- aics[!is.na(aics)]
  if (sum(!is.na(aics)) == 0) {
    # if all fits failed, use none for method
    fit_method <- "none"
    rmse <- NA_real_
  } else {
    # use nested chisq to choose between poly1 and poly2, remove poly2 if it fails.
    # pvalue hardcoded to .05
    aics <- nestselect(aics, "poly1", "poly2", dfdiff = 1, pval = .05)
    dfs <- dfs[names(dfs) %in% names(aics)]
    # it's useful to keep original aics so we can extract loglikelihoods for nested models (above) and mll (below)
    if (aicc) saics <- aics + 2 * dfs * (dfs + 1) / (length(resp) - dfs - 1) else saics <- aics
    if (conthits) {
      # get aikaike weight of winner (vs constant) for cont hitcalls
      # never choose constant as winner for cont hitcalls
      nocnstaics <- saics[names(saics) != "cnst"]
      fit_method <- names(nocnstaics)[which.min(nocnstaics)]
      caikwt <- exp(-saics["cnst"] / 2) / (exp(-saics["cnst"] / 2) + exp(-saics[fit_method] / 2))
      if (is.nan(caikwt)) {
        term <- exp(saics["cnst"] / 2 - saics[fit_method] / 2)
        if (term == Inf) {
          caikwt <- 0
        } else {
          caikwt <- 1 / (1 + term)
        }
        # caikwt <- 1
      }
    } else {
      fit_method <- names(saics)[which.min(saics)]
    }
    fitout <- params[[fit_method]]
    rmse <- fitout$rme
    # hacky way to get modpars without hardcoding the names or needing the list
    # basically each model (except cnst) has an sd for each parameter
    # we can use that to select all model parameters
    modpars <- fitout[gsub("_sd","",names(fitout)[grepl("_sd",names(fitout))])]
    if(fit_method == "cnst") modpars <- fitout["er"]  #since no sd for cnst
    list2env(fitout, envir = environment()) # put all parameters in environment
  }

  n_gt_cutoff <- sum(abs(resp) > cutoff)

  # compute discrete or continuous hitcalls
  if (fit_method == "none") {
    hitcall <- 0
  } else if (conthits) {
    mll <- length(modpars) - aics[[fit_method]] / 2
    hitcall <- hitcontinner(conc, resp, top, cutoff, er,
      ps = unlist(modpars), fit_method,
      caikwt = caikwt, mll = mll
    )
  } else {
    hitcall <- hitloginner(conc, resp, top, cutoff, ac50)
  }

  bmr <- onesd * bmr_magic 
  if (hitcall > 0) {

    # fill ac's; can put after hit logic
    ac5 <- acy(.05 * top, modpars, type = fit_method) # note: cnst model automatically returns NAs
    ac10 <- acy(.1 * top, modpars, type = fit_method)
    ac20 <- acy(.2 * top, modpars, type = fit_method)
    acc <- acy(sign(top) * cutoff, modpars, type = fit_method)
    ac1sd <- acy(sign(top) * onesd, modpars, type = fit_method)
    bmd <- acy(sign(top) * bmr, modpars, type = fit_method)

    # get bmdl and bmdu
    bmdl <- bmdbounds(fit_method,
      bmr = sign(top) * bmr, pars = unlist(modpars), conc, resp, onesidedp = .05,
      bmd = bmd, which.bound = "lower"
    )
    bmdu <- bmdbounds(fit_method,
      bmr = sign(top) * bmr, pars = unlist(modpars), conc, resp, onesidedp = .05,
      bmd = bmd, which.bound = "upper"
    )
  }

  top_over_cutoff <- abs(top) / cutoff
  conc <- paste(conc, collapse = "|")
  resp <- paste(resp, collapse = "|")

  # row contains the specified columns and any identifying, unused columns in the input
  name.list <- c(
    "n_gt_cutoff", "cutoff", "fit_method",
    "top_over_cutoff", "rmse", "a", "b", "tp", "p", "q", "ga", "la", "er", "bmr", "bmdl", "bmdu", "caikwt",
    "mll", "hitcall", "ac50", "ac50_loss", "top", "ac5", "ac10", "ac20", "acc", "ac1sd", "bmd", "conc", "resp"
  )
  row <- as.data.frame(c(identifiers, mget(name.list)), stringsAsFactors = F)
  return(row)
}
"""

