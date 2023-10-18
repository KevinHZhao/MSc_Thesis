# Run a simulation with the specified parameters.
# eligible effects are "M","C", "M:E", "D" (terms included in model)
# eligible MT models are "HW", "MS", "MaS"
# includeE, PopStrat, and Control means that these conditions are added to
#   the dataset. However, to include an environmental effect to the analysis
#   must have V not be all 0 and to include controls in analysis than a "D"
#   term will need to be added to the effects.
runTrioReps <- function(nreps = 100, ntrios = 1000, maf = 0.3,
                        R = c(1, 1, 1), S = c(1, 1, 1), V = c(1, 1, 1),
                        mtCoef = c(1, 1, 1), effects = c("M", "C"),
                        includeE = FALSE, envint = "Mother", prE = 0,
                        includePopStrat = FALSE, numPop = 1, Fst = 0.005, prCase.byPop = NULL, prControl.byPop = NULL,
                        includeControl = FALSE, prControl = 0, prE.control = prE,
                        Loglin = TRUE, loglinMT = c("HW", "MS"),
                        Haplin = TRUE, haplinControls = FALSE,
                        EMIM = TRUE, EMIMMT = c("HW", "MS"),
                        PStest = FALSE, HWtest = FALSE, MStest = FALSE) {
  # Set up results tables. There will be two tables for each method that
  # will be run. There may be objects for HW, mating symmetry and population
  # stratification tests.
  pvalPS <- NULL
  resLoglinPS <- NULL
  pvalLoglinPS <- NULL
  resLoglin <- NULL
  pvalLoglin <- NULL
  resHaplin <- NULL
  pvalHaplin <- NULL
  resEMIM <- NULL
  pvalEMIM <- NULL
  HWres <- NULL
  MSres <- NULL

  if (Loglin) {
    ncolLoglin <- length(effects) * length(loglinMT)
    temp <- expand.grid(effects, loglinMT)[, 2:1]
    namesLoglin <- paste(temp[, 1], temp[, 2], sep = ".")
    resLoglin <- matrix(nrow = nreps, ncol = ncolLoglin)
    colnames(resLoglin) <- namesLoglin
    pvalLoglin <- matrix(nrow = nreps, ncol = ncolLoglin)
    colnames(pvalLoglin) <- namesLoglin

    if (PStest == TRUE) {
      if (prControl == 0) {
        stop("PStest only possible when there are controls\n")
      } else {
        pvalPS <- matrix(nrow = nreps, ncol = length(loglinMT))
        colnames(pvalPS) <- loglinMT
        resLoglinPS <- matrix(nrow = nreps, ncol = ncolLoglin)
        colnames(resLoglinPS) <- namesLoglin
        pvalLoglinPS <- matrix(nrow = nreps, ncol = ncolLoglin)
        colnames(pvalLoglinPS) <- namesLoglin
      }
    }
  }

  if (Haplin) {
    if (is.element("E:M", effects)) {
      neweffects <- c(setdiff(effects, c("E:M", "M")), "M[E=0]", "M[E=1]", "M[E=1]/M[E=0]")
      ncolHaplin <- length(neweffects)
      resHaplin <- matrix(nrow = nreps, ncol = ncolHaplin)
      colnames(resHaplin) <- neweffects
      pvalHaplin <- matrix(nrow = nreps, ncol = length(effects))
      colnames(pvalHaplin) <- effects
    } else {
      ncolHaplin <- length(effects)
      namesHaplin <- effects
      resHaplin <- matrix(nrow = nreps, ncol = ncolHaplin)
      colnames(resHaplin) <- namesHaplin
      pvalHaplin <- matrix(nrow = nreps, ncol = ncolHaplin)
      colnames(pvalHaplin) <- namesHaplin
    }
  }

  if (EMIM) {
    if (is.element("E:M", effects)) {
      neweffects <- c(setdiff(effects, c("E:M", "M")), "M[E=0]", "M[E=1]", "M[E=1]/M[E=0]")
      ncolEMIM <- length(EMIMMT) * length(neweffects)
      temp.pval <- expand.grid(effects, EMIMMT)[, 2:1]
      temp.res <- expand.grid(neweffects, EMIMMT)[, 2:1]
      resEMIM <- matrix(nrow = nreps, ncol = ncolEMIM)
      colnames(resEMIM) <- paste(temp.res[, 1], temp.res[, 2], sep = ".")
      pvalEMIM <- matrix(nrow = nreps, ncol = nrow(temp.pval))
      colnames(pvalEMIM) <- paste(temp.pval[, 1], temp.pval[, 2], sep = ".")
    } else {
      ncolEMIM <- length(EMIMMT) * length(effects)
      temp <- expand.grid(effects, EMIMMT)[, 2:1]
      namesEMIM <- paste(temp[, 1], temp[, 2], sep = ".")
      resEMIM <- matrix(nrow = nreps, ncol = ncolEMIM)
      colnames(resEMIM) <- namesEMIM
      pvalEMIM <- matrix(nrow = nreps, ncol = ncolEMIM)
      colnames(pvalEMIM) <- namesEMIM
    }
  }

  if (HWtest == TRUE) {
    if (includeControl == TRUE) {
      HWres <- matrix(nrow = nreps, ncol = 6)
      colnames(HWres) <- c(
        "Mothers.Case", "Fathers.Case", "Both.Case",
        "Mothers.Control", "Fathers.Control", "Both.Control"
      )
    } else {
      # browser()
      HWres <- matrix(nrow = nreps, ncol = 3)
      colnames(HWres) <- c("Mothers", "Fathers", "Both")
    }
  }

  if (MStest == TRUE) {
    MSres <- vector(length = nreps)
  }


  # Run the simulation repetitions
  for (i in 1:nreps) {
    dat <- simulateData(
      ntrios = ntrios, maf = maf,
      R = R, S = S, V = V,
      mtCoef = mtCoef, mtmodel = mtmodel,
      includeE = includeE, envint = envint, prE = prE,
      includePopStrat = includePopStrat, numPop = numPop, Fst = Fst,
      includeControl = includeControl, prControl = prControl, prE.control = prE.control,
      prCase.byPop = prCase.byPop, prControl.byPop = prControl.byPop
    )

    if (HWtest == TRUE) {
      if (includeControl == TRUE) {
        HWres[i, 1:3] <- getHWE(subset(dat$dat4R, dat$dat4R$D == 1))
        HWres[i, 4:6] <- getHWE(subset(dat$dat4R, dat$dat4R$D == 0))
      } else {
        HWres[i, ] <- getHWE(dat$dat4R)
      }
    }

    if (MStest == TRUE) {
      MSres[i] <- MaSlrt(dat$dat4R)
    }


    # For the ith dataset, run the log linear approach if requested
    if (Loglin) {
      for (j in 1:length(loglinMT)) { # Loop over mating models selected
        res <- runLoglin(loglinMT[j], effects, dat$dat4R, PStest = PStest)
        colvals <- paste(loglinMT[j], names(res$effects), sep = ".")
        names(res$effects) <- colvals
        names(res$pvals) <- colvals
        resLoglin[i, colvals] <- res$effects[colvals]
        pvalLoglin[i, colvals] <- res$pvals[colvals]

        if (PStest == TRUE) {
          # browser()
          names(res$effectsPS) <- colvals
          names(res$pvalsPS) <- colvals
          resLoglinPS[i, colvals] <- res$effectsPS[colvals]
          pvalLoglinPS[i, colvals] <- res$pvalsPS[colvals]
          pvalPS[i, loglinMT[j]] <- res$PS.test
        }
      }
    }

    # For the ith dataset, run haplin if requested
    if (Haplin) {
      res <- runHaplin(effects, dat$dat4haplin, haplinControls = haplinControls)
      colvals.effects <- names(res$effects)
      colvals.pvals <- names(res$pvals)
      resHaplin[i, colvals.effects] <- res$effects[colvals.effects]
      pvalHaplin[i, colvals.pvals] <- res$pvals[colvals.pvals]
    }

    # For the ith dataset, run EMIM if requested
    if (EMIM) {
      for (j in 1:length(EMIMMT)) { # Loop over mating models selected
        res <- runEMIM(EMIMMT[j], effects, dat$dat4EMIM)
        colvals.res <- paste(EMIMMT[j], names(res$effects), sep = ".")
        colvals.pvals <- paste(EMIMMT[j], names(res$pvals), sep = ".")
        names(res$effects) <- colvals.res
        names(res$pvals) <- colvals.pvals
        resEMIM[i, colvals.res] <- res$effects[colvals.res]
        pvalEMIM[i, colvals.pvals] <- res$pvals[colvals.pvals]
      }
    }
  }

  return(list(
    resLoglin = resLoglin, pvalLoglin = pvalLoglin,
    resLoglinPS = resLoglinPS, pvalLoglinPS = pvalLoglinPS,
    resHaplin = resHaplin, pvalHaplin = pvalHaplin,
    resEMIM = resEMIM, pvalEMIM = pvalEMIM,
    pvalHW = HWres, pvalMaS = MSres, pvalPS = pvalPS
  ))
}
