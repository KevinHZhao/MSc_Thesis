#' Simulate data for TriLLIEM
#'
#' @param ntrios Number of trios to simulate.
#' @param maf Minor allele frequency in the population, a proportion between 0 and 1.
#' @param R Vector of 3 elements representing child effects for 0, 1, and 2
#' copies of the risk allele, respectively.
#' @param S Vector of 3 elements representing maternal effects for 0, 1, and 2
#' copies of the risk allele, respectively.
#' @param V Vector of 3 elements representing gene-environment effects
#' for 0, 1, and 2 copies of the risk allele, respectively.
#' @param mtCoef Mating type coefficients.
#' @param mtmodel Mating type of the population, can be "`HWE`" for Hardy-Weinberg
#' Equilibrium, "`MS`" for Mating Symmetry, and "`MaS`" for Mating Asymmetry.
#' @param Im Maternal imprinting effect.
#' @param If Paternal imprinting effect.
#' @param includeE A logical value indicating whether environmental effects
#' should be included in the simulation.
#' @param envint If set to "`Mother`", simulates maternal gene-environment
#' interactions, otherwise simulates child gene-environmnet interactions.
#' @param prE Probability of a trio to have the environmental effects.
#' @param includeControl A logical value indicating whether controls should
#' be included in the simulations.
#' @param prControl Probability of a trio to be a control trio.
#' @param prE.control Probablity of a trio to be a control trio with the
#' environmental effect.
#' @param includePopStrat A logical value indicating whether to include
#' population stratification in the simulation.
#' @param numPop
#' @param Fst
#' @param prCase.byPop
#' @param prControl.byPop
#'
#' @return A data frame of the same format as [example_dat4R]
#' @export
#'
#' @examples
#' ## Simulating data with multiplicative maternal effect of 2, and paternal imprinting of 3.
#' simulateData(S = c(1, 2, 4), If = 3)
simulateData <- function(ntrios = 1000, maf = 0.3,
                         R = c(1, 1, 1), S = c(1, 1, 1), V = c(1, 1, 1),
                         mtCoef = c(1, 1, 1), mtmodel = "MS",
                         Im = 1, If = 1,
                         includeE = FALSE, envint = "Mother", prE = 0,
                         includeControl = FALSE, prControl = 0, prE.control = prE,
                         includePopStrat = FALSE, numPop = 1, Fst = 0.005,
                         prCase.byPop = NULL, prControl.byPop = NULL) {
  # prE.control should be a random binomial parameter, NOT a ratio
  # Start off small with tests
  if(Im != 1 && If != 1 && !identical(R, c(1,1,1))){
    warning("Maternal and paternal imprinting included with child effects,
            resulting data cannot have all parameters simultaneously modelled.")
  }

  if (includePopStrat == TRUE) {
    # Ensure that the proportion of the cases from each population is provided.
    if (round(sum(prCase.byPop), 8) != 1) {
      stop("The sum of prCase.byPop must equal 1. Each element is the proportion
              of the case trios from each subpopulation.\n")
    }

    if ((includeControl == TRUE) & (round(sum(prControl.byPop), 8) != 1)) {
      stop("The sum of prControl.byPop must equal 1. Each element is the proportion
              of the case trios from each subpopulation.\n")
    }
  } else {
    prCase.byPop <- 1
    prControl.byPop <- 1
  }

  # Ensure we know proportion of families that are control trios
  if (includeControl == TRUE) {
    if ((prControl <= 0) || (prControl >= 1)) {
      stop("Proportion of trios that are control trios must be between 0 and 1\n")
    }
  } else {
    prControl <- 0 # In case this was not 0 by accident
  }


  # Check proportion of environmental variable
  if (includeE == TRUE) {
    if ((sum(prE <= 0) > 0) || (sum(prE >= 1) > 0)) {
      stop("Probabilities for environmental variables for cases must be between 0 and 1\n")
    }
    if ((sum(prE.control <= 0) > 0) || (sum(prE.control >= 1) > 0)) {
      stop("Probabilities for environmental variables for controls must be between 0 and 1\n")
    }
  } else {
    prE <- 0 # In case these were not 0 by accident
    prE.control <- prE
  }


  # If only a single frequency of environmental variable is given,
  # make all populations have the same frequency
  if ((includePopStrat == TRUE) && (length(prE) == 1)) {
    prE <- rep(prE, numPop)
    prE.control <- rep(prE.control, numPop)
  }


  # Use Balding-Nichols for MAF in subpopulations
  if ((includePopStrat == TRUE) && (length(maf) == 1)) {
    warning("MAF of subpopulations not provided. Will use Balding-Nichols model\n")
    alpha <- maf * (1 - Fst) / Fst
    beta <- (1 - maf) * (1 - Fst) / Fst

    # Generate allele frequencies for the two populations
    q <- rbeta(numPop, shape1 = alpha, shape2 = beta)
  } else {
    q <- maf
  }


  # Count tables across all populations
  caseE0.all <- NULL
  caseE1.all <- NULL
  controlE0.all <- NULL
  controlE1.all <- NULL
  full.tables <- NULL

  # Create all the datasets separately in each population
  for (i in 1:numPop) {
    ntrios.pop.case <- round(ntrios * (1 - prControl) * prCase.byPop[i])
    ntrios.pop.control <- round(ntrios * prControl * prControl.byPop[i])

    ntrios.pop.notE.case <- round(ntrios.pop.case * (1 - prE[i]))
    ntrios.pop.E.case <- round(ntrios.pop.case * prE[i])
    ntrios.pop.notE.control <- round(ntrios.pop.control * (1 - prE.control[i]))
    ntrios.pop.E.control <- round(ntrios.pop.control * prE.control[i])


    # Simulate "case" trios,
    caseE0 <- simulateDataSubset(
      ntrios = ntrios.pop.notE.case, maf = q[i], R = R, S = S, mtCoef = mtCoef,
      includeE = FALSE, Im = Im, If = If
    )
    if (i == 1) {
      caseE0.all <- caseE0
      full.tables <- list(condition = c(paste0("pop=", i), "case=1", "E=0"), table = caseE0$datFull)
    } else {
      caseE0.all <- mergeCounts(caseE0.all, caseE0)
      full.tables <- c(full.tables, list(
        condition = c(paste0("pop=", i), "case=1", "E=0"),
        table = caseE0$datFull
      ))
    }

    if (includeE == TRUE) {
      caseE1 <- simulateDataSubset(
        ntrios = ntrios.pop.E.case, maf = q[i], R = R, S = S, mtCoef = mtCoef,
        V = V, includeE = TRUE, envint = envint, Im = Im, If = If
      )
      if (i == 1) {
        caseE1.all <- caseE1
      } else {
        caseE1.all <- mergeCounts(caseE1.all, caseE1)
      }
      full.tables <- c(full.tables, list(
        condition = c(paste0("pop=", i), "case=1", "E=1"),
        table = caseE1$datFull
      ))
    }

    # Simulate control trios
    if (includeControl == TRUE) {
      controlE0 <- simulateDataSubset(
        ntrios = ntrios.pop.notE.control, maf = q[i],
        R = c(1, 1, 1), S = c(1, 1, 1),
        mtCoef = mtCoef, includeE = FALSE, includeControl = TRUE
      )
      if (i == 1) {
        controlE0.all <- controlE0
      } else {
        controlE0.all <- mergeCounts(controlE0.all, controlE0)
      }
      full.tables <- c(full.tables, list(
        condition = c(paste0("pop=", i), "case=0", "E=0"),
        table = controlE0$datFull
      ))

      if (includeE == TRUE) {
        controlE1 <- simulateDataSubset(
          ntrios = ntrios.pop.E.control, maf = q[i],
          R = c(1, 1, 1), S = c(1, 1, 1), mtCoef = mtCoef,
          V = c(1, 1, 1), includeE = TRUE, envint = envint,
          includeControl = TRUE
        )
        if (i == 1) {
          controlE1.all <- controlE1
        } else {
          controlE1.all <- mergeCounts(controlE1.all, controlE1)
        }
        full.tables <- c(full.tables, list(
          condition = c(paste0("pop=", i), "case=0", "E=1"),
          table = controlE1$datFull
        ))
      }
    }
  }

  # Create final dataset by stacking all the individual tables (E=T/F, control=T/F)
  finaldat <- list(dat4R = caseE0.all$dat4R, dat4haplin = caseE0.all$dat4haplin)

  if (includeE == TRUE) {
    finaldat <- stackCounts(finaldat, caseE1.all)
  }
  if (includeControl == TRUE) {
    finaldat <- stackCounts(finaldat, controlE0.all)
    if (includeE == TRUE) {
      finaldat <- stackCounts(finaldat, controlE1.all)
    }
  }

  peddat <- createPed(finaldat$dat4R)

  return(list(
    dat4R = finaldat$dat4R, dat4haplin = finaldat$dat4haplin, dat4EMIM = peddat,
    datAll = full.tables
  ))
}
