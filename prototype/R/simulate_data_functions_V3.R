## Functions needed for maternal effect analyses in R

createGenoMat <- function() {
  M <- c(rep(0, 3), rep(1, 2), 0, 2, rep(1, 5), rep(2, 3))
  F <- c(0, rep(1, 2), rep(0, 2), 2, 0, rep(1, 3), rep(2, 2), rep(1, 2), 2)
  C <- c(rep(0, 2), 1, 0, rep(1, 3), 0, 1, 2, 1, 2, 1, rep(2, 2))
  return(data.frame(M, F, C))
}

createHaplinGeno <- function() {
  M <- c(rep("1;1", 3), rep("1;2", 2), "1;1", "2;2", rep("1;2", 5), rep("2;2", 3))
  F <- c("1;1", rep("1;2", 2), rep("1;1", 2), "2;2", "1;1", rep("1;2", 3), rep("2;2", 2), rep("1;2", 2), "2;2")
  C <- c(rep("1;1", 2), "1;2", "1;1", rep("1;2", 3), "1;1", "1;2", "2;2", "1;2", "2;2", "1;2", rep("2;2", 2))
  return(data.frame(M, F, C))
}

mtmat <- function(maf = 0.4, C = c(1, 1, 1)) {
  genocat <- createGenoMat()
  mts <- unique(genocat[, 1:2])

  Mprobs <- dbinom(mts[, 1], 2, prob = maf)
  Fprobs <- dbinom(mts[, 2], 2, prob = maf)

  asymfactor <- c(1, C[1], 2 - C[1], C[2], 2 - C[2], 1, C[3], 2 - C[3], 1)
  mts <- cbind(mts, prMF = Mprobs * Fprobs * asymfactor)

  mt_MS <- c(1, 2, 2, 3, 3, 4, 5, 5, 6)
  mt_MaS <- 1:9

  mts <- cbind(mt_MS, mt_MaS, mts)

  prCGivenMF <- c(1, rep(1 / 2, 4), 1, 1, 1 / 4, 1 / 2, 1 / 4, rep(1 / 2, 4), 1)

  genocat <- cbind(genocat, prCGivenMF)
  temp <- merge(genocat, mts)
  temp$prMFC <- temp$prMF * temp$prCGivenMF

  return(temp[
    order(temp$mt_MaS),
    c("mt_MS", "mt_MaS", "M", "F", "C", "prMF", "prCGivenMF", "prMFC")
  ])
}



# A simplified function for simulating a subset of the full data
# under particular conditions
# - Case Trios
# - Control Trios
# - With E or without
# - With MaS or not
simulateDataSubset <- function(ntrios = 1000, maf = 0.3,
                               R = c(1, 1, 1), S = c(1, 1, 1),
                               mtCoef = c(1, 1, 1),
                               V = c(1, 1, 1), includeE = FALSE, envint = "Mother",
                               includeControl = FALSE) {
  genomat <- mtmat(maf, C = mtCoef)

  # Compute the values proportional to P(D|M,F,C,E) for the
  # selected model.

  # Add main genetic effects
  diseasefactor <- R[(genomat$C + 1)] * S[(genomat$M + 1)]

  # Add environmental effects for trios that will have E=1
  if (envint == "Mother") {
    diseasefactor <- diseasefactor * V[(genomat$M + 1)]
  } else {
    diseasefactor <- diseasefactor * V[(genomat$C + 1)]
  }

  # Compute the values proportional to P(M,F,C,E|D) = P(D|M,F,C,E)P(M,F)P(C|M,F)
  # and re-scale so that probabilities sum to 1 (if only cases this happens)
  diseaseprob <- genomat$prMFC * diseasefactor
  diseaseprob <- diseaseprob / sum(diseaseprob)


  # Sample the trios
  triocounts <- sample(1:length(diseaseprob), ntrios, prob = diseaseprob, replace = TRUE)
  counts <- data.frame(table(triocounts))
  colnames(counts) <- c("type", "count")

  # Prepare data for R log-linear
  genomat <- data.frame(type = 1:length(diseaseprob), genomat, PrCMFEgivenD = diseaseprob)
  triodat <- merge(genomat, counts, all = TRUE, by = "type")
  triodat[is.na(triodat[, "count"]), "count"] <- 0

  if (includeE == TRUE) {
    E <- rep(1, nrow(triodat))
  } else {
    E <- rep(0, nrow(triodat))
  }

  if (includeControl == TRUE) {
    D <- rep(0, nrow(triodat))
  } else {
    D <- rep(1, nrow(triodat))
  }

  triodat <- data.frame(triodat[, 1:6], E = E, D = D, triodat[, 7:ncol(triodat)])
  subdat <- data.frame(triodat[, c(1:8)], count = triodat$count)


  # Prepare data for Haplin. Haplin has 1 row per trio, with the row
  # consisting of M,F,C genotype columns in format 1;1, 1;2 or 2;2
  # The first column will be the environmental covariate
  # The second column will be the disease status
  # The third, fourth and fifth columns are the genotype for M, F, and C
  haplingeno <- createHaplinGeno()
  # if (includeControl==TRUE){ # Only parents so make child genotype NA
  #  MF=unique(haplingeno[,1:2])
  #  haplingeno=data.frame(type=1:nrow(MF),MF,C=rep(NA,nrow(MF)))
  # } else{
  haplingeno <- data.frame(type = 1:15, haplingeno)
  # }
  haplindat <- merge(haplingeno, triodat, by = "type")


  finaldat <- NULL
  for (j in 1:nrow(haplindat)) {
    if (haplindat[j, "count"] > 0) {
      genos <- matrix(unlist(rep(haplindat[j, 2:4], haplindat[j, "count"])), ncol = 3, byrow = TRUE)
      genos <- cbind(
        rep(haplindat[j, "E"], length = haplindat[j, "count"]),
        rep(haplindat[j, "D"], length = haplindat[j, "count"]),
        genos
      )
      finaldat <- rbind(finaldat, genos)
    }
  }
  finaldat <- as.data.frame(finaldat)

  # If D=0, remove the genotypes for the children since these
  # are not used for log-linear and EMIM
  # if (includeControl==TRUE){
  #  subdat$C=rep(-9,nrow(subdat))
  # }

  return(list(dat4R = subdat, dat4haplin = finaldat, datFull = triodat))
}




# To be used when there is hidden population substructure. This
# function merges the count column for two hidden populations.
mergeCounts <- function(dat1, dat2) {
  triodat1 <- dat1$dat4R
  triodat2 <- dat2$dat4R

  # Merge for log linear analysis
  alltriodat <- merge(triodat1, triodat2, by = "type", all = TRUE)
  counts <- alltriodat$count.x + alltriodat$count.y
  alltriodat <- data.frame(triodat1[, c(1:8)], count = counts)

  # Stack for haplin
  allhaplindat <- rbind(dat1$dat4haplin, dat2$dat4haplin)


  return(list(dat4R = alltriodat, dat4haplin = allhaplindat))
}

# To be used with environmental effects and control trios where
# the datasets with E=0 or 1 and D=0 or 1 are stacked.
stackCounts <- function(dat1, dat2) {
  triodat1 <- dat1$dat4R
  triodat2 <- dat2$dat4R

  # Stack for log linear analysis
  alltriodat <- rbind(triodat1, triodat2)

  # Stack for haplin
  allhaplindat <- rbind(dat1$dat4haplin, dat2$dat4haplin)

  return(list(dat4R = alltriodat, dat4haplin = allhaplindat))
}

# Make a Paternal origin & Maternal origin boolean column
# Can limit things (ex imprinting only for cases no controls)
# can't really have both Im and If in simulating (overparameterized, Ainsworth)

# This function creates a full dataset possibly with
# - control trios
# - case trios
# - environmental variable
# - population stratification
simulateData <- function(ntrios = 1000, maf = 0.3,
                         R = c(1, 1, 1), S = c(1, 1, 1), V = c(1, 1, 1),
                         mtCoef = c(1, 1, 1), mtmodel = "MS",
                         impM = 1, impF = 1,
                         includeE = FALSE, envint = "Mother", prE = 0,
                         includeControl = FALSE, prControl = 0, prE.control = prE,
                         includePopStrat = FALSE, numPop = 1, Fst = 0.005,
                         prCase.byPop = NULL, prControl.byPop = NULL) {
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
      includeE = FALSE
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
        V = V, includeE = TRUE, envint = envint
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

# Create a ped file
# Family ID
# Individual ID
# Paternal ID
# Maternal ID
# Sex (1=male; 2=female; other=unknown)
# Phenotype
createPed <- function(dat) {
  # Create the fixed columns
  ntrios <- sum(dat[, "count"])

  famid <- as.vector(t(matrix(rep(1:ntrios, 3), ncol = 3, byrow = F)))
  indid <- seq(1:length(famid))
  pid <- rep(0, length(famid))
  pid[3 * (1:ntrios)] <- indid[3 * (1:ntrios) - 2]
  mid <- rep(0, length(famid))
  mid[3 * (1:ntrios)] <- indid[3 * (1:ntrios) - 1]
  sex <- as.vector(rbind(
    rep(1, ntrios),
    rep(2, ntrios),
    sample(c(1, 2), ntrios, replace = TRUE)
  ))


  # Create the columns with data that changes from row
  # to row
  finaldat <- NULL

  for (i in 1:nrow(dat)) {
    ntrios <- dat[i, "count"]

    if (ntrios > 0) {
      phenotype <- rep(c(-9, -9, (dat[i, "D"] + 1)), ntrios)
      E <- rep(dat[i, "E"], 3 * ntrios)
      genotype1 <- unlist(rep(dat[i, c("F", "M", "C")], ntrios))
      genotype1[genotype1 == 2] <- 1
      genotype2 <- unlist(rep(dat[i, c("F", "M", "C")], ntrios)) - 1
      genotype2[genotype2 == -1] <- 0
      genotype1 <- genotype1 + 1
      genotype2 <- genotype2 + 1

      # Deal with missing genotypes. They are coded -9, so will be changed
      # to -8 by above code. Use "0" for missing
      # genotype1[genotype1==-8]=0
      # genotype2[genotype2==-8]=0

      tempdat <- data.frame(E, D = phenotype, genotype1, genotype2)

      finaldat <- rbind(finaldat, tempdat)
    }
    tempdat <- NULL
  }

  finaldat <- data.frame(famid, indid, pid, mid, sex, finaldat)

  return(finaldat)
}
