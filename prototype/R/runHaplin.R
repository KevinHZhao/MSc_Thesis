## Run a haplin analysis with the specified model
runHaplin <- function(effects, dat, haplinControls = FALSE) {
  # Setup results objects
  if (is.element("E:M", effects)) {
    neweffects <- c(setdiff(effects, c("E:M", "M")), "M[E=0]", "M[E=1]", "M[E=1]/M[E=0]")
    resVec <- vector(length = length(neweffects))
    pvalVec <- vector(length = length(effects))
    names(resVec) <- neweffects
    names(pvalVec) <- effects
  } else {
    resVec <- vector(length = length(effects))
    names(resVec) <- effects
    pvalVec <- vector(length = length(effects))
    names(pvalVec) <- effects
  }

  # Number of columns that are not genotype. First column is E, second is Phenotype
  nvars <- 2

  # Write out the data for haplin, then read it in using haplin
  write.table(dat, "haplin_temp.dat", quote = F, row = F, col = F, sep = " ")
  write(c("chrom snp a", "1 rs1 0"), "temp.map", ncol = 1)
  dat.raw <- invisible(genDataRead("haplin_temp.dat",
    format = "haplin", overwrite = TRUE, n.vars = nvars,
    map.file = "temp.map"
  ))

  if (is.element("E:M", effects)) {
    if (haplinControls) {
      dat.processed <- invisible(genDataPreprocess(
        data.in = dat.raw, design = "cc.triad",
        overwrite = TRUE, map.file = "temp.map"
      ))
      res <- invisible(haplinStrat(dat.processed,
        response = "mult", verbose = FALSE,
        design = "cc.triad", ccvar = 2,
        printout = FALSE, strata = 1, reference = "ref.cat", maternal = TRUE
      ))
    } else {
      dat.processed <- invisible(genDataPreprocess(
        data.in = dat.raw, design = "triad",
        overwrite = TRUE, map.file = "temp.map"
      ))
      res <- invisible(haplinStrat(dat.processed,
        response = "mult", verbose = FALSE,
        printout = FALSE, strata = 1, reference = "ref.cat", maternal = TRUE
      ))
    }
    GEtest <- gxe(res)
    resVec["M[E=1]"] <- haptable(res)[6, "RRm.est."] # RR for E=1, M=1
    resVec["M[E=0]"] <- haptable(res)[4, "RRm.est."] # RR for E=0, M=1
    resVec["M[E=1]/M[E=0]"] <- resVec["M[E=1]"] / resVec["M[E=0]"]
    pvalVec["E:M"] <- GEtest$gxe.test[3, "pval"] # pval for stratified test
    pvalVec["M"] <- haptable(res)[2, "RRm.p.value"] # pval for unstratified analysis

    if (is.element("C", effects)) { # Using the RR and p-value for the unstratified analysis
      resVec["C"] <- haptable(res)[2, "RR.est."]
      pvalVec["C"] <- haptable(res)[2, "RR.p.value"]
    }
  } else {
    if (is.element("M", effects)) {
      includeMaternal <- TRUE
    } else {
      includeMaternal <- FALSE
    }

    if (haplinControls == TRUE) {
      dat.processed <- invisible(genDataPreprocess(
        data.in = dat.raw, design = "cc.triad",
        overwrite = TRUE, map.file = "temp.map"
      ))
      res <- invisible(haplin(dat.processed,
        response = "mult", verbose = FALSE, design = "cc.triad",
        ccvar = 2, printout = FALSE, reference = "ref.cat", maternal = includeMaternal
      ))
    } else {
      dat.processed <- invisible(genDataPreprocess(
        data.in = dat.raw, design = "triad",
        overwrite = TRUE, map.file = "temp.map"
      ))
      res <- invisible(haplin(dat.processed,
        response = "mult", verbose = FALSE,
        printout = FALSE, reference = "ref.cat", maternal = includeMaternal
      ))
    }

    res <- haptable(res)[2, ]
    if (is.element("M", effects)) {
      resVec["M"] <- res$RRm.est.
      pvalVec["M"] <- res$RRm.p.value
    }
    if (is.element("C", effects)) {
      resVec["C"] <- res$RR.est.
      pvalVec["C"] <- res$RR.p.value
    }
  }

  system("rm temp.map haplin_temp.dat")
  return(list(effects = resVec, pvals = pvalVec))
}
