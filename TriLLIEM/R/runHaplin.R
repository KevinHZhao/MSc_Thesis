#' Title
#'
#' @param effects
#' @param dat
#' @param haplinControls
#' @param PoO
#' @param verbose
#'
#' @return
#' @export
#' @keywords internal
#'
#' @examples
runHaplin <- function(effects = c("C", "M"), dat, haplinControls = FALSE,
                      includeE = FALSE, PoO = FALSE, verbose = FALSE) {
  ## Test with pooxe
  ## Set up temp wd so Haplin files don't show up
  withr::local_dir(new = withr::local_tempdir())

  # Setup results objects
  if (includeE) {
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
  dat.raw <- invisible(Haplin::genDataRead("haplin_temp.dat",
    format = "haplin", overwrite = TRUE, n.vars = nvars,
    map.file = "temp.map"
  ))

  if (includeE) {
    if (haplinControls) {
      dat.processed <- invisible(Haplin::genDataPreprocess(
        data.in = dat.raw, design = "cc.triad",
        overwrite = TRUE, map.file = "temp.map"
      ))
      res <- invisible(Haplin::haplinStrat(dat.processed,
        response = "mult", verbose = verbose,
        design = "cc.triad", ccvar = 2,
        printout = FALSE, strata = 1, reference = "ref.cat", maternal = TRUE,
        poo = PoO
      ))
    } else {
      dat.processed <- invisible(Haplin::genDataPreprocess(
        data.in = dat.raw, design = "triad",
        overwrite = TRUE, map.file = "temp.map"
      ))
      res <- invisible(Haplin::haplinStrat(dat.processed,
        response = "mult", verbose = verbose,
        printout = FALSE, strata = 1, reference = "ref.cat", maternal = TRUE,
        poo = PoO
      ))
    }
    GEtest <- Haplin::gxe(res)
    resVec["M[E=1]"] <- Haplin::haptable(res)[6, "RRm.est."] # RR for E=1, M=1
    resVec["M[E=0]"] <- Haplin::haptable(res)[4, "RRm.est."] # RR for E=0, M=1
    resVec["M[E=1]/M[E=0]"] <- resVec["M[E=1]"] / resVec["M[E=0]"]
    pvalVec["E:M"] <- GEtest$gxe.test[3, "pval"] # pval for stratified test
    pvalVec["M"] <- Haplin::haptable(res)[2, "RRm.p.value"] # pval for unstratified analysis

    if (is.element("C", effects)) { # Using the RR and p-value for the unstratified analysis
      resVec["C"] <- Haplin::haptable(res)[2, "RR.est."]
      pvalVec["C"] <- Haplin::haptable(res)[2, "RR.p.value"]
    }
  } else {
    if (is.element("M", effects)) {
      includeMaternal <- TRUE
    } else {
      includeMaternal <- FALSE
    }

    if (haplinControls == TRUE) {
      dat.processed <- invisible(Haplin::genDataPreprocess(
        data.in = dat.raw, design = "cc.triad",
        overwrite = TRUE, map.file = "temp.map"
      ))
      res <- invisible(Haplin::haplin(dat.processed,
        response = "mult", verbose = verbose, design = "cc.triad",
        ccvar = 2, printout = FALSE, reference = "ref.cat", maternal = includeMaternal,
        poo = PoO
      ))
    } else {
      dat.processed <- invisible(Haplin::genDataPreprocess(
        data.in = dat.raw, design = "triad",
        overwrite = TRUE, map.file = "temp.map"
      ))
      res <- invisible(Haplin::haplin(dat.processed,
        response = "mult", verbose = verbose,
        printout = FALSE, reference = "ref.cat", maternal = includeMaternal,
        poo = PoO
      ))
    }

    # res <- Haplin::haptable(res)[2, ]
  }
  # system("rm temp.map haplin_temp.dat")
  return(res)
  # return(list(effects = resVec, pvals = pvalVec))
}
