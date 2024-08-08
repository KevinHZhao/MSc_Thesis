#' Title
#'
#' @param mtmodel
#' @param effects
#' @param peddat
#' @param emimpath
#' @param includeIm
#' @param includeIf
#' @param weinberg
#'
#' @return
#' @export
#' @keywords internal
#'
#' @examples
runEMIM <- function(mtmodel = "MS", effects = c("C", "M"), peddat,
                    emimpath = "C:/Users/Kevin/emim-v3.22-windows-x86_64/",
                    includeIm = FALSE, includeIf = FALSE, weinberg = FALSE) {
  ## Set up temp wd so EMIM files don't show up
  withr::local_dir(new = withr::local_tempdir())

  if(includeIm){
    effects <- c(effects, "Im")
  }
  if(includeIf){
    effects <- c(effects, "If")
  }

  if (all(c("C", "Im", "If") %in% effects)){
    stop("Cannot include maternal and paternal imprinting with child effects.")
  }

  ## This is for ensuring EMIM knows to use "2" as the risk allele (otherwise it
  ## will default to the least common allele)
  write(c("1", "A", "2"), "emim_rfile", ncol = 3)

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


  # Set up options for the parameter file
  options <- " -a -so -rfile emim_rfile"

  if (is.element("C", effects)) { # Multiplicative allele model for C effect
    options <- c(options, "-ct")
  }
  if (is.element("Im", effects)) { # Multiplicative allele model for C effect
    options <- c(options, "-im")
  }
  if (is.element("If", effects)) { # Multiplicative allele model for C effect
    options <- c(options, "-ip")
  }

  if (is.element("M", effects) || is.element("E:M", effects)) { # Multiplicative allele model for M effect
    options <- c(options, "-mt")
  }
  options <- paste0(options, " ", collapse = " ")


  # Note that ped data includes an E column, which must
  # must be removed
  if (is.element("E:M", effects)) {
    ## Subset the data into cases where E=0 and once with E=1
    subset0 <- subset(peddat, peddat$E == 0)
    subset1 <- subset(peddat, peddat$E == 1)

    # All
    peddat <- peddat[, c("famid", "indid", "pid", "mid", "sex", "D", "genotype1", "genotype2")]
    write.table(peddat, "temp_pedigree_all.ped", col.names = F, row.names = F, quote = F)
    write(c("1", "A", "0", "0"), "temp_pedigree_all.map", ncol = 4)

    # E=0 data
    peddat <- subset0[, c("famid", "indid", "pid", "mid", "sex", "D", "genotype1", "genotype2")]
    write.table(peddat, "temp_pedigree_0.ped", col.names = F, row.names = F, quote = F)
    write(c("1", "A", "0", "0"), "temp_pedigree_0.map", ncol = 4)

    # E=1 data
    peddat <- subset1[, c("famid", "indid", "pid", "mid", "sex", "D", "genotype1", "genotype2")]
    write.table(peddat, "temp_pedigree_1.ped", col.names = F, row.names = F, quote = F)
    write(c("1", "A", "0", "0"), "temp_pedigree_1.map", ncol = 4)

    ## Run PREMIM and EMIM
    for (i in c("all", 0, 1)) {
      # Run PREMIM
      command <- paste0(emimpath, "premim", options, "temp_pedigree_", i, ".ped temp_pedigree_", i, ".map")
      system(command, intern = TRUE)

      # Change parameter in options file for mating symmetry
      if (mtmodel == "MS") {
        params <- readLines("emimparams.dat")
        params[16] <- "0   << assume HWE and random mating (0=no=estimate 6 mu parameters, 1=yes)"
        write(params, "emimparams.dat", ncol = 1)
      }
      if (mtmodel == "MaS") {
        params <- readLines("emimparams.dat")
        params[16] <- "0   << assume HWE and random mating (0=no=estimate 6 mu parameters, 1=yes)"
        params[18] <- "1   << use CPG likelihood (estimate 9 mu parameters)"
        write(params, "emimparams.dat", ncol = 1)
      }

      # run EMIM and rename results
      command <- paste0(emimpath,
                        "emim",
                        ifelse(.Platform$OS.type == "unix",
                               "",
                               ".exe")
                        )
      system(command, intern = TRUE)
      system(paste0("mv emimsummary.out emimsummary_", i, ".out"))
    }

    ## Read in results and parse
    resAll <- read.table("emimsummary_all.out", header = T)
    res0 <- read.table("emimsummary_0.out", header = T)
    res1 <- read.table("emimsummary_1.out", header = T)

    if (is.element("C", effects)) {
      resVec["C"] <- exp(resAll$lnR1)
      pvalVec["C"] <- 2 * pnorm(abs(resAll$lnR1 / resAll$sd_lnR1), lower = F)
    }
    pvalVec["M"] <- 2 * pnorm(abs(resAll$lnS1 / resAll$sd_lnS1), lower = F)

    resVec["M[E=0]"] <- exp(res0$lnS1)
    resVec["M[E=1]"] <- exp(res1$lnS1)
    resVec["M[E=1]/M[E=0]"] <- resVec["M[E=1]"] / resVec["M[E=0]"]

    # Get a Wald-type GE test like Haplin
    z <- abs(res0$lnS1 - res1$lnS1) / sqrt(res0$sd_lnS1^2 + res1$sd_lnS1^2)
    pvalVec["E:M"] <- 2 * pnorm(z, lower = F)
  } else {
    peddat <- peddat[, c("famid", "indid", "pid", "mid", "sex", "D", "genotype1", "genotype2")]

    write.table(peddat, "temp_pedigree.ped", col.names = F, row.names = F, quote = F)
    write(c("1", "A", "0", "0"), "temp_pedigree.map", ncol = 4)


    # Run PREMIM
    command <- paste0(emimpath,
                      "premim",
                      ifelse(.Platform$OS.type == "unix",
                             "",
                             ".exe"
                             ),
                      options,
                      "temp_pedigree.ped temp_pedigree.map")
    system(command, intern = TRUE)

    # Change parameter for mating symmetry
    if (mtmodel == "MS") {
      params <- readLines("emimparams.dat")
      params[16] <- "0   << assume HWE and random mating (0=no=estimate 6 mu parameters, 1=yes)"
      write(params, "emimparams.dat", ncol = 1)
    }
    if (mtmodel == "MaS") {
      params <- readLines("emimparams.dat")
      params[16] <- "0   << assume HWE and random mating (0=no=estimate 6 mu parameters, 1=yes)"
      params[18] <- "1   << use CPG likelihood (estimate 9 mu parameters)"
      write(params, "emimparams.dat", ncol = 1)
    }

    # RUN EMIM
    command <- paste0(emimpath,
                      "emim",
                      ifelse(.Platform$OS.type == "unix",
                             "",
                             ".exe")
                      )
    system(command, intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)

    # Read in results
    res <- read.table("emimsummary.out", header = T)

    # if (is.element("M", effects)) {
    #   resVec["M"] <- exp(res$lnS1)
    #   pvalVec["M"] <- 2 * pnorm(abs(res$lnS1 / res$sd_lnS1), lower = F)
    # }
    # if (is.element("C", effects)) {
    #   resVec["C"] <- exp(res$lnR1)
    #   pvalVec["C"] <- 2 * pnorm(abs(res$lnR1 / res$sd_lnR1), lower = F)
    # }
    # if(includeI){
    #   if(MatImp){
    #     resVec["Im"] <- exp(res$lnIm)
    #     pvalVec["Im"] <- 2 * pnorm(abs(res$lnIm / res$sd_lnIm), lower = F)
    #   } else {
    #     resVec["Ip"] <- exp(res$lnIp)
    #     pvalVec["Ip"] <- 2 * pnorm(abs(res$lnIp / res$sd_lnIp), lower = F)
    #   }
    # }
  }
  #system("rm temp_pedigree*")
  return(res)
  #return(list(effects = resVec, pvals = pvalVec))
}
