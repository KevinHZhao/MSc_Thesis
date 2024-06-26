#' Running a loglinear analysis of trio data.
#'
#' @param mtmodel Mating type model to use in the analysis, can be "`HWE`" for
#' Hardy-Weinberg Equilibrium, "`MS`" for Mating Symmetry, and "`MaS`" for Mating Asymmetry.
#' @param effects A vector listing the effects, as strings, to include
#' in the model.  Example effects include:
#' \describe{
#'  \item{"`C`"}{Child effects.}
#'  \item{"`M`"}{Maternal effects.}
#'  \item{"`Im`"}{Maternal imprinting effects.}
#'  \item{"`If`"}{Paternal imprinting effects.}
#'  \item{"`E:M`"}{Maternal gene-environment effects.}
#' }
#' @param dat A data frame with triad data, with the formatting of [example_dat4R].
#' @param PStest A logical value indicating whether to performa population
#' stratification test on the data.
#' @param includeIm A logical value indicating whether to include maternal
#' imprinting in the model, equivalent to adding "`Im`" in the "`effects`" vector.
#' @param includeIf A logical value indicating whether to include paternal
#' imprinting in the model, equivalent to adding "`If`" in the "`effects`" vector.
#' @param Minit Initial proportion of maternal inheritence to split the triple
#' heterozygote cell by if the EM algorithm is necessary.
#' @param max.iter Maximum number of iterations for the EM algorithm.
#' @param EM.diag A logical value indicating whether to show diagnostic messages
#' for the EM algorithm.
#'
#' @return An object of class "`glm`".
#' @export
#'
#' @examples
#' res <- TriLLIEM(mtmodel = "HWE", effects = c("C", "M", "Im"), dat = example_dat4R)
#' summ_trill(res)
TriLLIEM <- function(mtmodel = "MS", effects = c("C", "M"), dat, PStest = FALSE,
                     includeIm = FALSE, includeIf = FALSE, Minit = 0.5, max.iter = 30,
                     EM.diag = FALSE) {
  if(includeIm){
    effects <- c(effects, "Im")
  }
  if(includeIf){
    effects <- c(effects, "If")
  }

  if (all(c("C", "Im", "If") %in% effects)){
    stop("Cannot include maternal and paternal imprinting with child effects.")
  }

  # Portion of model equation and offset depends on mating type model
  origDat <- add_PoO_data(dat, Mprop = Minit)

  # Portion of model equation depends on mating type model
  # No offset as we split the (1,1,1) case
  if (mtmodel == "HWE") {
    origDat$HWgeno <- origDat$M + origDat$F
    mteffect <- "HWgeno"
    modelformula <- "count~" # Must include intercept for HW model because of log(1-p) term
  } else if (mtmodel == "MS") {
    mteffect <- "as.factor(mt_MS)"
    modelformula <- "count~-1+" # I think I can remove the intercept for MS model
  } else if (mtmodel == "MaS") {
    if (length(unique(origDat$D)) == 1) {
      stop("Only 1 phenotype in the phenotype column. Mating asymmetry models require \n
            both cases and controls\n")
    }
    mteffect <- "as.factor(mt_MaS)"
    modelformula <- "count~-1+" # I think I can remove intercept
  }

  if (is.element("E:M", effects)) { # check gene-environment paper on adding other effects ex E:C
    if (sum(origDat$D == 0) > 0) { # There are controls; include environmental interaction
      origDat$C <- origDat$C * origDat$D # 1 if C=1, D=1; 2 if C=2, D=1; 0 OW
      origDat$M <- origDat$M * origDat$D # 1 if C=1, D=1; 2 if C=2, D=1; 0 OW
      # (Note that the E:M term in model will be from crossing this
      # variable with E=1, which is exactly what is needed.
      if (mtmodel == "HWE") { # Include main effect of E to give different intercept for HW+E case
        ## For imprinting, If:E
        modeleffects <- c(mteffect, paste0(mteffect, ":E"), effects, "E", "D", "E:D")
      } else {
        ## Look into if need E as well here
        modeleffects <- c(mteffect, paste0(mteffect, ":E"), effects, "D", "E:D")
      }
    } else { # No controls

      if (mtmodel == "HWE") { # Include main effect of E to give different intercept for HW+E case
        modeleffects <- c(mteffect, paste0(mteffect, ":E"), effects, "E")
      } else {
        modeleffects <- c(mteffect, paste0(mteffect, ":E"), effects)
      }
    }
  } else { # No E:M effect

    if (sum(origDat$D == 0) > 0) {
      origDat$C <- origDat$C * origDat$D # 1 if C=1, D=1; 2 if C=2, D=2; 0 OW
      origDat$M <- origDat$M * origDat$D # 1 if C=1, D=1; 2 if C=2, D=2; 0 O
      modeleffects <- c(mteffect, effects, "D")
    } else {
      modeleffects <- c(mteffect, effects)
    }
  }

  linpred <- paste(modeleffects,
                   collapse = "+")
  modelformula <- paste0(modelformula, linpred)

  # Setup results objects
  resVec <- vector(length = length(effects))
  names(resVec) <- effects
  pvalVec <- vector(length = length(effects))
  names(pvalVec) <- effects
  resVecPS <- NULL
  pvalVecPS <- NULL

  # Include test and results under population stratification
  ## try catch for PStest in case glm fails due to low counts
  if (PStest == TRUE) {
    resVecPS <- vector(length = length(effects))
    names(resVecPS) <- effects
    pvalVecPS <- vector(length = length(effects))
    names(pvalVecPS) <- effects

    if (sum(origDat$D == 0) == 0) {
      stop("Can only test for population stratification if there are control trios\n")
    } else {
      PSeffect <- paste0(mteffect, ":D")
      modelformula.PS <- paste(modelformula, PSeffect, sep = "+")
    }
  }

  # Run model and save results

  # EM for Imprinting, using same stopping criteria as Haplin...
  if(any(c("Im", "If") %in% effects)){
    counter <- 0
    if(EM.diag){
      message(paste0("Initial proportion for maternal inheritance cell = ", Minit))
    }
    prev_coeffs <- -1
    prev_dev <- -1
    repeat{
      counter <- counter + 1
      res <- suppressWarnings(
        glm(as.formula(modelformula), data = origDat, family = poisson(), x = TRUE)
      )
      class(res) <- c("TriLLIEM", "glm", "lm")

      Imhat <- ifelse(is.na(exp(res$coefficients["Im"])),
                      1,
                      exp(res$coefficients["Im"])
                      )
      Ifhat <- ifelse(is.na(exp(res$coefficients["If"])),
                      1,
                      exp(res$coefficients["If"])
                      )

      if(EM.diag){
        message(paste0("\nIteration ", counter, " of EM algorithm.
                        \nIm hat = ", Imhat,
                       "\nIf hat = ", Ifhat,
                       "\nProportion for maternal inheritance cell = ", Imhat/(Ifhat + Imhat)))
      }

      origDat <- add_PoO_data(dat,
                              Mprop = Imhat/(Ifhat + Imhat)
                              ) %>%
        dplyr::mutate(HWgeno = M + F)
      ## Use a proper deviance function for imprinting
      ## Check Haplin LogLik code
      if(abs(deviance(res) - prev_dev) < 2e-006 &&
         max(abs(coef(res) - prev_coeffs)) < 1e-006){
        break
      }
      else if(counter == max.iter){
        stop("Max iterations reached without convergence.")
        break
      }
      prev_coeffs <- coef(res)
      prev_dev <- deviance(res)
    }

    filled_inds <- which(dat$type == 9)
    res$aic <- -2 * (sum(dpois(x = res$y[-c(filled_inds, filled_inds+1)], lambda = res$fitted.values[-c(filled_inds, filled_inds+1)], log = TRUE),
                         dpois(x = dat$count[filled_inds], lambda = (Imhat + Ifhat) * exp(res$coefficients[["as.factor(mt_MS)4"]] + res$coefficients[["C"]] + res$coefficients[["M"]]), log = TRUE)
                         )
                     ) +
      2 * res$rank
    ## saturated model set lambda = counts (perfect)
    ## Instead of using subclass, return a Trilliem object only that has the glm returned as a part of it (with incorrect), and the relevant correct info, avoids the redundant functions that return nonsense
    ## Add EM coefficients
    ## EM works with 16-rows due to the "proof"
    ## Test out emax.glm
    ## Compare EM with emax.glm performance
    ## Data abstraction and problem solving with cpp Chpt 1
  } else {
    res <- glm(as.formula(modelformula), data = dat, family = poisson())
    class(res) <- c("TriLLIEM", "glm", "lm")
  }

  # R is not consistent about how interaction is specified. Even though it
  # is fit as E:M, sometimes R flips it to M:E in the output of results.
  # if (is.element("M:E", names(coef(res)))) {
  #   effects[effects == "E:M"] <- "M:E"
  # }
  #
  # for (j in 1:length(effects)) {
  #   resVec[j] <- exp(summary(res)$coef[effects[j], 1])
  #   pvalVec[j] <- summary(res)$coef[effects[j], 4]
  # }
  #
  # test.res <- NULL
  # if (PStest == TRUE) {
  #   res.PS <- glm(as.formula(modelformula.PS), data = origDat, family = poisson())
  #   test.res <- anova(res, res.PS, test = "LRT")
  #
  #   for (j in 1:length(effects)) {
  #     resVecPS[j] <- exp(summary(res.PS)$coef[effects[j], 1])
  #     pvalVecPS[j] <- summary(res.PS)$coef[effects[j], 4]
  #   }
  # }

  return(res)
}
