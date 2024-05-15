add_PoO_data <- function(dat, Mprop) {
  # Portion of model equation and offset depends on mating type model
  heteroInds <- with(dat, which(type == 9))
  M.count <- Mprop * dat$count[heteroInds]
  PoO_dat <- dat %>%
    dplyr::left_join(PoO_df, by = c("M", "F", "C")) %>%
    dplyr::mutate(count = replace(count, is.na(patOrg), M.count),
                  matOrg = replace(matOrg, is.na(matOrg), 1),
                  patOrg = replace(patOrg, is.na(patOrg), 0)) %>%
    dplyr::add_row(dat %>%
                     dplyr::filter(dplyr::row_number() %in% heteroInds) %>%
                     dplyr::mutate(count = dat$count[heteroInds] - M.count,
                                   matOrg = 0,
                                   patOrg = 1)) %>%
    dplyr::arrange(desc(D), desc(E), type, desc(matOrg)) %>%
    dplyr::mutate(typeOrig = rep(1:16, dplyr::n()/16),
                  Im = matOrg * D,
                  If = patOrg * D) %>%
    dplyr::relocate(typeOrig)
  return(PoO_dat)
}

add_PoO_data_15 <- function(dat) {
  # Portion of model equation and offset depends on mating type model
  heteroInds <- with(dat, which(type == 9))
  PoO_dat <- dat %>%
    dplyr::left_join(PoO_df, by = c("M", "F", "C")) %>%
    dplyr::mutate(matOrg = replace(matOrg, is.na(matOrg), 0),
                  patOrg = replace(patOrg, is.na(patOrg), 0)) %>%
    dplyr::mutate(Im = matOrg * D,
                  If = patOrg * D) %>%
    dplyr::arrange(desc(D), desc(E), type, desc(matOrg))
  return(PoO_dat)
}

## elements of covariance matrix function, a = row index, b = col index
## mu = complete poisson means
## L must be a matrix of 0s and 1s with no more than one 1 per column
cov_trill <- function(y, L, a, b, mu){
  if (a == b) {
    if (1 %in% L[,a]) {
      r_a <- which(L[,a] == 1)
      return(y[[r_a]] *
               mu[[a]] / (t(L[r_a,]) %*% mu) *
               (1 - mu[[a]] / (t(L[r_a,]) %*% mu))
      )
    } else {
      return(mu[[a]])
    }
  } else {
    r_a <- which(L[,a] == 1)
    r_b <- which(L[,b] == 1)
    if(r_a == r_b){
      return(- y[[r_a]] *
               mu[[a]] / (t(L[r_a,]) %*% mu) *
               mu[[b]] / (t(L[r_b,]) %*% mu)
      )
    } else {
      return(0)
    }
  }
}

summ_trill <- function(res){
  resVec <- exp(coef(res))
  if (any(c("Im", "If") %in% as.character(attributes(res$terms)$variables))) {
    seVec <- TriLLIEM_SE(res)
  } else {
    seVec <- coef(summary(res))[,"Std. Error"]
  }
  pvalVec <- 2 * pnorm(abs(log(resVec) / seVec), lower = F)
  return(list(effects = resVec, se = seVec, pvals = pvalVec))
}

summ_haplin <- function(res, PoO = FALSE){
  resVec <- c()
  pvalVec <- c()
  if(PoO){
    resVec["M"] <- res$RRm.est
    pvalVec["M"] <- res$RRm.p.value
    resVec["RRcm"] <- res$RRcm.est
    pvalVec["RRcm"] <- res$RRcm.p.value
    resVec["RRcf"] <- res$RRcf.est
    pvalVec["RRcf"] <- res$RRcf.p.value
  }
  else{
    resVec["C"] <- res$RR.est.
    pvalVec["C"] <- res$RR.p.value
    resVec["M"] <- res$RRm.est.
    pvalVec["M"] <- res$RRm.p.value
  }
  return(list(effects = resVec, pvals = pvalVec))
}

summ_emim <- function(res){
  resVec <- c()
  sdVec <- c()
  pvalVec <- c()
  resVec["C"] <- exp(res$lnR1)
  sdVec["C"] <- res$sd_lnR1
  pvalVec["C"] <- 2 * pnorm(abs(res$lnR1 / res$sd_lnR1), lower = F)
  resVec["M"] <- exp(res$lnS1)
  sdVec["M"] <- res$sd_lnS1
  pvalVec["M"] <- 2 * pnorm(abs(res$lnS1 / res$sd_lnS1), lower = F)
  resVec["Im"] <- exp(res$lnIm)
  sdVec["Im"] <- res$sd_lnIm
  pvalVec["Im"] <- 2 * pnorm(abs(res$lnIm / res$sd_lnIm), lower = F)
  resVec["Ip"] <- exp(res$lnIp)
  sdVec["Ip"] <- res$sd_lnIp
  pvalVec["Ip"] <- 2 * pnorm(abs(res$lnIp / res$sd_lnIp), lower = F)
  return(list(effects = resVec, se = sdVec, pvals = pvalVec))
}
