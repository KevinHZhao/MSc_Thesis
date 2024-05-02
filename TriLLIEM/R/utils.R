add_PoO_data <- function(dat, Mprop) {
  # Portion of model equation and offset depends on mating type model
  heteroInds <- with(dat, which((M == 1) & (F == 1) & (C == 1)))
  PoO_dat <- dat %>%
    dplyr::left_join(PoO_df, by = c("M", "F", "C")) %>%
    dplyr::mutate(count = ifelse(is.na(patOrg), ceiling(Mprop * count), count),
                  patOrg = ifelse(is.na(patOrg), 0, patOrg),
                  matOrg = ifelse(is.na(matOrg), 1, matOrg)) %>%
    dplyr::add_row(dat %>%
                     dplyr::filter(dplyr::row_number() == heteroInds) %>%
                     dplyr::mutate(count = floor((1 - Mprop) * count),
                                   patOrg = 1,
                                   matOrg = 0)) %>%
    dplyr::arrange(desc(D), desc(E), type, desc(matOrg)) %>%
    dplyr::mutate(typeOrig = rep(1:16, dplyr::n()/16),
                  If = patOrg * D,
                  Im = matOrg * D) %>%
    dplyr::relocate(typeOrig)
  return(PoO_dat)
}

summ_trill <- function(res, effects){
  resVec <- c()
  sdVec <- c()
  pvalVec <- c()
  for (j in effects) {
    resVec[j] <- exp(summary(res)$coef[j, 1])
    sdVec[j] <- summary(res)$coef[j, 2]
    pvalVec[j] <- summary(res)$coef[j, 4]
  }
  return(list(effects = resVec, se = sdVec, pvals = pvalVec))
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
  sdVec["Im"] <- res$sd_lnIp
  pvalVec["Ip"] <- 2 * pnorm(abs(res$lnIp / res$sd_lnIp), lower = F)
  return(list(effects = resVec, se = sdVec, pvals = pvalVec))
}
