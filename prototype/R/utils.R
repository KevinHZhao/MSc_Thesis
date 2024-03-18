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
                  patOrgInf = patOrg * D,
                  matOrgInf = matOrg * D) %>%
    dplyr::relocate(typeOrig)
  return(PoO_dat)
}
