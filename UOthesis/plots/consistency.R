library(TriLLIEM)
library(tidyverse);theme_set(theme_bw())

dat <- readRDS("../data/consistencydat.RDS")

paramsets <-
  list(
    c("C", "M"),
    c("M", "Im", "If"),
    c("C", "M", "Im")
  )

paramsets_nohap <-
  list(
  )

paramsets_strat <-
  list(
    c("C", "M", "E:M"),
    c("C", "M", "Im", "E:Im")
  )

runall <- function(effects, mtmodel, dat, hapon, group, Einter) {
  trill <- TriLLIEM(effects = effects, mtmodel = mtmodel, dat = dat$dat4R, includeD = TRUE, Estrat = TRUE, includeE = any(grep("E", effects)), Einteraction = Einter) %>% summary() %>% coef()
  emim <- runEMIM(effects = effects, mtmodel = mtmodel, peddat = dat$dat4EMIM, includeE = any(grep("E", effects)), Einteraction = Einter) %>% TriLLIEM:::summ_emim()
  if (hapon) {
    hap <- runHaplin(effects = effects, dat = dat$dat4haplin, includeD = TRUE, PoO = "Im" %in% effects, includeE = any(grep("E", effects))) %>% TriLLIEM:::summ_haplin(PoO = "Im" %in% effects, includeE = any(grep("E", effects)))
  }

  if (any(grep("E", effects))){
    effects <- effects[grep("E", effects)]
  }

  emimeff <- emim$effects[names(emim$effects) %in% effects]
  emimp <- emim$pvals[names(emim$pvals) %in% effects]

  if (hapon) {
    hapeff <- hap$effects[names(hap$effects) %in% effects]
    happ <- hap$pvals[names(hap$pvals) %in% effects]
  } else {
    hapeff <- rep(NA, length(effects))
    happ <- rep(NA, length(effects))
  }

  point_ests <-
    data.frame(trill = exp(trill[effects, 1]), EMIM = emimeff, Haplin = hapeff) %>%
    rownames_to_column(var = "Effect") %>%
    mutate(cat = "est",
           mt = mtmodel,
           group = group)
  pvals <-
    data.frame(trill = trill[effects, 4], EMIM = emimp, Haplin = happ) %>%
    rownames_to_column(var = "Effect") %>%
    mutate(cat = "pval",
           mt = mtmodel,
           group = group)
  return(rbind(point_ests, pvals))
}

results <-
  lapply(
    X = c("HWE", "MS", "MaS"),
    FUN = function(x){
      mapply(
        FUN = runall,
        effects = list(
          c("C", "M"),
          c("M", "Im", "If"),
          c("C", "M", "Im"),
          c("C", "M", "E:M"),
          c("C", "M", "Im", "E:Im")
        ),
        mtmodel = list(x),
        dat = list(dat),
        Einter = list("M", "M", "M", "M", "Im"),
        hap = if (x == "HWE") list(TRUE, TRUE, FALSE, TRUE, FALSE) else FALSE,
        group = c(1:4, 4),
        SIMPLIFY = FALSE
      ) %>%
        bind_rows()
    }
  ) %>%
  bind_rows() %>%
  mutate(mt = factor(mt, levels = c("HWE", "MS", "MaS")))

pdf(file = "consistency-est.pdf", width = 6.5, height = 5)
ggplot(results %>% filter(cat == "est") %>% pivot_longer(cols = c("EMIM", "Haplin"), names_to = "Method", values_to = "value")) +
  geom_point(aes(x = trill, y = value, shape = Method, colour = Effect), alpha = 0.5) +
  geom_abline(col = "red", alpha = 0.25) +
  facet_wrap(mt ~ group, scales = "free") +
  xlab("TriLLIEM estimated risk") +
  ylab("Method estimated risk") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  theme(axis.text.x = element_text(angle = 70, hjust=1))
dev.off()

pdf(file = "consistency-p.pdf", width = 6.5, height = 5)
ggplot(results %>% filter(cat == "pval") %>% pivot_longer(cols = c("EMIM", "Haplin"), names_to = "Method", values_to = "value")) +
  geom_point(aes(x = trill, y = value, shape = Method, colour = Effect), alpha = 0.5) +
  geom_abline(col = "red", alpha = 0.25) +
  facet_wrap(mt ~ group, scales = "free") +
  xlab("TriLLIEM estimated p-value") +
  ylab("Method estimated p-value") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  theme(axis.text.x = element_text(angle = 70, hjust=1))
dev.off()

