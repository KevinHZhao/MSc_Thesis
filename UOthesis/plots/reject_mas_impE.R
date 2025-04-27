library(tidyverse); theme_set(theme_bw())
library(knitr); options(knitr.kable.NA = '')
library(kableExtra)

trill_res <- readRDS("../data/TriLLIEM_res.RDS")
#emim_res <- readRDS("../data/EMIM_res.RDS")
#hap_res <- readRDS("../data/Haplin_res.RDS")

conditions <- expand.grid(
  maf = 0.3,
  R = c(1, 1.15),
  S = 1,
  V = c(1, 1.5, 1.6),
  mtCoef = c(0.85, 1, 1.15),
  Im = c(1, 1.2, 1.4),
  If = c(1, 1.2, 1.4),
  includeE = c(FALSE, TRUE),
  Einteraction = c("Im", "If"),
  ntrios = 2000,
  propE = c(0.3, 0.5),
  includeControl = c(FALSE, TRUE),
  nControl = 1000,
  stringsAsFactors = FALSE
) %>%
  filter(
    !(!includeE & !(Einteraction == "Im" & V == 1 & propE == 0.3)),
    (R == 1) + (S == 1) + (V == 1) + (Im == 1) + (If == 1) >= 4,
    !(If != 1 & (Im != 1 | Einteraction == "Im")),
    !(Im != 1 & Einteraction == "If")
  ) %>%
  rowid_to_column("rowid") %>%
  mutate(rowid = (rowid - 1) %% (nrow(.)/2) + 1)

source("../utils/get_reject.R")

reject_table <- get_trill_reject_table(trill_res, conditions) %>%
  left_join(conditions, by = c("includeControl", "rowid"), suffix = c("", ".T"))

pdf("reject_mas_mimpE.pdf", width = 8.5, height = 4.5)
reject_table %>%
  filter(rowid %in% c(13,17,21), strat == "nostrat") %>%
  filter(includeE, Einteraction == "Im") %>%
  arrange(rowid) %>%
  mutate(effects = get_true_effects(conditions = conditions, id = rowid),
         strat = ifelse(strat == "strat", "Y", "N")) %>%
  select(C, M, "E:Im", Im, mt, mtCoef, includeControl) %>%
  pivot_longer(cols = c("C", "M", "E:Im", "Im"), values_to = "rate", names_to = "Effect") %>%
  mutate(mt = factor(mt, levels = c("HWE", "MS", "MaS")),
         mtCoef = factor(mtCoef, levels = c(0.85, 1, 1.15)),
         Effect = factor(Effect, levels = c("C", "M", "E:Im", "Im"), labels = c("C", "M", "E:Im", "Im")),
         includeControl = factor(includeControl, levels = c(FALSE, TRUE), labels = c("Case-triad data only", "Case-triad & control-triad data"))) %>%
  ggplot() +
  geom_bar(aes(x = mt, y = rate, fill = mtCoef), stat = "identity", position = position_dodge2(), colour = "black") +
  scale_fill_brewer(palette="Set1") +
  facet_grid(rows = vars(Effect), cols = vars(includeControl)) +
  geom_hline(yintercept = 0.05, col = "red", alpha = 0.7) +
  xlab("Mating type model") +
  ylab("Rejection rate") +
  labs(fill = "MaS coefficients") +
  theme(
    legend.position = c(.14, .985),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    legend.key.size = unit(0.5, "lines")
  )
dev.off()

pdf("reject_mas_fimpE.pdf", width = 8.5, height = 4.5)
reject_table %>%
  filter(rowid %in% c(31,35,39), strat == "nostrat") %>%
  filter(includeE, Einteraction == "If") %>%
  arrange(rowid) %>%
  mutate(effects = get_true_effects(conditions = conditions, id = rowid),
         strat = ifelse(strat == "strat", "Y", "N")) %>%
  select(C, M, "E:If", If, mt, mtCoef, includeControl) %>%
  pivot_longer(cols = c("C", "M", "E:If", "If"), values_to = "rate", names_to = "Effect") %>%
  mutate(mt = factor(mt, levels = c("HWE", "MS", "MaS")),
         mtCoef = factor(mtCoef, levels = c(0.85, 1, 1.15)),
         Effect = factor(Effect, levels = c("C", "M", "E:If", "If"), labels = c("C", "M", "E:If", "If")),
         includeControl = factor(includeControl, levels = c(FALSE, TRUE), labels = c("Case-triad data only", "Case-triad & control-triad data"))) %>%
  ggplot() +
  geom_bar(aes(x = mt, y = rate, fill = mtCoef), stat = "identity", position = position_dodge2(), colour = "black") +
  scale_fill_brewer(palette="Set1") +
  facet_grid(rows = vars(Effect), cols = vars(includeControl)) +
  geom_hline(yintercept = 0.05, col = "red", alpha = 0.7) +
  xlab("Mating type model") +
  ylab("Rejection rate") +
  labs(fill = "MaS coefficients") +
  theme(
    legend.position = c(.14, .985),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    legend.key.size = unit(0.5, "lines")
  )
dev.off()
