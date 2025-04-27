library(tidyverse); theme_set(theme_bw())
library(knitr); options(knitr.kable.NA = '')
library(kableExtra)
library(ggpattern)

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
source("../utils/get_bias.R")

bias_table <- get_trill_table(trill_res, conditions)
reject_table <- get_trill_reject_table(trill_res, conditions) %>%
  left_join(conditions, by = c("includeControl", "rowid"), suffix = c("", ".T"))

pdf("stratcomp_Im.pdf", width = 14, height = 6)
reject_table %>%
  filter(mtCoef == 1, mt == "HWE", R == 1, S == 1, Im.T == 1, If.T == 1, includeE, Einteraction == "Im")%>%
  arrange(rowid) %>%
  select(propE, includeControl, rowid, "E:Im", V, strat) %>%
  left_join(bias_table %>% filter(mt == "HWE") %>% select(rowid, includeControl, "E:Im", strat), by = c("rowid", "includeControl", "strat"), suffix = c(".Rejection rate", ".Relative bias")) %>%
  pivot_longer(cols = starts_with("E:Im."), values_to = "value", names_to = "type", names_prefix = "E:Im.") %>%
  mutate(propE = factor(propE , levels = c(0.3, 0.5), labels = c("30% exposures", "50% exposures")),
         V = factor(V, levels = c(1, 1.5, 1.6)),
         includeControl = factor(includeControl, levels = c(FALSE, TRUE), labels = c("Case-triad data only", "Case-triad & control-triad data"))) %>%
  ggplot() +
  geom_bar_pattern(
    aes(x = V, y = value, fill = includeControl, pattern = strat),
    stat = "identity",
    position = position_dodge2(),
    colour = "black",
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.025,
    pattern_key_scale_factor = 0.6) +
  facet_grid(cols = vars(propE), rows = vars(type), scales = "free_y", switch = "y") +
  scale_fill_brewer(palette="Dark2") +
  scale_pattern_manual(values = c(nostrat = "none", strat = "stripe"), guide = "none") +
  xlab("True maternal imprinting by environment risk") +
  ylab("Rejection rate") +
  labs(fill = "Data set composition") +
  theme(
    legend.position = c(.175, .985),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    axis.title.y = element_blank(),
    strip.placement = "outside"
  ) +
  guides(fill = guide_legend(override.aes = list(pattern = c("none", "none")))) +
  geom_hline(data = . %>% filter(type == "Rejection rate"), aes(yintercept = 0.05), colour = "red", alpha = 0.75)
dev.off()

pdf("stratcomp_If.pdf", width = 14, height = 6)
reject_table %>%
  filter(mtCoef == 1, mt == "HWE", R == 1, S == 1, Im.T == 1, If.T == 1, includeE, Einteraction == "If")%>%
  arrange(rowid) %>%
  select(propE, includeControl, rowid, "E:If", V, strat) %>%
  left_join(bias_table %>% filter(mt == "HWE") %>% select(rowid, includeControl, "E:If", strat), by = c("rowid", "includeControl", "strat"), suffix = c(".Rejection rate", ".Relative bias")) %>%
  pivot_longer(cols = starts_with("E:If."), values_to = "value", names_to = "type", names_prefix = "E:If.") %>%
  mutate(propE = factor(propE , levels = c(0.3, 0.5), labels = c("30% exposures", "50% exposures")),
         V = factor(V, levels = c(1, 1.5, 1.6)),
         includeControl = factor(includeControl, levels = c(FALSE, TRUE), labels = c("Case-triad data only", "Case-triad & control-triad data"))) %>%
  ggplot() +
  geom_bar_pattern(
    aes(x = V, y = value, fill = includeControl, pattern = strat),
    stat = "identity",
    position = position_dodge2(),
    colour = "black",
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.025,
    pattern_key_scale_factor = 0.6) +
  facet_grid(cols = vars(propE), rows = vars(type), scales = "free_y", switch = "y") +
  scale_fill_brewer(palette="Dark2") +
  scale_pattern_manual(values = c(nostrat = "none", strat = "stripe"), guide = "none") +
  xlab("True paternal imprinting by environment risk") +
  ylab("Rejection rate") +
  labs(fill = "Data set composition") +
  theme(
    legend.position = c(.175, .985),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    axis.title.y = element_blank(),
    strip.placement = "outside"
  ) +
  guides(fill = guide_legend(override.aes = list(pattern = c("none", "none")))) +
  geom_hline(data = . %>% filter(type == "Rejection rate"), aes(yintercept = 0.05), colour = "red", alpha = 0.75)
dev.off()

pdf("stratcomp_C.pdf", width = 14, height = 6)
reject_table %>%
  filter(mtCoef == 1, mt == "HWE", S == 1, V == 1, Im.T == 1, If.T == 1, Einteraction == "Im", strat == "nostrat", includeE)%>%
  arrange(rowid) %>%
  select(propE, includeE, includeControl, rowid, C, R) %>%
  left_join(bias_table %>% filter(strat == "nostrat", mt == "HWE") %>% select(rowid, includeControl, C), by = c("rowid", "includeControl"), suffix = c(".Rejection rate", ".Relative bias")) %>%
  pivot_longer(cols = starts_with("C."), values_to = "value", names_to = "type", names_prefix = "C.") %>%
  mutate(propE = factor(propE * includeE, levels = c(0.3, 0.5), labels = c("30% exposures", "50% exposures")),
         R = factor(R, levels = c(1, 1.15)),
         includeControl = factor(includeControl, levels = c(FALSE, TRUE), labels = c("Case-triad data only", "Case-triad & control-triad data"))) %>%
  ggplot() +
  geom_bar(
    aes(x = R, y = value, fill = includeControl),
    stat = "identity",
    position = position_dodge2(),
    colour = "black") +
  facet_grid(cols = vars(propE), rows = vars(type), scales = "free_y", switch = "y") +
  scale_fill_brewer(palette="Dark2") +
  xlab("True child effect") +
  ylab("Rejection rate") +
  labs(fill = "Data set composition") +
  theme(
    legend.position = c(.175, .985),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    axis.title.y = element_blank(),
    strip.placement = "outside"
  ) +
  guides(fill = guide_legend(override.aes = list(pattern = c("none", "none")))) +
  geom_hline(data = . %>% filter(type == "Rejection rate"), aes(yintercept = 0.05), colour = "red", alpha = 0.75)
dev.off()
