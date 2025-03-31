library(tidyverse)
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

reject_table <- get_trill_reject_table(trill_res, conditions)

sink("reject_noe_noc.tex")
cat(
  kbl(reject_table %>%
        left_join(conditions, by = c("includeControl", "rowid"), suffix = c("", ".T")) %>%
        filter(!includeE, !includeControl) %>%
        arrange(rowid) %>%
        mutate(effects = get_true_effects(conditions = conditions, id = rowid)) %>%
        select(rowid, C, M, Im, mt, effects, mtCoef),
      format = "latex",
      align = "l",
      col.names = c("ID", "$\\CC$", "$\\MM$", "$I_{\\MM}$", "Mating type", "True effects", "mtCoef"),
      caption = "Rejection rates at $\\alpha < 0.05$ for simulated data sets without controls or environmental exposures.  To see the conditions of each simulated data set, please use the ID column to refer to the appropriate row of \\cref{tab:conditions}.",
      digits = 5,
      label = "reject_noe_noc",
      escape = FALSE,
      booktabs = TRUE,
      linesep = c("", "\\addlinespace"),
      position = "H") %>%
    kable_styling(latex_options="scale_down")#%>%
  #column_spec(column = c(2, 5, 6), width = c("1.5in", "0.8in", "0.6in"))
)
sink()

sink("reject_noe_c.tex")
cat(
  kbl(reject_table %>%
        left_join(conditions, by = c("includeControl", "rowid"), suffix = c("", ".T")) %>%
        filter(!includeE, includeControl) %>%
        arrange(rowid) %>%
        mutate(effects = get_true_effects(conditions = conditions, id = rowid)) %>%
        select(rowid, C, M, Im, mt, effects, mtCoef),
      format = "latex",
      align = "l",
      col.names = c("ID", "$\\CC$", "$\\MM$", "$I_{\\MM}$", "Mating type", "True effects", "mtCoef"),
      caption = "Rejection rates at $\\alpha < 0.05$ for simulated data sets with controls but without environmental exposures.  To see the conditions of each simulated data set, please use the ID column to refer to the appropriate row of \\cref{tab:conditions}.",
      digits = 5,
      label = "reject_noe_c",
      escape = FALSE,
      booktabs = TRUE,
      linesep = c("", "", "\\addlinespace"),
      position = "H") %>%
    kable_styling(latex_options="scale_down")#%>%
  #column_spec(column = c(2, 5, 6), width = c("1.5in", "0.8in", "0.6in"))
)
sink()

sink("reject_eif_noc.tex")
cat(
  kbl(reject_table %>%
        left_join(conditions, by = c("includeControl", "rowid"), suffix = c("", ".T")) %>%
        filter(includeE, !includeControl, Einteraction == "If") %>%
        arrange(rowid) %>%
        mutate(effects = get_true_effects(conditions = conditions, id = rowid),
               strat = ifelse(strat == "strat", "Y", "N")) %>%
        select(rowid, C, M, "E:C", "E:M", "E:If", If, mt, strat, effects, mtCoef),
      format = "latex",
      align = "l",
      booktabs = TRUE,
      longtable = TRUE,
      col.names = c("ID", "$\\CC$", "$\\MM$", "$\\EE:\\CC$", "$\\EE:\\MM$", "$\\EE:I_{\\FF}$", "$I_{\\FF}$", "Mating type", "Stratified", "True effects", "mtCoef"),
      caption = "Rejection rates at $\\alpha < 0.05$ for simulated data sets without controls but with environmental exposures, where the interaction parameter is $I_{\\FF}$.  To see the conditions of each simulated data set, please use the ID column to refer to the appropriate row of \\cref{tab:conditions}.",
      digits = 5,
      label = "reject_eif_noc",
      escape = FALSE,
      linesep = c("", "", "", "\\addlinespace"),
      position = "H") %>%
    kable_styling(latex_options=c("repeat_header", "hold_position"))#%>%
  #column_spec(column = c(2, 5, 6), width = c("1.5in", "0.8in", "0.6in"))
)
sink()

sink("reject_eim_noc.tex")
cat(
  kbl(reject_table %>%
        left_join(conditions, by = c("includeControl", "rowid"), suffix = c("", ".T")) %>%
        filter(includeE, !includeControl, Einteraction == "Im") %>%
        arrange(rowid) %>%
        mutate(effects = get_true_effects(conditions = conditions, id = rowid),
               strat = ifelse(strat == "strat", "Y", "N")) %>%
        select(rowid, C, M, "E:C", "E:M", "E:Im", Im, mt, strat, effects, mtCoef),
      format = "latex",
      align = "l",
      col.names = c("ID", "$\\CC$", "$\\MM$", "$\\EE:\\CC$", "$\\EE:\\MM$", "$\\EE:I_{\\MM}$", "$I_{\\MM}$", "Mating type", "Stratified", "True effects", "mtCoef"),
      caption = "Rejection rates at $\\alpha < 0.05$ for simulated data sets without controls but with environmental exposures, where the interaction parameter is $I_{\\MM}$.  To see the conditions of each simulated data set, please use the ID column to refer to the appropriate row of \\cref{tab:conditions}.",
      digits = 5,
      label = "reject_eim_noc",
      escape = FALSE,
      booktabs = TRUE,
      longtable = TRUE,
      linesep = c("", "", "", "\\addlinespace"),
      position = "H") %>%
    kable_styling(latex_options=c("repeat_header", "hold_position"))#%>%
  #column_spec(column = c(2, 5, 6), width = c("1.5in", "0.8in", "0.6in"))
)
sink()

sink("reject_eif_c.tex")
cat(
  kbl(reject_table %>%
        left_join(conditions, by = c("includeControl", "rowid"), suffix = c("", ".T")) %>%
        filter(includeE, includeControl, Einteraction == "If") %>%
        arrange(rowid) %>%
        mutate(effects = get_true_effects(conditions = conditions, id = rowid),
               strat = ifelse(strat == "strat", "Y", "N")) %>%
        select(rowid, C, M, "E:C", "E:M", "E:If", If, mt, strat, effects, mtCoef),
      format = "latex",
      align = "l",
      col.names = c("ID", "$\\CC$", "$\\MM$", "$\\EE:\\CC$", "$\\EE:\\MM$", "$\\EE:I_{\\FF}$", "$I_{\\FF}$", "Mating type", "Stratified", "True effects", "mtCoef"),
      caption = "Rejection rates at $\\alpha < 0.05$ for simulated data sets with controls and environmental exposures, where the environmental interaction paramter is $I_{\\FF}$.  To see the conditions of each simulated data set, please use the ID column to refer to the appropriate row of \\cref{tab:conditions}.",
      digits = 5,
      label = "reject_eif_c",
      escape = FALSE,
      booktabs = TRUE,
      longtable = TRUE,
      linesep = c("", "", "", "", "", "\\addlinespace"),
      position = "H") %>%
    kable_styling(latex_options=c("repeat_header", "hold_position"))#%>%
  #column_spec(column = c(2, 5, 6), width = c("1.5in", "0.8in", "0.6in"))
)
sink()

sink("reject_eim_c.tex")
cat(
  kbl(reject_table %>%
        left_join(conditions, by = c("includeControl", "rowid"), suffix = c("", ".T")) %>%
        filter(includeE, includeControl, Einteraction == "Im") %>%
        arrange(rowid) %>%
        mutate(effects = get_true_effects(conditions = conditions, id = rowid),
               strat = ifelse(strat == "strat", "Y", "N")) %>%
        select(rowid, C, M, "E:C", "E:M", "E:Im", Im, mt, strat, effects, mtCoef),
      format = "latex",
      align = "l",
      col.names = c("ID", "$\\CC$", "$\\MM$", "$\\EE:\\CC$", "$\\EE:\\MM$", "$\\EE:I_{\\MM}$", "$I_{\\MM}$", "Mating type", "Stratified", "True effects", "mtCoef"),
      caption = "Rejection rates at $\\alpha < 0.05$ for simulated data sets with controls and environmental exposures, where the environmental interaction paramter is $I_{\\MM}$.  To see the conditions of each simulated data set, please use the ID column to refer to the appropriate row of \\cref{tab:conditions}.",
      digits = 5,
      label = "reject_eim_c",
      escape = FALSE,
      booktabs = TRUE,
      longtable = TRUE,
      linesep = c("", "", "", "", "", "\\addlinespace"),
      position = "H") %>%
    kable_styling(latex_options=c("repeat_header", "hold_position"))#%>%
  #column_spec(column = c(2, 5, 6), width = c("1.5in", "0.8in", "0.6in"))
)
sink()
