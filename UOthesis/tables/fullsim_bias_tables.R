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
  includeE = TRUE,
  Einteraction = c("Im", "If"),
  ntrios = 2000,
  propE = c(0.3, 0.5),
  includeControl = c(FALSE, TRUE),
  nControl = 1000,
  stringsAsFactors = FALSE
) %>%
  filter(
    (R == 1) + (S == 1) + (V == 1) + (Im == 1) + (If == 1) >= 4,
    !(If != 1 & (Im != 1 | Einteraction == "Im")),
    !(Im != 1 & Einteraction == "If")
  ) %>%
  rowid_to_column("rowid") %>%
  mutate(rowid = (rowid - 1) %% (nrow(.)/2) + 1)

source("../utils/get_bias.R")

bias_table <- get_trill_table(trill_res, conditions)

sink("bias_eif_noc.tex")
cat(
  kbl(bias_table %>%
        left_join(conditions, by = c("includeControl", "rowid"), suffix = c("", ".T")) %>%
        filter(!includeControl, Einteraction == "If", propE == 0.3, mtCoef != 1.15) %>%
        arrange(rowid) %>%
        mutate(effects = get_true_effects(conditions = conditions, id = rowid),
               strat = case_when(strat == "noE" ~ "no $\\EE$",
                                 strat == "strat" ~ "S",
                                 strat == "nostrat" ~ "NS")) %>%
        select(rowid, C, M, "E:If", If, mt, strat, effects, mtCoef),
      format = "latex",
      align = "l",
      booktabs = TRUE,
      longtable = TRUE,
      col.names = c("ID", "$\\CC$", "$\\MM$", "$\\EE:I_{\\FF}$", "$I_{\\FF}$", "MT", "$\\EE$ approach", "True effect", "MaS"),
      caption = "Relative biases for simulated data sets without controls where $\\EE = 0.3$, where the interaction parameter is $I_{\\FF}$.  S and NS stand for the stratified and non-stratified approaches for gene-environment interactions, respectively.  The MaS column represents the value of the mating asymmetry coefficients used during simulation.  To see the scenario conditions, please use the ID column to refer to the appropriate row of \\cref{tab:conditions}.  Results from data simulated with MaS coefficients of $1.15$ are omitted for brevity.",
      digits = 5,
      label = "bias_eif_noc",
      escape = FALSE,
      linesep = c("", "", "\\addlinespace"),
      caption.short = "Average relative biases for $I_{\\FF}$ models without controls.",
      position = "H") %>%
    kable_styling(latex_options=c("repeat_header", "hold_position"),
                  repeat_header_method = "replace")#%>%
  #column_spec(column = c(2, 5, 6), width = c("1.5in", "0.8in", "0.6in"))
)
sink()

sink("bias_eim_noc.tex")
cat(
  kbl(bias_table %>%
        left_join(conditions, by = c("includeControl", "rowid"), suffix = c("", ".T")) %>%
        filter(!includeControl, Einteraction == "Im", propE == 0.3, mtCoef != 1.15) %>%
        arrange(rowid) %>%
        mutate(effects = get_true_effects(conditions = conditions, id = rowid),
               strat = case_when(strat == "noE" ~ "no $\\EE$",
                                 strat == "strat" ~ "S",
                                 strat == "nostrat" ~ "NS")) %>%
        select(rowid, C, M, "E:Im", Im, mt, strat, effects, mtCoef),
      format = "latex",
      align = "l",
      col.names = c("ID", "$\\CC$", "$\\MM$", "$\\EE:I_{\\MM}$", "$I_{\\MM}$", "MT", "$\\EE$ approach", "True effect", "MaS"),
      caption = "Relative biases for simulated data sets without controls where $\\EE = 0.3$, where the interaction parameter is $I_{\\MM}$.  S and NS stand for the stratified and non-stratified approaches for gene-environment interactions, respectively.  The MaS column represents the value of the mating asymmetry coefficients used during simulation.  To see the conditions scenario conditions, please use the ID column to refer to the appropriate row of \\cref{tab:conditions}.  Results from data simulated with MaS coefficients of $1.15$ are omitted for brevity.",
      digits = 5,
      label = "bias_eim_noc",
      escape = FALSE,
      booktabs = TRUE,
      longtable = TRUE,
      linesep = c("", "", "\\addlinespace"),
      caption.short = "Average relative biases for $I_{\\MM}$ models without controls.",
      position = "H") %>%
    kable_styling(latex_options=c("repeat_header", "hold_position"),
                  repeat_header_method = "replace")#%>%
  #column_spec(column = c(2, 5, 6), width = c("1.5in", "0.8in", "0.6in"))
)
sink()

sink("bias_eif_c.tex")
cat(
  kbl(bias_table %>%
        left_join(conditions, by = c("includeControl", "rowid"), suffix = c("", ".T")) %>%
        filter(includeControl, Einteraction == "If", propE == 0.3, mtCoef != 1.15) %>%
        arrange(rowid) %>%
        mutate(effects = get_true_effects(conditions = conditions, id = rowid),
               strat = case_when(strat == "noE" ~ "no $\\EE$",
                                 strat == "strat" ~ "S",
                                 strat == "nostrat" ~ "NS")) %>%
        select(rowid, C, M, "E:If", If, mt, strat, effects, mtCoef),
      format = "latex",
      align = "l",
      col.names = c("ID", "$\\CC$", "$\\MM$", "$\\EE:I_{\\FF}$", "$I_{\\FF}$", "MT model", "$\\EE$ approach", "True effect", "MaS"),
      caption = "Relative biases for simulated data sets with controls where $\\EE = 0.3$, where the environmental interaction paramter is $I_{\\FF}$.  S and NS stand for the stratified and non-stratified approaches for gene-environment interactions, respectively.  The MaS column represents the value of the mating asymmetry coefficients used during simulation.  To see the scenario conditions, please use the ID column to refer to the appropriate row of \\cref{tab:conditions}.  Results from data simulated with MaS coefficients of $1.15$ are omitted for brevity.",
      digits = 5,
      label = "bias_eif_c",
      escape = FALSE,
      booktabs = TRUE,
      longtable = TRUE,
      linesep = c("", "", "\\addlinespace"),
      caption.short = "Average relative biases for $I_{\\F}$ models with controls.",
      position = "H") %>%
    kable_styling(latex_options=c("repeat_header", "hold_position"),
                  repeat_header_method = "replace")#%>%
  #column_spec(column = c(2, 5, 6), width = c("1.5in", "0.8in", "0.6in"))
)
sink()

sink("bias_eim_c.tex")
cat(
  kbl(bias_table %>%
        left_join(conditions, by = c("includeControl", "rowid"), suffix = c("", ".T")) %>%
        filter(includeControl, Einteraction == "Im", propE == 0.3, mtCoef != 1.15) %>%
        arrange(rowid) %>%
        mutate(effects = get_true_effects(conditions = conditions, id = rowid),
               strat = case_when(strat == "noE" ~ "no $\\EE$",
                                 strat == "strat" ~ "S",
                                 strat == "nostrat" ~ "NS")) %>%
        select(rowid, C, M, "E:Im", Im, mt, strat, effects, mtCoef),
      format = "latex",
      align = "l",
      col.names = c("ID", "$\\CC$", "$\\MM$", "$\\EE:I_{\\MM}$", "$I_{\\MM}$", "MT model", "$\\EE$ approach", "True effect", "MaS"),
      caption = "Relative biases for simulated data sets with controls where $\\EE = 0.3$, where the environmental interaction paramter is $I_{\\MM}$.  S and NS stand for the stratified and non-stratified approaches for gene-environment interactions, respectively.  The MaS column represents the value of the mating asymmetry coefficients used during simulation.  To see the scenario conditions, please use the ID column to refer to the appropriate row of \\cref{tab:conditions}.  Results from data simulated with MaS coefficients of $1.15$ are omitted for brevity.",
      digits = 5,
      label = "bias_eim_c",
      escape = FALSE,
      booktabs = TRUE,
      longtable = TRUE,
      linesep = c("", "", "\\addlinespace"),
      caption.short = "Average relative biases for $I_{\\MM}$ models with controls.",
      position = "H") %>%
    kable_styling(latex_options=c("repeat_header", "hold_position"),
                  repeat_header_method = "replace")#%>%
  #column_spec(column = c(2, 5, 6), width = c("1.5in", "0.8in", "0.6in"))
)
sink()

write.csv2(
  bias_table %>%
    left_join(conditions, by = c("includeControl", "rowid"), suffix = c("", ".T")) %>%
    arrange(rowid) %>%
    mutate(effects = get_true_effects(conditions = conditions, id = rowid),
           strat = case_when(strat == "noE" ~ "no $\\EE$",
                             strat == "strat" ~ "S",
                             strat == "nostrat" ~ "NS")),
  "fulltables/full_bias_table.csv"
)

