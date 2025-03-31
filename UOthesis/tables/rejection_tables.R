library(tidyverse)
library(knitr); options(knitr.kable.NA = '')
library(kableExtra)

conditions <- expand.grid(
  R = 1,
  S = 1,
  V = 1,
  mtCoef = 1,
  Im = 1,
  If = 1,
  ntrios = 2000,
  includeE = TRUE,
  Einteraction = c("M", "Im", "If"),
  includeControl = c(TRUE, FALSE),
  nControl = 1000,
  stringsAsFactors = FALSE
) %>%
  rowid_to_column("rowid") %>%
  mutate(rowid = (rowid - 1) %% (nrow(.)/2) + 1)
trill_PS_results <- readRDS("../data/TriLLIEM_res_PS.RDS")
source("../utils/get_reject.R")

reject_table <- get_trill_reject_table(trill_PS_results, conditions)

sink(paste0("reject_PS_c.tex"))
cat(
  kbl(reject_table %>% filter(strat == "nostrat", includeControl) %>% select(-strat, -includeControl, -rowid, -"E:C"),
      format = "latex",
      align = "l",
      col.names = c("$\\CC$", "$\\MM$", "$\\EE:\\MM$", "$\\EE:I_{\\MM}$", "$\\EE:I_{\\FF}$", "$I_{\\MM}$", "$I_{\\FF}$", "Mating type"),
      caption = "Rejection rates at $\\alpha < 0.05$ for various models in our simulation study on
                    population stratification, where the data includes both cases and controls.
                    Given rates are the average across our 2,000 simulations.",
      digits = 5,
      label = "reject_PS_c",
      escape = FALSE,
      booktabs = TRUE,
      linesep = c("", "", "\\addlinespace"),
      position = "H") %>%
    kable_styling(latex_options="scale_down")#%>%
  #column_spec(column = c(2, 5, 6), width = c("1.5in", "0.8in", "0.6in"))
)
sink()

sink(paste0("reject_PS_noc.tex"))
cat(
  kbl(reject_table %>% filter(strat == "nostrat", !includeControl) %>% select(-strat, -includeControl, -rowid, -"E:C"),
      format = "latex",
      align = "l",
      col.names = c("$\\CC$", "$\\MM$", "$\\EE:\\MM$", "$\\EE:I_{\\MM}$", "$\\EE:I_{\\FF}$", "$I_{\\MM}$", "$I_{\\FF}$", "Mating type"),
      caption = "Rejection rates at $\\alpha < 0.05$ for various models in our simulation study on
                    population stratification, where the data only has cases.
                    Given rates are the average over 2,000 simulations.",
      digits = 5,
      label = "reject_PS_noc",
      escape = FALSE,
      booktabs = TRUE,
      linesep = c("", "", "\\addlinespace"),
      position = "H") %>%
    kable_styling(latex_options="scale_down")#%>%
  #column_spec(column = c(2, 5, 6), width = c("1.5in", "0.8in", "0.6in"))
)
sink()
