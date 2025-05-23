library(tidyverse)
library(knitr); options(knitr.kable.NA = '')
library(kableExtra)

source("../utils/get_bias.R")

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

bias_table <- get_trill_table(trill_PS_results, conditions) %>% mutate(strat = rep(c("strat", "nostrat"), 15))

sink(paste0("bias_PS_c.tex"))
cat(
  kbl(bias_table %>% filter(strat == "nostrat", includeControl) %>% select(-strat, -includeControl, -rowid, -"E:C"),
      format = "latex",
      align = "l",
      col.names = c("$\\CC$", "$\\MM$", "$\\EE:\\MM$", "$\\EE:I_{\\MM}$", "$\\EE:I_{\\FF}$", "$I_{\\MM}$", "$I_{\\FF}$", "Mating type"),
      caption = "Average relative biases across $2,000$ datasets for various models in our simulation study on
                    population stratification, where the data has both cases and controls.",
      digits = 5,
      label = "bias_PS_c",
      escape = FALSE,
      booktabs = TRUE,
      linesep = c("", "", "\\addlinespace"),
      caption.short = "Average relative biases for various gene-environment interaction models for hybrid data with population stratification.",
      position = "H") %>%
    kable_styling(latex_options="scale_down")#%>%
  #column_spec(column = c(2, 5, 6), width = c("1.5in", "0.8in", "0.6in"))
)
sink()

sink(paste0("bias_PS_noc.tex"))
cat(
  kbl(bias_table %>% filter(strat == "nostrat", !includeControl) %>% select(-strat, -includeControl, -rowid, -"E:C"),
      format = "latex",
      align = "l",
      col.names = c("$\\CC$", "$\\MM$", "$\\EE:\\MM$", "$\\EE:I_{\\MM}$", "$\\EE:I_{\\FF}$", "$I_{\\MM}$", "$I_{\\FF}$", "Mating type"),
      caption = "Average relative biases across $2,000$ datasets for various models in our simulation study on
                    population stratification, where the data only has cases.",
      digits = 5,
      label = "bias_PS_noc",
      escape = FALSE,
      booktabs = TRUE,
      linesep = c("", "", "\\addlinespace"),
      caption.short = "Average relative biases for various gene-environment interaction models for case-triad data with population stratification.",
      position = "H") %>%
    kable_styling(latex_options="scale_down")#%>%
  #column_spec(column = c(2, 5, 6), width = c("1.5in", "0.8in", "0.6in"))
)
sink()

## In methods, can justify 2000 repetitions because of maximum confidence interval width around binomial
## Get emim hap results for PS
