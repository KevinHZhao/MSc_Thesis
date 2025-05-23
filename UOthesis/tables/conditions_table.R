library(tidyverse)
library(knitr); options(knitr.kable.NA = '')
library(kableExtra)

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
  select(-includeE) %>%
  rowid_to_column("rowid") %>%
  mutate(rowid = (rowid - 1) %% (nrow(.)/2) + 1)

sink(paste0("conditions_table.tex"))
cat(
  kbl(conditions %>%
        filter(!includeControl) %>%
        mutate(Einteraction = ifelse(Einteraction == "If", "$I_{\\FF}$", "$I_{\\MM}$")) %>%
        select(rowid, R, S, V, mtCoef, Im, If, propE, Einteraction),
      format = "latex",
      align = "l",
      col.names = c("ID", "$R$", "$S$", "$V$", "$C_i$", "$I_{\\MM}$", "$I_{\\FF}$", "$\\varepsilon$", "Interaction"),
      caption = "Conditions used for simulating data in our simulation study.",
      digits = 5,
      label = "conditions",
      escape = FALSE,
      booktabs = TRUE,
      longtable = TRUE,
      linesep = "",
      position = "H") %>%
    kable_styling(latex_options=c("hold_position"))#%>%
  #column_spec(column = c(2, 5, 6), width = c("1.5in", "0.8in", "0.6in"))
)
sink()
