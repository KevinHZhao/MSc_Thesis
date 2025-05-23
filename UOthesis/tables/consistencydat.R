library(tidyverse)
library(knitr); options(knitr.kable.NA = '')
library(kableExtra)

df <-
  readRDS("../data/consistencydat.RDS")$dat4R %>%
  select(mt_MS, M, F, C, E, D, count)

sink("consistencydat.tex")
df %>%
  pivot_wider(names_from = E, values_from = count, names_prefix = "E=") %>%
  pivot_wider(names_from = D, values_from = c("E=0","E=1"), names_prefix = "D=") %>%
  select(mt_MS, M, F, C, "E=0_D=0", "E=1_D=0", "E=0_D=1", "E=1_D=1") %>%
  kbl(
    format = "latex",
    align = "l",
    col.names = c("Mating type", "$\\MM$", "$\\FF$", "$\\CC$", "$\\EE=0$", "$\\EE=1$", "$\\EE=0$", "$\\EE=1$"),
    caption = "Simulated null data set for comparison between \\textbf{EMIM}, \\textbf{Haplin}, and \\textbf{TriLLIEM}.",
    label = "consistencydat",
    escape = FALSE,
    booktabs = TRUE,
    linesep = "",
    position = "H"
  ) %>%
  add_header_above(
    c(" " = 4, "$\\\\DD=0$" = 2, "$\\\\DD=1$" = 2),
    escape = FALSE
  ) %>%
  cat()
sink()

