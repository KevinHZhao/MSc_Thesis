library(tidyverse)
library(knitr); options(knitr.kable.NA = '')
library(kableExtra)

df <-
  readRDS("../data/consistencydat.RDS")$dat4R %>%
  select(mt_MS, M, F, C, E, D, count)

n <- nrow(df)
half <- ceiling(n / 2)
df1 <- df[1:half, ]
df2 <- df[(half+1):n, ]

tab1 <- kbl(df1,
            format = "latex",
            booktabs = TRUE,
            row.names = FALSE,
            align = "l",
            col.names = c("Mating type", "$\\MM$", "$\\FF$", "$\\CC$", "$\\EE$", "$\\DD$", "count"),
            escape = FALSE,
            linesep = c(rep("", 14), "\\addlinespace"))

tab2 <- kbl(df2,
            format = "latex",
            booktabs = TRUE,
            row.names = FALSE,
            align = "l",
            col.names = c("Mating type", "$\\MM$", "$\\FF$", "$\\CC$", "$\\EE$", "$\\DD$", "count"),
            escape = FALSE,
            linesep = c(rep("", 14), "\\addlinespace"))

tex_output <- paste0(
  "\\begin{table}[htb]\n",
  "\\centering\n",
  "\\caption{Simulated null data set for comparison between EMIM, Haplin, and \\texttt{TriLLIEM}.}\n",
  "\\label{tab:consistencydat}\n",
  "\\begin{minipage}{0.48\\textwidth}\n",
  tab1,
  "\\end{minipage}\n",
  "\\hfill\n",
  "\\begin{minipage}{0.48\\textwidth}\n",
  tab2,
  "\\end{minipage}\n",
  "\\end{table}"
)

writeLines(tex_output, "consistencydat.tex")

# sink(paste0("consistencydat.tex"))
#  %>%
#   kbl(
#     format = "latex",
#     align = "l",
#     col.names = c("Mating type", "$\\MM$", "$\\FF$", "$\\CC$", "$\\EE$", "$\\DD$", "count"),
#     caption = "Simulated null data set for comparison between EMIM, Haplin, and \\texttt{TriLLIEM}.",
#     label = "consistencydat",
#     escape = FALSE,
#     booktabs = TRUE,
#     linesep = c(rep("", 14), "\\addlinespace"),
#     position = "H"
#   ) %>%
#   cat()
# sink()
