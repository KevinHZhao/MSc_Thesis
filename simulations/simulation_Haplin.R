library(TriLLIEM)
library(tidyverse)
library(parallel)

## MII = 1: 11.3 s data, 1.40 min trill, 3.69 min hap
## MII = 2: 10.78 s data, 1.48 min trill, 3.68 min hap
## MIf = 1, Im = 2: 10.55 s data, 1.94 min, 4.15 min hap

nsim <- 2000
cores <- 10

my.cluster <- makeCluster(cores)
clusterEvalQ(my.cluster, {
  library(TriLLIEM)
  library(tidyverse)
})
# registerDoSNOW(my.cluster)
# pb <- txtProgressBar(max = nsim, style = 3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts <- list(progress = progress)

startT <- Sys.time()

dat <- parLapply(
  my.cluster,
  rep(1000, nsim),
  function(x){
    simulateData(
      ntrios = x, 
      S = c(1, 1, 1),
      Im = 2,
      If = 1)
  }
)

endT <- Sys.time()
print("sim:")
print(endT - startT)

TriLL_res <- parLapply(
  my.cluster,
  dat,
  function(x){
    TriLLIEM(dat = x$dat4R,
             effects = c("M", "Im", "If"),
             max.iter = 30) %>%
      summ_trill(effects = c("M", "Im", "If")) %>%
      unlist()
  }
) %>%
  do.call(rbind, .)

endT <- Sys.time()
print("trill:")
print(endT - startT)

haplin_res <- parLapply(
  my.cluster,
  dat,
  function(x){
    runHaplin(dat = x$dat4haplin, 
              PoO = TRUE) %>% 
      summ_haplin(PoO = TRUE) %>%
      unlist()
  }
) %>%
  do.call(rbind, .)

endT <- Sys.time()
print("haplin:")
print(endT - startT)
stopCluster(cl = my.cluster)

# sim_results <-
#   foreach(
#     k = 1:nsim,
#     .packages = c("TriLLIEM", "dplyr"),
#     .combine = rbind,
#     .options.snow = opts
#   ) %dopar% {
#     dat <-
#       simulateData(
#         ntrios = 1000,
#         R = c(1, 1, 1),
#         S = c(1, 1, 1)
#       )
#     TriLL_res <-
#       TriLLIEM(dat = dat$dat4R,
#                includeI = TRUE,
#                max.iter = 30)
#     emim_res <- runEMIM(peddat = dat$dat4EMIM, includeI = TRUE) %>% summ_emim()
#     res <- data.frame(
#       T.C = exp(coef(TriLL_res))["C"],
#       T.M = exp(coef(TriLL_res))["M"],
#       T.I = exp(coef(TriLL_res))["Im"],
#       E.C = emim_res$effects["C"],
#       E.M = emim_res$effects["M"],
#       E.I = emim_res$effects["Im"],
#       p.T.C = coef(summary(TriLL_res))["C", 4],
#       p.T.M = coef(summary(TriLL_res))["M", 4],
#       p.T.I = coef(summary(TriLL_res))["Im", 4],
#       p.E.C = emim_res$pvals["C"],
#       p.E.M = emim_res$pvals["M"],
#       p.E.I = emim_res$pvals["Im"]
#     )
#     return(res)
#   }
# close(pb)