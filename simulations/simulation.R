library(TriLLIEM)
library(tidyverse)
library(parallel)

## CMI = 1: 12.8 s data, 1.56 min trill, 2.86 min EMIM
## CMI = 2: 12.9 s data, 2.2 min trill, 3.51 min EMIM
## CM = 1, I = 2: 13 s data, 2.18 min trill, 3.5 min EMIM

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
      S = c(1, 1, 1))
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
             includeIm = TRUE,
             max.iter = 30) %>%
      summ_trill(effects = c("C", "M", "Im")) %>%
      unlist()
  }
) %>%
  do.call(rbind, .)

endT <- Sys.time()
print("trill:")
print(endT - startT)

emim_res <- parLapply(
  my.cluster,
  dat,
  function(x){
    runEMIM(peddat = x$dat4EMIM, 
            includeI = TRUE) %>% 
      summ_emim() %>%
      unlist()
  }
) %>%
  do.call(rbind, .)

endT <- Sys.time()
print("emim:")
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