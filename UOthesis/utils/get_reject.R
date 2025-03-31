library(tidyverse)

get_reject <- function(results){
  if(is.null(results)) return(NULL)

  pvals <-
    mapply(
      FUN = function(mat){
        mat[,2]
      },
      mat = results
    ) %>%
    t() %>%
    apply(
      MARGIN = 2,
      FUN = function(x) (x < 0.05) %>% mean()
    )

  res <- c("C" = NA, "M" = NA, "E:C" = NA, "E:M" = NA, "E:Im" = NA, "E:If" = NA, "Im" = NA, "If" = NA)
  elements <- intersect(c("C", "M", "E:C", "E:M", "E:Im", "E:If", "Im", "If"), names(pvals))
  res[elements] <- pvals[elements]

  return(res)
}

## Works for haplin too
get_bias_emim <- function(results, C, M0, M1, V, Im, If){
  if(is.null(results)) return(NULL)

  point_est_means <-
    mapply(
      FUN = function(mat){
        eff <- mat$effects
        return(eff)
      },
      mat = results
    )

  ## Sometimes EMIM produces values exactly equal to 1, this ought to cover that almost always...
  ## Also making it work for haplin
  elements <- intersect(names(point_est_means[point_est_means[,1] != 1, 1]), c("C", "M[E=1]", "M[E=0]", "M[E=1]/M[E=0]", "Im", "If", "RRcm", "RRcf"))

  point_est_means <- point_est_means %>%
    t() %>%
    colMeans() %>%
    .[elements]

  true_vals <- c(C = C, "M[E=0]" = M0, "M[E=1]" = M1, "M[E=1]/M[E=0]" = V, Im = Im, If = If, "RRcm" = C, "RRcf" = C)[elements]
  bias <- (true_vals - point_est_means)/true_vals
  return(bias)
}

get_trill_reject_table <- function(res, conditions){
  mt <- c("HWE", "MS", "MaS")

  lapply(
    X = 1:3,
    FUN = function(b){
      lapply(
        X = 1:nrow(conditions),
        FUN = function(a, b){
          with(
            conditions[a,],
            rbind(get_reject(res[[b]][[a]]$strat), get_reject(res[[b]][[a]]$nostrat)) %>%
              as.data.frame() %>%
              mutate(strat = if (nrow(.) < 2) "nostrat" else c("strat", "nostrat"),
                     includeControl = includeControl,
                     rowid = rowid
              )
          )
        },
        b = b
      ) %>%
        bind_rows() %>%
        mutate(mt = mt[[b]])
    }
  ) %>%
    bind_rows()
}

get_true_effects_nonvec <- function(conditions, id) {
  effects_tex <- c("$\\CC$", "$\\MM$", "$\\EE:I_{\\MM}$", "$\\EE:I_{\\FF}$", "$I_{\\MM}$", "$I_{\\FF}$")
  conditions %>%
    filter(rowid == id, !includeControl) %>%
    with(
      {
        effects <- c(R, S, ifelse(Einteraction == "Im" & includeE, V, 1), ifelse(Einteraction == "If" & includeE, V, 1), Im, If)
        tex_output <- paste(effects_tex[effects != 1], sep = ", ")
        if (identical(tex_output, character(0)))
          tex_output <- "None"
        return(tex_output)
      }
    )
}

get_true_effects <- Vectorize(get_true_effects_nonvec, vectorize.args = "id")
