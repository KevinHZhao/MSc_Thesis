library(tidyverse)

get_bias <- function(results, C, M, V, Im, If, Einteraction){
  if(is.null(results)) return(NULL)# return(c(C = NA, M = NA, "E:M" = NA, "E:Im" = NA, "E:If" = NA, Im = NA, If = NA))

  point_est_means <-
    mapply(
      FUN = function(mat){
        mat[,1]
      },
      mat = results
    )

  elements <- intersect(c("C", "M", "E:C", "E:M", "E:Im", "E:If", "Im", "If"), row.names(point_est_means))

  point_est_means <- point_est_means %>%
    t() %>%
    colMeans() %>%
    exp() %>%
    .[elements]

  true_vals <- c(C = C, M = M, "E:C" = 1, "E:M" = if (Einteraction == "M") V else 1, "E:Im" = if (Einteraction == "Im") V else 1, "E:If" = if (Einteraction == "If") V else 1, Im = Im, If = If)
  true_vals[!(names(true_vals) %in% elements)] <- NA
  filtered_tv <- true_vals[elements]
  if (!identical(names(point_est_means), names(filtered_tv)))
    stop(paste("Somehow the names of the values are out of order!\nTrue:", paste(names(point_est_means)), "\nFiltered:", paste(names(filtered_tv))))
  bias <- true_vals
  bias[elements] <- (point_est_means - filtered_tv)/filtered_tv
  return(bias)
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

get_trill_table <- function(res, conditions){
  mt <- c("HWE", "MS", "MaS")

  lapply(
    X = 1:3,
    FUN = function(j) {
      df <- lapply(
        X = 1:nrow(conditions),
        FUN = function(i, j) {
          with(conditions[i,],
               rbind(
                 get_bias(
                   results = res[[j]][[i]]$strat,
                   C = R,
                   M = S,
                   V = V,
                   Im = Im,
                   If = If,
                   Einteraction = Einteraction
                 ),
                 get_bias(
                   results = res[[j]][[i]]$nostrat,
                   C = R,
                   M = S,
                   V = V,
                   Im = Im,
                   If = If,
                   Einteraction = Einteraction
                 )
               ) %>%
                 as.data.frame(check.names = FALSE) %>%
                 mutate(strat = if (nrow(.) < 2) "nostrat" else c("strat", "nostrat"),
                        includeControl = includeControl,
                        rowid = rowid)
               )
        },
        j = j
      ) %>%
        lapply(as.data.frame.list, check.names = FALSE) %>%
        bind_rows() %>%
        mutate(effect = colnames(df),
               mt = mt[[j]]) %>%
        return()
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
