## Calculates the Fisher info matrix, only for the imprinting case

#' Title
#'
#' @param res
#'
#' @return
#' @export
#' @keywords internal
#'
#' @examples
vcov.TriLLIEM <- function(res) {
  ## Calculating SE's
  ## IX is the 16-row info matrix
  ## IY is the 15-row info matrix

  nsets <- nobs(res)/16

  if(nsets %% 1 != 0) {
    return(NextMethod())
  }

  ## I_X <- vcov(res) is slightly less precise
  Z <- res$x
  I_X <- t(Z) %*% diag(res$fitted.values) %*% Z

  ## Matrix that transforms complete data (X) into observed data (Y)
  L <- diag(16 * nsets)
  L[9L + 16L * (1:nsets - 1),] <- L[9L + 16L * (1:nsets - 1),] + L[10L + 16L * (1:nsets - 1),]
  L[c(outer(10:15,16L * (1:nsets - 1), FUN = "+")),] <- L[c(outer(11:16,16L * (1:nsets - 1), FUN = "+")),]
  L <- L[-(16L * (1:nsets)),]

  ## Try to clean up so don't need to calculate all the zeros
  var_XgY <- vapply(X = 1:(nsets * 16),
                    FUN = function(x){
                      vapply(X = 1:(nsets * 16),
                             cov_trill,
                             a = x,
                             y = L %*% res$y,
                             L = L,
                             mu = fitted(res),
                             FUN.VALUE = double(1)
                      )
                    },
                    FUN.VALUE = double(nsets * 16)
  )
  I_XgY <- t(Z) %*% var_XgY %*% Z
  I_Y <- I_X - I_XgY
  return(solve(I_Y))
}
