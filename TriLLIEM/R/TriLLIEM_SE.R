## Calculates the Fisher info matrix, only for the imprinting case

TriLLIEM_SE <- function(res) {
  ## Calculating SE's
  ## IX is the 16-row info matrix
  ## IY is the 15-row info matrix

  nsets <- nobs(res)/16

  if(nsets %% 1 != 0) {
    stop("Number of observations is not a factor of 16.")
  }

  ## I_X <- vcov(res) is slightly less precise
  Z <- res$x
  I_X <- t(Z) %*% diag(res$fitted.values) %*% Z

  ## Matrix that transforms complete data (X) into observed data (Y)
  L <- diag(16 * nsets)
  L[9L + 16L * (1:nsets - 1),] <- L[9L + 16L * (1:nsets - 1),] + L[10L + 16L * (1:nsets - 1),]
  L[c(outer(10:15,16L * (1:nsets - 1), FUN = "+")),] <- L[c(outer(11:16,16L * (1:nsets - 1), FUN = "+")),]
  L <- L[-(16L * (1:nsets)),]

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
  se <- sqrt(diag(solve(I_Y)))
  return(se)
}
