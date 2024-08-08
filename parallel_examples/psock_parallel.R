library(parallel) # Lets us set up parallel clusters for either method, comes installed with R

## With PSOCK, each core creates a new R environment.
## Loading packages in THIS session WON'T let the cores use their functions.
## Our cores also can't see any of the variables in this environment.
## With the same example as before...

library(dplyr) # Load dplyr package in THIS R session

n.cores <- 10 # Use 10 cores
cl <- makeCluster(n.cores, type = "PSOCK")

## We can load libraries using clusterEvalQ
clusterEvalQ(cl, {
  library(dplyr)
})
## We can export variables using clusterExport
message <- "test"
clusterExport(cl, "message")
clusterEvalQ(cl, print(message))

data <- parLapply( # Simulates a list of 1000, 1000-row data frames that can each be fit to a lm
  cl = cl,
  X = 1:1000, # plugging in x = 1, 2, 3, ... to FUN, won't do anything
  fun = function(x){ # function that produces a 1000-row data frame simulating the lm y = 2 + 3x1 + 5x2
    data.frame(
      x1 = runif(1000),
      x2 = runif(1000),
      eps = rnorm(1000)
    ) %>% # We use the magrittr pipe here, as well as a dplyr function!
      mutate(
        y = 2 + 3*x1 + 5*x2 + eps
      )
  }
) # FUN uses functions from dplyr, errors if dplyr isn't exported to our cluster!

# Now we can fit the lm 1000 times, with parallelization!
system.time({
  fits <- parSapply(
    cl = cl,
    X = data, ## Each data frame in the list data will be given to FUN
    FUN = function(x){ # x = our 1000-row data frame
      lm(y ~ x1 + x2, data = x) %>%
        coef() # Extracts the fitted coefficients from the lm
    }
  )
})
stopCluster(cl = cl)

# We can compare the speed to if we don't parallelize
system.time({
  fits_nopar <- sapply(
    X = data, ## Each data frame in the list data will be given to FUN
    FUN = function(x){ # x = our 1000-row data frame
      lm(y ~ x1 + x2, data = x) %>%
        coef() # Extracts the fitted coefficients from the lm
    }
  )
})

identical(fits, fits_nopar) # Same value, but time saved!
## Moral of the story, FORK > PSOCK, unless you want to use Windows
