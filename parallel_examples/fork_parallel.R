library(parallel) # Lets us set up parallel clusters for either method, comes installed with R

## With FORK, each core works on the same, original R environment.
## If we load a package, all the cores can access it, if we create a
## variable/function, all cores can use it.
## Example: Trying to fit 1000 random data sets to a linear model.

library(dplyr) # Load dplyr package in THIS R session

n.cores <- 10 # Use 10 cores
data <- mclapply( # Simulates a list of 1000, 1000-row data frames that can each be fit to a lm
  X = 1:1000, # plugging in x = 1, 2, 3, ... to FUN, won't do anything
  FUN = function(x){ # function that produces a 1000-row data frame simulating the lm y = 2 + 3x1 + 5x2
    data.frame(
      x1 = runif(1000),
      x2 = runif(1000),
      eps = rnorm(1000)
    ) %>% # We use the magrittr pipe here, as well as a dplyr function!
      mutate(
        y = 2 + 3*x1 + 5*x2 + eps
      )
  },
  mc.cores = n.cores
) # FUN uses functions from dplyr, but it works since we FORK!

# Now we can fit the lm 1000 times, with parallelization!
system.time({
  fits <- mcmapply(
    FUN = function(x){ # x = our 1000-row data frame
      lm(y ~ x1 + x2, data = x) %>%
        coef() # Extracts the fitted coefficients from the lm
    },
    x = data, ## Each data frame in the list data will be given to FUN
    mc.cores = n.cores
  )
})

# We can compare the speed to if we don't parallelize
system.time({
  fits_nopar <- mapply(
    FUN = function(x){ # x = our 1000-row data frame
      lm(y ~ x1 + x2, data = x) %>%
        coef() # Extracts the fitted coefficients from the lm
    },
    x = data ## Each data frame in the list data will be given to FUN
  )
})

identical(fits, fits_nopar) # Same value, but time saved!
