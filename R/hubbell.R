#' @title Hubbell community simulation
#' @description Neutral species abundances simulation according to the Hubbell model.
#' @details For more theoretical background see: Rosindell, James et al. “The unified neutral theory of biodiversity and biogeography at age ten.” Trends in ecology & evolution vol. 26,7 (2011).
#' @param I Fixed amount of individuals in the local community
#' @param N Amount of different species initially in the local community
#' @param M Amount of different species in the metacommunity, including those of the local community
#' @param d Fixed amount of deaths of local community individuals in each generation
#' @param mprob Species proportions in the metacommunity, or migration probability vector
#' for each species of the metacommunity to replace a death in the local community.
#' When given, length of the vector must be equal to M. By default set to 'NA':
#' normalised random deviates from the uniform distribution are used.
#' @param m Immigration rate: the probability that a death in the local community is
#' replaced by a migrant of the metacommunity rather than by the birth of a local community member
#' @param lprob Species proportions in the local community, or birth probability vector
#' for each species of the local community to replace a death in the local community.
#' When given, length of the vector must be equal to N.
#' By default set to 'NA': normalised random deviates from the uniform distribution
#' are used.
#' @param lpop Vector of the species indexes initially in the local community.
#' When given, length of the vector must be equal to I. By default set to 'NA':
#' I species indexes are initially randomly sampled from 1 to N.
#' @param tskip Nr of generations that should not be included in the outputted species
#' abundance matrix.
#' @param tend Nr of simulations to be simulated.
#' @param norm logical to indicate whether the time series should be returned
#' with the abundances as proportions (norm = TRUE) or the raw counts
#' (norm = FALSE, default)
#' @return matrix with species abundances as rows and time points as columns
#' @export

hubbell <- function(
  I, # Total individuals to be considered in the local community
  N, # Amount of species initially in the local community
  M, # Amount of metacomunnity species (INCLUDING LOCAL COMMUNITY SPECIES?)
  d = 10, # fixed amount of deaths in each generation
  mprob = NA, #Species proportions in metacommunity - probability for a species to migrate to the local community
  m = 0.02, # immigration rate (probability death individuum is replaced by individuum of local community)
  lprob = NA, #Species proportions in the local community - probabilities for the local species for 'birth' in the local community
  lpop = NA, # Contains the species index that are in the local community
  tskip = 0, # Timepoints to not include in the returned time series
  tend = 1000, # Timesteps
  norm = FALSE
){
  options(warn = -1) # supress warning messages (annoying in loops)

  # Getting the parameters correct:

  # If mprob / lprob are not given as an argument, use normalised random deviates from uniform distribution for probabilities
  if (sum(is.na(mprob)) > 0){ # equivalent to if set to NA (type conversion prob fixed)
    mprob =  runif(M, min = 0, max = 1)
    mprob = mprob/sum(mprob) # normalising: need to sum up to 1
  } else if (length(mprob) != M){
    stop("mprob vector must contain M elements, for each species of the metacommunity a migration probability")
  }

  if (sum(is.na(lprob)) > 0){
    lprob = runif(N, min = 0, max = 1)
    lprob = lprob/sum(lprob)
  } else if (length(lprob) != N){
    stop("lprob vector must contain N elements, for each species of the local community a birth probability")
  }

  if (sum(is.na(lpop)) > 0){
    lpop = sample(1:N, I, replace = TRUE)
  } else if (length(lpop) != I){
    stop("lpop vector must contain I elements, a species index for each individual initially in the local community")
  }

  # Initialize output time series matrix
  tseries=matrix(NA,nrow=M,ncol=(tend-tskip))

  # Add count of initial local community to tseries matrix
  ## only if timestep not to be skipped
  if(tskip == 0){
    tseries[,1] = countSpecies(lpop, M)
  }

  # Now loop over remaining of time series
  for (t in 2:tend){
    # lpop death coordinates of the species to be killed
    dcoor = sample(1:I, d)

    # Determine replacement of each death:
    ## 1 = birth (replaced by birth of local community species)
    ## 0 = migration (replaced by migration of a metacommunity species to the local community)
    event = rbinom(n = d, size = 1, prob = m)
    event[which(is.na(event))] = 0 # omitting warning message

    # Which species will replace the deaths?
    births = ((1:N)%*%rmultinom(n = d, size=1, prob=lprob))*event
    migrants = ((1:M)%*%rmultinom(n = d, size = 1, prob = mprob))*as.numeric(!event)
    replacers = births + migrants

    lpop[dcoor] = replacers

    if(t > tskip){
      tseries[,(t-tskip)] = countSpecies(lpop, M)
    }
  }
  colnames(tseries) = seq(from = tskip+1, to = tend)
  rownames(tseries) = 1:nrow(tseries)

  options(warn = 0) # turn warning messages back on
  if(norm){
    tseries <- tseries/colSums(tseries)
  }

  return(tseries)
}

# countSpecies function also reviewed: much shorter possible
countSpecies <- function(pop, N){
  counts = rep(0, N)
  for(s in pop){
    counts[s] = length(which(pop==s))
  }
  return(counts)
}
