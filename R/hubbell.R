#' @title Hubbell community simulation
#' @description Neutral species abundances simulation according to the Hubbell model.
#' @param N Amount of different species initially in the local community
#' @param M Amount of different species in the metacommunity, including those of the local community
#' @param I Fixed amount of individuals in the local community
#' @param d Fixed amount of deaths of local community individuals in each generation
#' @param pbirth probabilities of birth to construct the initial abundances,
#' will be updated in all next generations with previous abundances
#' @param m Immigration rate: the probability that a death in the local community is
#' replaced by a migrant of the metacommunity rather than by the birth of a local community member
#' @param pmigr probabilities of migration
#' @param tskip Nr of generations that should not be included in the outputted species
#' abundance matrix.
#' @param tend Nr of simulations to be simulated.
#' @param norm logical to indicate whether the time series should be returned
#' with the abundances as proportions (norm = TRUE) or the raw counts
#' (norm = FALSE, default)
#' @examples hubbell(N = 8, M = 10, I = 1000, d = 50, m = 0.02, tend = 100)
#' @return matrix with species abundances as rows and time points as columns
#' @references Rosindell, James et al. “The unified neutral theory of biodiversity and biogeography at age ten.” Trends in ecology & evolution vol. 26,7 (2011).
#' @export

hubbell <- function(
  N, # amount of local species
  M, # amount of meta species, incl local species (total species)
  I = 1000, # community size (nr of individuals)
  d = 10, # nr of deaths per generation
  m = 0.02, # immigration rate (probability dead indiv replaced by meta-indiv)
  pbirth = runif(N, min = 0, max = 1), # probabilities of birth
  pmigr = runif(M, min = 0, max = 1), # probabilities of migration
  tskip = 0, # amount of timepoints to not return
  tend, # amount of time points
  norm = FALSE
){
  #####################################################################################
  # First setting the function arguments right
  if(length(pbirth)!=N | length(pmigr)!=M){
    stop("Either length of pbirth vector does not match with N or length of pmigr vector does not match with M")
  }
  pbirth <- c(pbirth, rep(0, times = (M-N)))
  pbirth <- pbirth/sum(pbirth)
  pmigr <- pmigr/sum(pmigr)
  com <- ceiling(I*pbirth)
  if(sum(com)<I){
    ind <- sample(1:M, size = I-sum(com), prob = pbirth)
    com[ind] <- com[ind] +1
  } else if(sum(com)>I){
    ind <- sample(1:M, size = sum(com)-I, prob = 1-pbirth)
    com[ind] <- com[ind] -1
  }

  #################################################################################
  # The simulation

  tseries <- matrix(0, nrow = M, ncol = tend)
  colnames(tseries) <- paste0("t", 1:tend)
  rownames(tseries) <- 1:M

  com[which(com < 0)] <- 0
  tseries[,1] <- com
  for (t in 2:tend){
    # Each iteration the probability of births is updated by the counts
    pbirth <- com/sum(com)
    pbirth[which(pbirth < 0)] <- 0
    # Probability of births is used to pick the species that will die
    # because species with count 0 will have probability 0 and species not present in the community can also not die
    deaths <- rmultinom(n = 1, size = d, prob = pbirth)
    while(sum(com-deaths <0) >0){
      neg_sp <- which(com-deaths <0)
      pbirth[neg_sp] <- 0
      deaths <- rmultinom(n = 1, size = d, prob = pbirth)
    }

    event <- rbinom(d, 1, prob = m) # immigration rate m: probability death replaced by immigrant
    # immigration 1, birth 0
    n_migrants <- sum(event)
    n_births <- length(event) - n_migrants

    births <- rmultinom(1, n_births, prob = pbirth)
    migr <- rmultinom(1, n_migrants, prob = pmigr)
    com <- com - deaths + births + migr
    com[which(com < 0)] <- 0
    tseries[,t] <- com
  }
  if(norm){
    tseries <- t(t(tseries)/colSums(tseries))
  }

  return(tseries[, (tskip +1):tend])
}
