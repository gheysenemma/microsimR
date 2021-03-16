#' @title Hubbell community simulation
#' @description Neutral species abundances simulation according to the Hubbell model.
#' @details For more theoretical background see: Rosindell, James et al. “The unified neutral theory of biodiversity and biogeography at age ten.” Trends in ecology & evolution vol. 26,7 (2011).
#' @param I Fixed amount of individuals in the local community
#' @param N Amount of different species initially in the local community
#' @param M Amount of different species in the metacommunity, including those of the local community
#' @param d Fixed amount of deaths of local community individuals in each generation
#' @param pbirth probabilities of birth
#' @param m Immigration rate: the probability that a death in the local community is
#' replaced by a migrant of the metacommunity rather than by the birth of a local community member
#' @param pmigr probabilities of migration
#' @param com initial community vector of length M where the ith element is the amount of individuals of the ith species
#' @param tskip Nr of generations that should not be included in the outputted species
#' abundance matrix.
#' @param tend Nr of simulations to be simulated.
#' @param norm logical to indicate whether the time series should be returned
#' with the abundances as proportions (norm = TRUE) or the raw counts
#' (norm = FALSE, default)
#' @examples hubbell(N = 8, M = 10, I = 1000, d = 50, m = 0.02, tend = 100)
#' @return matrix with species abundances as rows and time points as columns
#' @export

hubbell <- function(
  N = NULL, # amount of local species
  M = NULL, # amount of meta species, incl local species (total species)
  I, # community size (nr of individuals)
  d, # nr of deaths per generation
  m, # immigration rate (probability dead indiv replaced by meta-indiv)
  com = NULL, # initial community vector of length M where the ith element is the amount of individuals of the ith species
  pbirth = NULL, # probabilities of birth
  pmigr = NULL, # probabilities of migration
  tskip = 0, # amount of timepoints to not return
  tend, # amount of time points
  norm = FALSE
){
  #####################################################################################
  # First setting the function arguments right

  if(is.null(N) & is.null(pbirth)){
    stop("Please give the amount of local species N or a birth probability vector of length N")
  } else if (is.null(N)){
    if(sum(pbirth)==1){
      N = length(pbirth[pbirth>0])
    } else {
      stop("pbirth probability vector must be values between 0 and 1, and sum up to 1")
    }
  } else if (is.null(pbirth)){
    pbirth = c(
      # Local species: probabilities of birth from uniform distribution
      runif(N, min = 0, max = 1),
      # Meta species: probabilities of birth initially set to 0
      rep(0, times = (M - N))
    )
    pbirth = pbirth/sum(pbirth)
  } else if (sum(pbirth>0)>N){
    stop("do not give more nonzero probabilties in the pbirth vector than there are local species")
  } else if (!is.null(pbirth)){
    pbirth = c(pbirth, rep(0, times = (M-N)))
    pbirth = pbirth/sum(pbirth)
  }
  sp_names <- paste0("sp_", 1:M)
  names(pbirth) <- sp_names

  if(is.null(M) & is.null(pmigr)){
    stop("Please give the amount of meta species M (incl local) or a migration probability vector of length M")
  } else if (is.null(M)){
    if(sum(pmigr)==1){
      M = length(pmigr)
    } else {
      stop("pmigr probability vector must be values between 0 and 1, and sum up to 1")
    }
  } else if (is.null(pmigr)){
    pmigr = runif(M, min = 0, max = 1)
    pmigr = pmigr/sum(pmigr)
  } else if (length(pmigr)!=M){
    stop("pmigr vector must be of length M")
  }
  names(pmigr) <- sp_names

  if(is.null(com)){
    com = round(I*pbirth)
  } else if (length(com)!=M){
    stop("length of com vector must be equal to M")
  }
  names(com) <- sp_names

  #################################################################################
  # The simulation

  tseries <- matrix(0, nrow = M, ncol = tend)
  colnames(tseries) <- paste0("t", 1:tend)
  rownames(tseries) <- sp_names

  tseries[,1] <- com
  for (t in 2:tend){
    # Each iteration the probability of births is updated by the counts
    pbirth <- com/sum(com)
    pbirth[pbirth < 0] <- 0
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
    tseries[,t] <- com
  }
  if(norm){
    tseries <- tseries/colSums(tseries)
  }

  return(tseries[, (tskip +1):tend])
}
