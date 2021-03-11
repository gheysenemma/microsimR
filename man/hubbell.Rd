\name{hubbell}
\alias{hubbell}
\title{
Hubbell community simulation
}
\description{
Neutral species abundances simulation according to the Hubbell model.
}
\usage{
hubbell(N,
  M,
  I,
  d,
  pmigr = NULL,
  m,
  pbirth = NULL,
  com = NULL,
  tskip = 0,
  tend,
  norm = FALSE)
}
\arguments{
  \item{N}{
Amount of different species initially in the local community
}
  \item{M}{
Amount of different species in the metacommunity, including those of the local community
}
  \item{I}{
Fixed amount of individuals in the local community
}
  \item{d}{
Fixed amount of deaths of local community individuals in each generation
}
  \item{pbirth}{
probabilities of birth
}
  \item{m}{
Immigration rate: the probability that a death in the local community is replaced by a migrant of the metacommunity rather than by the birth of a local community member
}
  \item{pmigr}{
probabilities of migration
}
  \item{com}{
initial community vector of length M where the ith element is the amount of individuals of the ith species
}
  \item{tskip}{
Number of generations that should not be included in the outputted species abundance matrix.
}
  \item{tend}{
Number of simulations to be simulated.
}
  \item{norm}{
logical to indicate whether the time series should be returned with the abundances as proportions (norm = TRUE) or the raw counts (norm = FALSE, default)
}
}

\value{
matrix with species abundances as rows and time points as columns
}

\author{
Gheysen E.
}

\examples{
mprob = c(0.56,0.31,0.37,0.49,0.75,0.67,0.70,0.36,0.75,0.15)
lprob = c(0.55,0.12,0.03,0.58,0.12)
lpop = sample(1:5, 1000, replace = TRUE)

hubbell(I = 1000, N = 5, M = 10, d = 12,
    pmigr = mprob, pbirth = lprob, com = lpop,
    m = 0.02, tskip = 800, tend = 1000, norm = TRUE)
}

