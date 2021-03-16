\name{soi}
\alias{soi}
\title{
Self-Organised Instability Model Community Simulation
}
\description{
Generate time-series with the self-organised instability (SOI) model.
Implements a K-leap method for accelerating stochastic simulation.
}
\usage{
soi(
  N,
  I,
  A = generateA(N, c = 0.1),
  migr_rates = runif(N, min = 0.1, max = 0.8),
  death_rates = runif(N, min = 0.01, max = 0.08),
  com = NULL,
  tend,
  k = 5,
  perturb = NULL,
  norm = FALSE
)
}
\arguments{
  \item{N}{
number of species
}
  \item{I}{
community size, number of available sites (individuals)
}
  \item{A}{
interaction matrix of dimension NxN
}
  \item{migr_rates}{
species-specific immigration probabilities (also determine initial abundances)
}
  \item{death_rates}{
species-specific extinction probabilities
}
  \item{com}{
Initial community abundances vector. If NULL (default), based on migration rates.
}
  \item{tend}{
number of timepoints to be returned in the time series (nr of generations)
}
  \item{k}{
the number of transition events that are allowed to take place during one leap.
By default set to 5. Higher values reduce runtime, but also accuracy of the simulation.
}
  \item{perturb}{
perturbation object
}
  \item{norm}{
logical to indicate whether the time series should be returned
with the abundances as proportions (norm = TRUE) or the raw counts (norm = FALSE, default)
}

}

\references{
Solé RV, Alonso D, McKane A. Self-organized instability in complex ecosystems. Philos Trans R Soc Lond B Biol Sci. 2002;357(1421):667-671. doi:10.1098/rstb.2001.0992

Cai X, Xu Z. K-leap method for accelerating stochastic simulation of coupled chemical reactions. J Chem Phys. 2007;126(7):074102. doi:10.1063/1.2436869
}

\value{
matrix with species abundances as rows and time points as columns
}
\author{
Gheysen E.
}
\examples{
soi(N = 10, I = 1000, A = generateA(10, c = 0.2), tend = 1000)
}
