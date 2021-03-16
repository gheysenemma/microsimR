\name{rmse_sample}
\alias{rmse_sample}
\title{Lowest-RMSE sample from time series}
\description{
Computes the Root Mean Square Error between subsequent generations (columns) in given species abundances matrix (time series) as a measure of convergence. The timepoint with the lowest RMSE-difference with respect to the previous timepoint is then extracted and a normalised sample is returned (proportions).
}
\usage{
rmse_sample(spab)
}

\arguments{
  \item{spab}{Species abundances matrix with OTUs in the rows and the timepoints as columns.}
}

\value{A sample vector with the length equal to the number of rows of given input \code{spab} species abundances matrix.}

\author{Gheysen, E.}

\details{Throws warning message if given simulation did not converge: lowest RMSE still larger than \eqn{1e-04}.}

\seealso{
\code{\link{asPhyloseq}}, \code{\link{glv}}, \code{\link{hubbell}}, \code{\link{soi}}
}
\examples{
n <- 10 # number of species
tend <- 10000 # simulations must be run long enough for convergence

A1 <- powerlawA(n = n, alpha = 1.2) # interaction matrix
spab1 <- glv(N = n, A = A1, tend = tend) # first community matrix: glv
sample1 <- rmse_sample(spab1)

A2 <- generateA(n = n, c = 0.2) # interaction matrix
spab2 <- soi(N = n, I = 1000, tend = tend) # second community matrix: soi
sample2 <- rmse_sample(spab2)

physeq <- asPhyloseq(sample1, sample2)
}
