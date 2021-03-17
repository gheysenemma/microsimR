#' @title Simulate time series with the generalized Lotka-Volterra model
#'
#' @description Simulate a community time series using the generalized Lotka-Volterra model, defined as
#' dx/dt = x(b+Ax), where x is the vector of species abundances, A is the interaction matrix
#' and growth_rates the vector of growth rates.
#'
#' @param N species number
#' @param A interaction matrix
#' @param b growth rates
#' @param x initial abundances
#' @param tend timepoints
#' @param norm return normalised abundances (proportions in each generation)
#' @return a matrix with species abundances as rows and time points as columns, column names give time points
#' @examples glv(N = 4, A = powerlawA(n = 4, alpha = 2), tend = 1000)
#' @note Calls upon the lsoda function from package deSolve to solve the ODE. If function crashes, consider running again with while loop and tryCatch handlings.
#' @seealso deSolve::lsoda
#' @export

glv <- function(
  N,
  A,
  b = runif(N),
  x = runif(N),
  tend = 1000,
  norm = FALSE
){
  # 1 Model specification
  # Model parameters
  parameters <- cbind(b, A)

  # 2 Model application
  # Time specification
  times <- seq(0, tend, by = 1)
  # Model integration
  out <- ode(
    y = x,
    times = times,
    func = dxdt,
    parms = parameters
  )
  spab <- t(out[,2:ncol(out)])
  if(norm){
    spab <- spab/colSums(spab)
  }
  return(spab)
}

# Model equations
dxdt <- function(t, x, parameters){
  b <- parameters[,1]
  A <- parameters[,2:ncol(parameters)]
  # rate of change
  dx <- x*(b+A %*% x)
  # return rate of change
  list(dx)
}
