#' @title Simulate time series with the generalized Lotka-Volterra model
#'
#' @description Simulate a community time series using the generalized Lotka-Volterra model, defined as
#' \eqn{\frac{dx}{dt}=x(growth_rates+Ax)}, where x is the vector of species abundances, A is the interaction matrix
#' and growth_rates the vector of growth rates.
#'
#' @param N species number
#' @param A interaction matrix
#' @param growth_rates growth rates
#' @param com initial abundances
#' @param tend final time point
#' @return a matrix with species abundances as rows and time points as columns, column names give time points
#' @export

glv<-function(N=4,A,growth_rates=runif(N),com=runif(N),tend=100){
  # checks perturb object and includes in parms
  # parms as matrix
  parms = cbind(rep(N,N), growth_rates, A)
  times<-seq(0, tend, by=1)
  # run the simulation
  commtime<-lsoda(
    y = com, # initial state values for the ODE system
    times = times, # times at which the explicit estimates for com are desired
    func = glvsolve, # # function that computes the values of the derivatives in the ODE system
    parms = parms, # parameters used in func
    hmax = 0, # 0: no maximal integration stepsize specified
    maxsteps = 500000, # max nr of steps per output interval taken by the solver
  )
  #time=commtime[,1]
  commtime=t(commtime[,2:ncol(commtime)])
  colnames(commtime)=paste0("t", 1:ncol(commtime))
  rownames(commtime)=paste0("sp", 1:N)
  return(commtime)
}

# ==============================================================
# Equations (generalized Lotka-Volterra)
# ==============================================================

# matrix formulation of the ODE set
# t: current simulation time
# com: vector with current values of state variables (initial conditions)
# parms: parameter values
#
glvsolve<-function(t, com, parms){
  N=parms[1,1]  # species number
  growth_rates=parms[,2]   # vector of growth rates
  A=parms[,3:(N+2)] # interaction matrix
  dydt <- com*(growth_rates+A %*% com)
  list(dydt)
}
