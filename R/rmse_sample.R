#' @title Lowest-RMSE sample from time series
#' @description Computes the Root Mean Square Error between subsequent generations (columns)
#' in given species abundances matrix (time series) as a measure of convergence.
#' The timepoint with the lowest RMSE-difference with respect to the previous timepoint
#' is then extracted and a normalised sample is returned (proportions).
#' @param spab Species abundances matrix with OTUs in the rows and the timepoints as columns.
#' @return A sample vector with the length equal to the number of rows of given input species abundances matrix.
#' @export


rmse_sample <- function(spab){
  spab <- spab/colSums(spab)
  rmse <- rmse_t(spab)
  if(min(rmse) > (1e-04)){
    warning(paste0("WARNING -  simulation did not converge. ",
                   "Run the simulation longer or change parameters."))
  }
  sample <- spab[,which.min(rmse)]
  return(sample)
}

# Helper function to compute rmse vector over all timepoints of a given species abundances x time points matrix
rmse_t = function(spab){
  N = nrow(spab)
  K = ncol(spab)

  rmse <- rep.int(0, times = (K-1)) # initialize
  rmse <- sapply(X = 1:(K-1), FUN = function(j){
    rmse[j] = Metrics::rmse(spab[,j], spab[,(j+1)])
  })

  return(rmse)
}
