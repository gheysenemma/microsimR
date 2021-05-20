#' @title Lowest-RMSE sample from time series
#' @description Computes the Root Mean Square Error between subsequent generations (columns)
#' in given species abundances matrix (time series) as a measure of convergence.
#' The timepoint with the lowest RMSE-difference with respect to the previous timepoint
#' is then extracted and a normalised sample is returned (proportions).
#' @param spab Species abundances matrix with OTUs in the rows and the timepoints as columns.
#' @param warn default TRUE, set to FALSE to suppress convergence warning
#' @param norm default TRUE: compute RMSE on compositional time series (abundances per time point sum to 1)
#' @param cutoff The value the minimum RMSE cannot exceed in order to declare convergence.
#' @return A sample vector with the length equal to the number of rows of given input species abundances matrix.
#' @examples
#' spab <- glv(N = 10, A = powerlawA(n = 10, alpha = 1.2), tend = 10000)
#' sample <- rmse_sample(spab)
#' @export


rmse_sample <- function(spab, warn = TRUE, norm = TRUE, cutoff = 1e-04){
  if(norm){
    rmse_vec <- rmse_t(spab = t(t(spab)/colSums(spab)))
  } else {
    rmse_vec <- rmse_t(spab)
  }
  if(warn){
    if(min(rmse_vec) >= (cutoff)){
      warning(paste0("WARNING -  simulation did not converge. ",
                     "Run the simulation longer or change parameters."))
    }
  }
  sample <- spab[,which.min(rmse_vec)]
  return(sample)
}

# Helper function to compute rmse vector over all timepoints of a given species abundances x time points matrix
rmse_t = function(spab){
  N = nrow(spab)
  K = ncol(spab)

  rmse_vec <- rep.int(0, times = (K-1)) # initialize
  rmse_vec <- sapply(X = 1:(K-1), FUN = function(j){
    rmse_vec[j] = rmse(spab[,j], spab[,(j+1)])
  })

  return(rmse_vec)
}
