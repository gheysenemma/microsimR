#' @title Generate random uniform interaction matrix
#' @description Generate random simplified interaction matrix from a uniform distribution.
#'
#' @param N number of species
#' @param d diagonal values (should be negative)
#' @param min.strength minimal off-diagonal interaction strength
#' @param max.strength maximal off-diagonal interaction strength
#' @param c connectance (interaction probability)
#' @examples
#' high_inter_A <- randomA(10, d = -0.4, min.strength = -0.8, max.strength = 0.8, c = 0.5)
#' low_inter_A <- randomA(10, c = 0.01)
#' @note Simplified version of generateA() function from seqtime package
#' @seealso seqtime generateA()
#' @export
#'

randomA <- function(
  N = 100, # nr of species
  d = -0.5, # diagonal values, should be negative
  min.strength = -0.5, # minimal off-diagonal interaction strength
  max.strength = 0.5, # maximal off-diagonal interaction strength
  c = 0.02 # connectance = interaction probability
){
  A = matrix(
    data = runif(N^2, min = min.strength, max = max.strength),
    nrow = N,
    ncol = N
  )
  diag(A) <- d
  c_obs <- (length(A[A!=0])-N)/(ncol(A)*ncol(A)-N)
  if(c_obs < c){
    edgeNumAdded = 0
    while(c_obs < c){
      # source node of edge
      xpos = sample(1:N, size = 1)
      # target node of edge
      ypos = sample(1:N, size = 1)
      # avoid diagonal
      if(xpos != ypos){
        if(A[xpos,ypos]==0){
          edgeNumAdded = edgeNumAdded +1
        }
        A[xpos,ypos] = 1
        c_obs = getConnectance(A)
      }
    }
  } else if (c_obs > c){
    edgeNumRemoved = 0
    while(c_obs > c){
      xpos = sample(1:N, size = 1)
      ypos = sample(1:N, size = 1)
      if(xpos != ypos){
        if(A[xpos, ypos] != 0){
          edgeNumRemoved = edgeNumRemoved +1
        }
        A[xpos,ypos]=0
        c_obs = getConnectance(A)
      }
    }
  }
  return(A)
}

getConnectance <- function(A){
  N <- nrow(A)
  c_obs <- (length(A[A!=0])-N)/(ncol(A)*ncol(A)-N)
  return(c_obs)
}

