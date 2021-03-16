#' @title Interaction matrix with Power-Law network adjacency matrix
#' @description Generate an interaction matrix A that can be decomposed as NH.Gs .
#' Where N is the Normal Interspecific Interaction matrix, H the interaction strength heterogeneity
#' drawn from a power-law distribution with given parameter alpha, and G the adjacency matrix of
#' power-law out-degree digraph ecological network, and s a scaling factor. Diagonal elements of A are subsequently set to -1.
#' @references Gibson TE, Bashan A, Cao HT, Weiss ST, Liu YY (2016)
#' On the Origins and Control of Community Types in the Human Microbiome.
#' PLOS Computational Biology 12(2): e1004688. https://doi.org/10.1371/journal.pcbi.1004688
#' @param n the number of species
#' @param alpha the power-law distribution parameter. Should be > 1.
#' Larger values will give lower interaction strength heterogeneity,
#' whereas values closer to 1 give strong heterogeneity in interaction strengths between the species.
#' In other words, values of alpha close to 1 will give Strongly Interacting Species (SIS).
#' @param stdev the standard deviation parameter of the normal distribution with mean 0 from which
#' the elements of the nominal interspecific interaction matrix N are drawn
#' @param s scaling parameter with which the final global interaction matrix A is multiplied
#' @return The global interaction matrix A with n rows and n columns.
#' @examples
#' # Low interaction heterogeneity
#' A_low <- powerlawA(n = 10, alpha = 3)
#' # Strong interaction heterogeneity
#' A_strong <- powerlawA(n = 10, alpha = 1.01)
#' @export

powerlawA <- function(
  n, # number of species
  alpha, # power-law distribution parameter
  stdev = 1, # sd normal distribution
  s = 1 # scaling parameter
){
  # Nominal Interspecific Interaction matrix N
  N <- matrix(
    data = rnorm(n^2, mean = 0, sd = stdev),
    nrow = n,
    ncol = n
  )
  diag(N) <- 0

  # Interaction strength heterogeneity H
  H <- rplcon(n = n, xmin = 1, alpha = alpha)
  H <- H / mean(H)
  H <- diag(H)

  # Adjacency matrix G of power-law out-degree digraph ecological network
  h_bar <- rplcon(n = n, xmin = 1, alpha = alpha)
  h <- sapply(1:n, FUN = function(i){
    min(
      ceiling(h_bar[i]/mean(h_bar)),
      n
    )
  })
  G <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    index <- sample(x = 1:n, size = h[i])
    G[index, i] <- i
  }

  A <- N %*% H * G * s
  diag(A) <- -1
  colnames(A) <- 1:n
  rownames(A) <- 1:n
  return(A)
}
