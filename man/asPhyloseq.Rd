\name{asPhyloseq}
\alias{asPhyloseq}
\title{
Microbiome Phyloseq wrapper for community samples
}
\description{
Takes in samples in vector or matrix format and merges them to a microbiome phyloseq class object.
}
\usage{
asPhyloseq(S, ...)
}

\arguments{
  \item{S, ...}{
Each parameter given to the function is a sample vector or matrix to be merged into one phyloseq object of samples. The vectors/matrices must all contain the same number of species in order to be merged. A phyloseq object requires at least a matrix of min. 2 columns, therefore an error is thrown if only a vector is given without additional vector/matrix objects.
}
}

\value{
A microbiome phyloseq class object with species/OTUs as rows and samples as columns.
}

\author{
Gheysen E.
}

\seealso{
microbiome R package
}
\examples{
# all samples normalised
glv_samples <- cbind(rmse_sample(glv(N = 8, A = generateA(8), tend = 10000)),
    rmse_sample(glv(N = 8, A = generateA(8, c = 0.3), tend = 5000))) # matrix
hubbell_sample <- hubbell(N = 4, M = 8, I = 1000, d = 10, m = 0.02, norm = TRUE, tend = 800)[,800] # vector
soi_sample <- soi(N = 8, I = 1000, tend = 5000, norm = TRUE)[,5000] # vector

physeq_samples <- asPhyloseq(glv_samples, hubbell_sample, soi_sample)
}

