% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rhm_wrapper_functions.R
\name{rhm}
\alias{rhm}
\title{Regional Heritability Mapping (RHM) Wrapper Function}
\usage{
rhm(
  y,
  X,
  chromosome_list,
  position_list,
  block_size,
  niter = NULL,
  burnin = NULL,
  thinin = NULL,
  similarity_method = "simple",
  regularization_method = "fused",
  svd_inertia = 0.95,
  verbose = TRUE,
  alpha = 0.05,
  save_Z = FALSE
)
}
\arguments{
\item{y}{Numeric vector, the response variable (phenotype) from the study.}

\item{X}{Matrix, the genotype data where rows represent individuals and columns represent SNPs.}

\item{chromosome_list}{Numeric or character vector, indicating the chromosome each SNP belongs to.}

\item{position_list}{Numeric vector, indicating the position of each SNP on the chromosome.}

\item{block_size}{Integer, the size of the genomic blocks to consider for similarity matrix computation.}

\item{niter}{Null or Integer, total number of Gibbs sampling iterations to perform.}

\item{burnin}{Null or Integer, number of initial samples to discard to allow the chain to reach stability.}

\item{thinin}{Null or Integer, thinning interval to reduce autocorrelation in the sampled results.}

\item{similarity_method}{Character, method to compute similarity matrices, must be one of "yang", "van_raden", "simple".}

\item{regularization_method}{Character, method for regularization in RHM, must be one of "horseshoe", "fused", "fusion", "fast".}

\item{svd_inertia}{Numeric, the proportion of inertia to keep during Singular Value Decomposition (SVD).}

\item{verbose}{Logical, if TRUE, the function will print progress messages.}

\item{alpha}{Numeric, the significance level for the selection of random effects coefficients.}

\item{save_Z}{Logical, if TRUE, the function will save the random effects coefficients.}
}
\value{
A list containing the results from the selected RHM method. The structure of the list depends on the method used:
\describe{
  \item{\code{horseshoe}}{A list with elements lambda, mu, sigma2, tau2, W2, selected_lambdas, lambda2_median.}
  \item{\code{fused}}{A list with elements lambda, mu, sigma2, tauO2, W2, Omega2, selected_lambdas, lambda2_median}
  \item{\code{fusion}}{A list with elements lambda, mu, sigma2, tauO2, W2, Omega2, selected_lambdas, lambda2_median.}
  \item{\code{fast}}{A list with elements R2, R2_weighted.}
}
}
\description{
This function performs Regional Heritability Mapping (RHM) on genomic data using various regularization techniques.
It supports Horseshoe, Fused Horseshoe, Fusion, and a fast R-squared based method. The function processes the genomic
data into similarity matrices and then runs the selected RHM method.
}
\examples{
dat = generate_snp_data(by_chr=FALSE)
result = rhm(dat$y, dat$X, dat$chromosome_list, dat$position_list,30,100,5,10)

}
