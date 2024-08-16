#' Regional Heritability Mapping via Gibbs Sampling with Horseshoe Regularization
#'
#' This function implements Gibbs sampling for Bayesian estimation of random effects
#' in regional heritability mapping. It uses Horseshoe techniques for regularization
#' of random effects coefficients. The function handles the random effects by utilizing singular
#' value decomposition components of given covariance matrices derived from genomic sub-blocks.
#'
#' @param y Numeric vector, the response variable from the study.
#' @param V_list List of matrices, containing the left singular vectors of the covariance matrices.
#' @param D_list List of numeric vectors, containing the normalized singular values of the covariance matrices.
#' @param niter Integer, total number of Gibbs sampling iterations to perform.
#' @param thinin Integer, thinning interval to reduce autocorrelation in the sampled results.
#' @param burnin Integer, number of initial samples to discard to allow the chain to reach stability.
#' @param a Numeric, hyperparameter for the prior distribution of the residual variance.
#' @param b Numeric, hyperparameter for the prior distribution of the residual variance.
#' @param alpha Numeric, the significance level for the selection of random effects coefficients.
#' @param save_Z Logical, if TRUE, the function will save the random effects coefficients.
#'
#' @return A list containing sampled values post burn-in for:
#' \itemize{
#'   \item{\code{lambda}}{Matrix of random effects coefficients for each group and each iteration.}
#'   \item{\code{mu}}{Global mean effects across iterations.}
#'   \item{\code{sigma2}}{Global variance across iterations.}
#'   \item{\code{tau2}}{Global control of random effects across iterations.}
#'   \item{\code{W2}}{Local control of random effects across iterations.}
#'   \item{\code{Z}}{Random effects across iterations.}
#'   \item{\code{lambda2_median}}{List of numeric vectors, containing the median of random effects coefficients for each group.}
#' }
#'
#' @examples
#' dat = generate_test_data()
#' rhm_horseshoe(dat$y, dat$V_list, dat$D_list, niter = 2, thinin = 1, burnin = 0)
#'
#' @importFrom stats rnorm
#' @export
rhm_horseshoe <- function(y, V_list, D_list, niter=5000, thinin=5, burnin=2000,a=1,b=1,alpha=0.05,save_Z=FALSE) {
  k = length(V_list) # Number of groups of similarity matrices
  q = Reduce("+",lapply(D_list,function(D_g){length(D_g)})) # Number of matrices
  n = length(y) # Number of Observations
  n_result = (niter - burnin) %/% thinin
  Z = lapply(D_list,function(D_g){
    lapply(D_g,function(D_l){
      sqrt(D_l)*rnorm(length(D_l))
    })
  })
  mu = mean(y)
  lambda = lapply(Z, function(z) rep(0,length(z)) ) # Random Effects Coefficients
  sigma2 = mean(y^2) -(mean(y))^2 # Global Variance
  tau2 = 1  # Global Control Rate of random effects
  xi_tau = 2
  W2 = lapply(lambda, function(g) rep(1, length(g)))  # Local Control of Random effects
  xi_W = lapply(W2, function(W2g){rep(2,length(W2g))})
  lambda_result = lapply(lambda,function(g){
    matrix(0,n_result,length(g))
  })
  mu_result = numeric(n_result)
  sigma2_result = numeric(n_result)
  tau2_result = numeric(n_result)
  W2_result = lapply(W2,function(g){
    matrix(0,n_result,length(g))
  })
  if(save_Z){
    Z_result = lapply(Z,function(g){
      lapply(g,function(z){
        matrix(0,n_result,length(z))
      })
    })
  }

  residuals = y - mu
  # Main Gibbs sampling loop
  for (iter in 1:niter) {
    residuals = residuals + mu
    mu = simulate_mu(sigma2,residuals)
    residuals = residuals - mu
    # RESIDUALS, Z AND LAMBDA UPDATED BY CPP SIDE-EFFECT
    simulate_Z_lambda_horseshoe_cpp(y,Z,mu,lambda,W2,sigma2,tau2,V_list,D_list,residuals)
    W2 = simulate_W2_horseshoe(lambda, sigma2, tau2, xi_W, W2)
    tau2 = simulate_tau2_horseshoe(lambda, sigma2, W2, xi_tau,q)
    xi_W = lapply(W2,simulate_xi)
    xi_tau = simulate_xi(tau2)
    sigma2 = simulate_sigma2_horseshoe(y, Z, lambda, mu, tau2, W2,  a, b,V_list,q,residuals)

    # Store samples post burn-in and per thinning interval
    if (iter > burnin && ((iter - burnin) %% thinin == 0)) {
      index = (iter - burnin) %/% thinin
      for(g in 1:length(lambda)) {
        lambda_result[[g]][index,] = lambda[[g]]
        W2_result[[g]][index,] = W2[[g]]
        if(save_Z){
          for(l in 1:length(Z[[g]])){
            Z_result[[g]][[l]][index,] = Z[[g]][[l]]
          }
        }

      }
      mu_result[index] = mu
      sigma2_result[index] = sigma2
      tau2_result[index] = tau2
    }
  }

  # Return list of samples for analysis
  lambda2_median = lapply(lambda_result, function(lambda_g) apply(lambda_g^2, 2, median))
  result = list(lambda = lambda_result, mu = mu_result, sigma2 = sigma2_result, tau2 = tau2_result, W2 = W2_result,lambda2_median = lambda2_median)
  if(save_Z){
    result$Z = Z_result
  }
  return(result)
}

#' Regional Heritability Mapping via Gibbs Sampling with Group Fused Horseshoe Regularization
#'
#' This function implements Gibbs sampling for Bayesian estimation of random effects
#' in regional heritability mapping. It uses group-fused Horseshoe techniques for regularization
#' of random effects coefficients. The function handles the random effects by utilizing singular
#' value decomposition components of given covariance matrices derived from genomic sub-blocks.
#'
#' @param y Numeric vector, the response variable from the study.
#' @param V_list List of matrices, containing the left singular vectors of the covariance matrices.
#' @param D_list List of numeric vectors, containing the normalized singular values of the covariance matrices.
#' @param niter Integer, total number of Gibbs sampling iterations to perform.
#' @param thinin Integer, thinning interval to reduce autocorrelation in the sampled results.
#' @param burnin Integer, number of initial samples to discard to allow the chain to reach stability.
#' @param a Numeric, hyperparameter for the prior distribution of the residual variance.
#' @param b Numeric, hyperparameter for the prior distribution of the residual variance.
#' @param alpha Numeric, the significance level for the selection of random effects coefficients.
#' @param save_Z Logical, if TRUE, the function will save the random effects coefficients.
#'
#' @return A list containing sampled values post burn-in for:
#' \itemize{
#'   \item{\code{lambda}}{Matrix of random effects coefficients for each group and each iteration.}
#'   \item{\code{mu}}{Global mean effects across iterations.}
#'   \item{\code{sigma2}}{Global variance across iterations.}
#'   \item{\code{tau2}}{Global control of random effects across iterations.}
#'   \item{\code{Omega2}}{Local control of consecutive differences across iterations.}
#'   \item{\code{W2}}{Local control of random effects across iterations.}
#'   \item{\code{Z}}{Random effects across iterations.}
#'   \item{\code{lambda2_median}}{List of numeric vectors, containing the median of random effects coefficients for each group.}
#' }
#'
#' @examples
#' dat = generate_test_data()
#' rhm_fused(dat$y, dat$V_list, dat$D_list, niter = 2, thinin = 1, burnin = 0)
#'
#' @importFrom stats rnorm
#' @export
rhm_fused <- function(y, V_list, D_list, niter=5000, thinin=5, burnin=2000,a=1,b=1,alpha=0.05,save_Z=FALSE) {
  k = length(V_list) # Number of groups of similarity matrices
  q = Reduce("+",lapply(D_list,function(D_g){length(D_g)}))
  n = length(y) # Number of Observations
  n_result = (niter - burnin) %/% thinin
  Z = lapply(D_list,function(D_g){
    lapply(D_g,function(D_l){
      sqrt(D_l)*rnorm(length(D_l))
    })
  })
  mu = mean(y)
  lambda = lapply(Z, function(z) rep(0,length(z))) # Random Effects Coefficients
  sigma2 = mean(y^2) - mean(y)^2 # Global Variance
  tauO2 = 1
  xi_tauO = simulate_xi(tauO2)
  W2 = lapply(lambda, function(g) rep(1, length(g)))  # Local Control of Random effects
  xi_W = lapply(W2, simulate_xi)
  Omega2 = lapply(lambda, function(lambdag) rep(1,length(lambdag)-1))  # Local Control of Consecutive differences
  xi_Omega = lapply(Omega2, simulate_xi)
  lambda_result = lapply(lambda,function(g){
    matrix(0,n_result,length(g))
  })
  mu_result = numeric(n_result)
  sigma2_result = numeric(n_result)
  tauO2_result = numeric(n_result)
  Omega2_result = lapply(Omega2,function(g){
    matrix(0,n_result,length(g))
  })
  W2_result = lapply(W2,function(g){
    matrix(0,n_result,length(g))
  })
  if(save_Z){
    Z_result = lapply(Z,function(g){
      lapply(g,function(z){
        matrix(0,n_result,length(z))
      })
    })
  }
  residuals = y - mu
  # Main Gibbs sampling loop
  for (iter in 1:niter) {
    residuals = residuals + mu
    mu = simulate_mu(sigma2,residuals)
    residuals = residuals - mu
    # RESIDUALS, Z AND LAMBDA UPDATED BY CPP SIDE-EFFECT
    simulate_Z_lambda_fused_cpp(y,Z,mu,lambda,W2,Omega2,sigma2,tauO2,V_list,D_list,residuals)
    W2 = simulate_W2_fused(lambda, sigma2, xi_W, W2)
    Omega2 = simulate_Omega2(lambda, sigma2, tauO2, xi_Omega,Omega2)
    tauO2 = simulate_tauO2_fused(lambda, sigma2, Omega2, xi_tauO,q)
    xi_tauO = simulate_xi(tauO2)
    xi_Omega = lapply(Omega2,simulate_xi)
    xi_W = lapply(W2,simulate_xi)
    sigma2 = simulate_sigma2_fused(y, Z, lambda, mu,tauO2, W2, Omega2, a, b,V_list,q,residuals)
    # Store samples post burn-in and per thinning interval
    if (iter > burnin && ((iter - burnin) %% thinin == 0)) {
      index = (iter - burnin) %/% thinin
      for(g in 1:length(lambda)) {
        lambda_result[[g]][index,] = lambda[[g]]
        Omega2_result[[g]][index,] = Omega2[[g]]
        W2_result[[g]][index,] = W2[[g]]
        if(save_Z){
          for(l in 1:length(Z[[g]])){
            Z_result[[g]][[l]][index,] = Z[[g]][[l]]
          }
        }
      }
      mu_result[index] = mu
      sigma2_result[index] = sigma2
      tauO2_result[index] = tauO2
    }
  }

  # Return list of samples for analysis
  lambda2_median = lapply(lambda_result, function(lambda_g) apply(lambda_g^2, 2, median))
  result = list(lambda = lambda_result, mu = mu_result, sigma2 = sigma2_result,tauO2 = tauO2_result,W2=W2_result, Omega2 = Omega2_result,lambda2_median = lambda2_median)
  if(save_Z){
    result$Z = Z_result
  }
  return(result)
}

#' Regional Heritability Mapping via Gibbs Sampling with a Fusion Regularization Technique
#'
#' This function performs Gibbs sampling for Bayesian estimation of random effects, using a horseshoe fusion regularization approach.
#' It handles the random effects by utilizing singular value decomposition components of given covariance matrices derived from genomic sub-blocks.
#'
#' @param y Numeric vector, the response variable from the study.
#' @param V_list List of matrices, containing the left singular vectors of the covariance matrices.
#' @param D_list List of numeric vectors, containing the normalized singular values of the covariance matrices.
#' @param niter Integer, total number of Gibbs sampling iterations to perform.
#' @param thinin Integer, thinning interval to reduce autocorrelation in the sampled results.
#' @param burnin Integer, number of initial samples to discard to allow the chain to reach stability.
#' @param a Numeric, hyperparameter for the prior distribution of the residual variance.
#' @param b Numeric, hyperparameter for the prior distribution of the residual variance.
#' @param alpha Numeric, the significance level for the selection of random effects coefficients.
#' @param save_Z Logical, if TRUE, the function will save the random effects coefficients.
#'
#' @return A list containing sampled values post burn-in for:
#' \itemize{
#'   \item{\code{lambda}}{Matrix of random effects coefficients for each group and each iteration.}
#'   \item{\code{mu}}{Global mean effects across iterations.}
#'   \item{\code{sigma2}}{Global variance across iterations.}
#'   \item{\code{tau2}}{Global control of random effects across iterations.}
#'   \item{\code{Omega2}}{Local control of consecutive differences across iterations.}
#'   \item{\code{W2}}{Local control of random effects across iterations.}
#'   \item{\code{Z}}{Random effects across iterations.}
#'   \item{\code{lambda2_median}}{List of numeric vectors, containing the median of random effects coefficients for each group.}
#' }
#' @examples
#' dat = generate_test_data()
#' rhm_fusion(dat$y, dat$V_list, dat$D_list, niter = 2, thinin = 1, burnin = 0)
#'
#' @importFrom stats rnorm
#' @export
rhm_fusion <- function(y, V_list, D_list, niter=5000, thinin=5, burnin=2000,a=1,b=1,alpha=0.05,save_Z=FALSE) {
  k = length(V_list) # Number of groups of similarity matrices
  q = Reduce("+",lapply(D_list,function(D_g){length(D_g)}))
  n = length(y) # Number of Observations
  n_result = (niter - burnin) %/% thinin
  Z = lapply(D_list,function(D_g){
    lapply(D_g,function(D_l){
      sqrt(D_l)*rnorm(length(D_l))
    })
  })
  mu = mean(y)
  lambda = lapply(Z, function(z) rep(0,length(z))) # Random Effects Coefficients
  sigma2 = mean(y^2) - mean(y)^2  # Global Variance
  tauO2 = 1
  xi_tauO = simulate_xi(tauO2)
  Omega2 = lapply(lambda, function(g) rep(1,length(g)-1))  # Local Control of Consecutive differences
  W2 = lapply(lambda,function(g) rep(1,2)) # Local control of limit-conditions for each group
  xi_Omega = lapply(Omega2, simulate_xi)
  xi_W = lapply(W2,simulate_xi)
  lambda_result = lapply(lambda,function(g){
    matrix(0,n_result,length(g))
  })
  mu_result = numeric(n_result)
  sigma2_result = numeric(n_result)
  tauO2_result = numeric(n_result)
  Omega2_result = lapply(Omega2,function(g){
    matrix(0,n_result,length(g))
  })
  W2_result= lapply(W2,function(g){
    matrix(0,n_result,2)
  })
  if(save_Z){
    Z_result = lapply(Z,function(g){
      lapply(g,function(z){
        matrix(0,n_result,length(z))
      })
    })
  }
  residuals = y - mu
  # Main Gibbs sampling loop
  for (iter in 1:niter) {
    residuals = residuals + mu
    mu = simulate_mu(sigma2,residuals)
    residuals = residuals - mu
    # RESIDUALS, Z AND LAMBDA UPDATED BY CPP SIDE-EFFECT
    simulate_Z_lambda_fusion_cpp(y,Z,mu,lambda,W2,Omega2,sigma2,tauO2,V_list,D_list,residuals)
    W2 = simulate_W2_fusion(lambda, sigma2, xi_W, W2)
    Omega2 = simulate_Omega2(lambda, sigma2, tauO2, xi_Omega,Omega2)
    tauO2 = simulate_tauO2_fusion(lambda, sigma2, Omega2, xi_tauO,q)
    xi_tauO = simulate_xi(tauO2)
    xi_W = lapply(W2,simulate_xi)
    xi_Omega = lapply(Omega2,simulate_xi)
    sigma2 = simulate_sigma2_fusion(y, Z, lambda, mu,tauO2,W2,Omega2, a, b,V_list,q,residuals)
    # Store samples post burn-in and per thinning interval
    if (iter > burnin && ((iter - burnin) %% thinin == 0)) {
      index = (iter - burnin) %/% thinin
      for(g in 1:length(lambda)) {
        lambda_result[[g]][index,] = lambda[[g]]
        Omega2_result[[g]][index,] = Omega2[[g]]
        W2_result[[g]][index,] = W2[[g]]
        if(save_Z){
            for(l in 1:length(Z[[g]])){
            Z_result[[g]][[l]][index,] = Z[[g]][[l]]
          }
        }
      }
      mu_result[index] = mu
      sigma2_result[index] = sigma2
      tauO2_result[index] = tauO2
    }
  }

  # Return list of samples for analysis
  lambda2_median = lapply(lambda_result, function(lambda_g) apply(lambda_g^2, 2, median))
  result = list(lambda = lambda_result, mu = mu_result, sigma2 = sigma2_result,tauO2 = tauO2_result,W2=W2_result, Omega2 = Omega2_result,lambda2_median = lambda2_median)
  if(save_Z){
    result$Z = Z_result
  }
  return(result)
}



#' Fast Regional Heritability Mapping Analysis
#'
#' Performs Regional Heritability Mapping (RHM) analysis on a list of individual similarity matrices. Each similarity matrix is represented by its Singular Value Decomposition (SVD).
#'
#' @param y A numeric vector representing the phenotype data.
#' @param V_list A list of matrices, where each matrix contains the left singular vectors from the SVD of a genomic data subblock.
#' @param D_list A list of numeric vectors, where each vector contains the singular values from the SVD of a genomic data subblock.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{R2}}{A matrix containing the R-squared values for each similarity matrix and number of axis considered.}
#'   \item{\code{R2_weighted}}{A matrix containing the weighted R-squared values for each similarity matrix and number of axis considered.}
#' }
#'
#' @details
#' The function computes the R-squared values and weighted R-squared values for each combination of similarity matrices and number of axis to consider.
#' The phenotype data is first scaled and centered before performing the calculations.
#'
#' @examples
#' dat = generate_test_data()
#' result = fast_rhm(dat$y, dat$V_list, dat$D_list)
#' R2 = result$R2
#' R2_weighted = result$R2_weighted
#'
#' @export
fast_rhm = function(y,V_list,D_list) {
  n = length(y)
  q = Reduce("+",lapply(D_list,function(D_g){length(D_g)}))
  p = max(sapply(D_list,length))
  R2 = matrix(0,q,p)
  R2_weighted = matrix(0,q,p)
  V_list = Reduce(c,V_list)
  D_list = Reduce(c,D_list)
  y = scale(y)*sqrt(n/(n-1))
  for(i in 1:q) {
    for(j in 1:min(p,length(D_list[[i]]))) {
      R2[i,j] = tcrossprod(t(y)%*%V_list[[i]][,1:j])/n
      R2_weighted[i,j] = tcrossprod((t(y)%*%V_list[[i]][,1:j])*sqrt(D_list[[i]][1:j]))/n
    }
  }
  list(R2=R2,R2_weighted=R2_weighted)
}

#' Regional Heritability Mapping (RHM) Wrapper Function
#'
#' This function performs Regional Heritability Mapping (RHM) on genomic data using various regularization techniques.
#' It supports Horseshoe, Fused Horseshoe, Fusion, and a fast R-squared based method. The function processes the genomic
#' data into similarity matrices and then runs the selected RHM method.
#'
#' @param y Numeric vector, the response variable (phenotype) from the study.
#' @param X Matrix, the genotype data where rows represent individuals and columns represent SNPs.
#' @param chromosome_list Numeric or character vector, indicating the chromosome each SNP belongs to.
#' @param position_list Numeric vector, indicating the position of each SNP on the chromosome.
#' @param block_size Integer, the size of the genomic blocks to consider for similarity matrix computation.
#' @param niter Null or Integer, total number of Gibbs sampling iterations to perform.
#' @param thinin Null or Integer, thinning interval to reduce autocorrelation in the sampled results.
#' @param burnin Null or Integer, number of initial samples to discard to allow the chain to reach stability.
#' @param similarity_method Character, method to compute similarity matrices, must be one of "yang", "van_raden", "simple".
#' @param regularization_method Character, method for regularization in RHM, must be one of "horseshoe", "fused", "fusion", "fast","horseshoe_max","fused_max","fusion_max".
#' @param svd_inertia Numeric, the proportion of inertia to keep during Singular Value Decomposition (SVD).
#' @param verbose Logical, if TRUE, the function will print progress messages.
#' @param alpha Numeric, the significance level for the selection of random effects coefficients.
#' @param save_Z Logical, if TRUE, the function will save the random effects coefficients.
#'
#' @return A list containing the results from the selected RHM method. The structure of the list depends on the method used:
#' \describe{
#'   \item{\code{horseshoe}}{A list with elements lambda, mu, sigma2, tau2, W2, selected_lambdas, lambda2_median.}
#'   \item{\code{fused}}{A list with elements lambda, mu, sigma2, tauO2, W2, Omega2, selected_lambdas, lambda2_median}
#'   \item{\code{fusion}}{A list with elements lambda, mu, sigma2, tauO2, W2, Omega2, selected_lambdas, lambda2_median.}
#'   \item{\code{fast}}{A list with elements R2, R2_weighted.}
#' }
#' @examples
#' dat = generate_snp_data(by_chr=FALSE)
#' result = rhm(dat$y, dat$X, dat$chromosome_list, dat$position_list,30,100,5,10)
#'
#' @importFrom stats rnorm
#' @export
rhm <- function(y, X, chromosome_list, position_list, block_size,niter=NULL,burnin=NULL,thinin=NULL, similarity_method = "simple", regularization_method = "fused", svd_inertia = 0.95, verbose = TRUE,alpha=0.05,save_Z=FALSE) {
  SIMILARITY_METHODS = c("yang","van_raden","simple")
  REGULARIZATION_METHODS = c("horseshoe","fused","fusion","fast")
  if (!(similarity_method %in% SIMILARITY_METHODS)) stop("Invalid similarity method")
  if (!(regularization_method %in% REGULARIZATION_METHODS)) stop("Invalid RHM method")
  if (regularization_method %in%REGULARIZATION_METHODS) {
    if(is.null(niter)){niter=5000}
    if(is.null(burnin)){burnin=2000}
    if(is.null(thinin)){thinin=5}
  }
  X <- lapply(unique(chromosome_list), function(chr) {
    X[, chromosome_list == chr]
  })
  position_list <- lapply(unique(chromosome_list), function(chr) {
    position_list[chromosome_list == chr]
  })

  if (verbose) print("Cutting SNP matrices...")
  X <- cut_snp_matrices(X, position_list, block_size, block_size / 2)

  if (verbose) print("Computing similarity matrices...")
  duration_similarity_computing = system.time({
    A_list <- get_snp_similarity_matrices(X, svd_inertia, similarity_method)
    V_list <- lapply(A_list, function(A_g) {
      lapply(A_g, function(A) A$V)
    })
    D_list <- lapply(A_list, function(A_g) {
      lapply(A_g, function(A) A$D)
    })
  })


  if (verbose) print("Running RHM...")
  duration = system.time({
    result = switch(regularization_method,
                    "horseshoe" = rhm_horseshoe(y, V_list, D_list,niter,thinin,burnin ,alpha=alpha),
                    "fused" = rhm_fused(y, V_list, D_list,niter,thinin,burnin  ,alpha=alpha),
                    "fusion" = rhm_fusion(y, V_list, D_list,niter,thinin,burnin  ,alpha=alpha),
                    "fast" = fast_rhm(y, V_list, D_list)
    )
  })
  result$duration = duration
  if(verbose) print("Done!")
  return(result)
}


