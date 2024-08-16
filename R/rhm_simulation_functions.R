
#' Simulate Global Mean
#'
#' Simulate global mean following the conditional complete posterior distribution.
#'
#' @param y Response vector.
#' @param Z List of matrices representing random effects.
#' @param lambda List of coefficients.
#' @param mu Mean value.
#' @param sigma2 Residual variance.
#' @param V_list List of matrices for variance calculations.
#' @param residuals Residuals.
#'
#' @return A simulated mean value.
#' @keywords internal
simulate_mu = function(sigma2,residuals){
  n = length(residuals)
  center = mean(residuals)
  sd  = sqrt(sigma2/n)
  return(rnorm(1,center,sd))
}

#' Simulate Control Variable of Horseshoe Prior
#'
#' Simulate control variable of the Horseshoe prior.
#'
#' @param x2 Squared coefficients.
#'
#' @return A simulated control variable.
#' @keywords internal
simulate_xi = function(x2) {
  n = length(x2)
  a = 1
  b = 1+ 1/x2
  return(1 / rgamma(n,a,b))
}

#' Simulate Local Control of Coefficients(Horseshoe)
#'
#' Simulate local control of coefficients following the conditional complete posterior distribution (fused version).
#'
#' @param lambda List of coefficients.
#' @param sigma2 Residual variance.
#' @param tau2 Global variance parameter.
#' @param xi_W Control variable of the Horseshoe prior.
#' @param W2 List of local variance parameters.
#'
#' @return A list of updated local control variables for coefficients.
#' @importFrom stats rgamma
#' @keywords internal
simulate_W2_horseshoe = function(lambda,sigma2,tau2,xi_W,W2){
  for(g in 1:length(W2)) {
    n = length(W2[[g]])
    a = 1
    b = 1/xi_W[[g]] + lambda[[g]]^2/(2*sigma2*tau2)
    W2[[g]] = 1 / rgamma(n,a,b)
  }
  return(W2)
}

#' Simulate Local Control of Coefficients(Fused)
#'
#' Simulate local control of coefficients following the conditional complete posterior distribution (fused version).
#'
#' @param lambda List of coefficients.
#' @param sigma2 Residual variance.
#' @param xi_W Control variable of the Horseshoe prior.
#' @param W2 List of local shrinkage parameters.
#'
#' @return A list of updated local control variables for coefficients.
#' @importFrom stats rgamma
#' @keywords internal
simulate_W2_fused = function(lambda,sigma2,xi_W,W2){
  for(g in 1:length(W2)) {
    n = length(W2[[g]])
    a = 1
    b = 1/xi_W[[g]] + lambda[[g]]^2/(2*sigma2)
    W2[[g]] = 1 / rgamma(n,a,b)
  }
  return(W2)
}

#' Simulate Local Control of Coefficients
#'
#' Simulate local control of coefficients following the conditional complete posterior distribution (fused version).
#'
#' @param lambda List of coefficients.
#' @param sigma2 Residual variance.
#' @param xi_W Control variable of the Horseshoe prior.
#' @param W2 List of local variance parameters.
#'
#' @return A list of updated local control variables for coefficients.
#' @importFrom stats rgamma
#' @keywords internal
simulate_W2_fusion = function(lambda,sigma2,xi_W,W2){
  for(g in 1:length(W2)) {
    n = length(lambda[[g]])
    a = 1
    b = c(lambda[[g]][1]^2,lambda[[g]][n]^2)/(2*sigma2) + 1/xi_W[[g]]
    W2[[g]] = 1 / rgamma(2,a,b)
  }
  return(W2)
}


#' Simulate Residual Variance (horseshoe)
#'
#' Simulate residual variance following the conditional complete posterior distribution (horseshoe version).
#'
#' @param y Response vector.
#' @param Z List of matrices representing random effects.
#' @param lambda List of coefficients.
#' @param mu Mean value.
#' @param tau2 Global variance parameter.
#' @param W2 List of local variance parameters.
#' @param a Shape parameter for the prior distribution.
#' @param b Rate parameter for the prior distribution.
#' @param V_list List of matrices for variance calculations.
#' @param q Degrees of freedom.
#' @param residuals Residuals.
#'
#' @return A simulated residual variance value.
#' @keywords internal
simulate_sigma2_horseshoe = function(y,Z,lambda,mu,tau2,W2,a,b,V_list,q,residuals){
  a = a + q/2 + length(y)/2
  b = b + sum(residuals^2)/2 + (1/(2*tau2))*Reduce("+",lapply(1:length(lambda),function(g){
    sum(lambda[[g]]^2/W2[[g]])
  }))
  return(1 / rgamma(1,a,b))
}

#' Simulate Residual Variance (Fusion)
#'
#' Simulate residual variance following the conditional complete posterior distribution.
#'
#' @param y Response vector.
#' @param Z List of matrices representing random effects.
#' @param lambda List of coefficients.
#' @param mu Mean value.
#' @param tauO2 Global shrinkage parameter for differences.
#' @param W2 List of local shrinkage parameters for limit-conditions of each group.
#' @param Omega2 List of local shrinkage parameters for differences.
#' @param a Shape parameter for the prior distribution.
#' @param b Rate parameter for the prior distribution.
#' @param V_list List of matrices for variance calculations.
#' @param q Degrees of freedom.
#' @param residuals Residuals.
#'
#' @return A simulated residual variance value.
#' @keywords internal
simulate_sigma2_fusion <- function(y, Z, lambda, mu,tauO2,W2, Omega2, a, b, V_list, q,residuals) {
  G = length(lambda)
  a = a + q/2 + G/2 + length(y)/2
  b = b + sum(residuals^2) / 2 + (1/2) * Reduce("+", lapply(1:length(lambda), function(g) {
    n = length(lambda[[g]])
    sum(diff((lambda[[g]]))^2 / Omega2[[g]]/tauO2) + sum(c(lambda[[g]][1]^2,lambda[[g]][n]^2)/W2[[g]])
  }))
  return(1 / rgamma(1, a, b))
}


#' Simulate Residual Variance (Fused)
#'
#' Simulate residual variance following the conditional complete posterior distribution (fused version).
#'
#' @param y Response vector.
#' @param Z List of matrices representing random effects.
#' @param lambda List of coefficients.
#' @param mu Mean value.
#' @param tauO2 Global shrinkage parameter for differences.
#' @param W2 List of local shrinkage parameters.
#' @param Omega2 List of local shrinkage parameters for differences.
#' @param a Shape parameter for the prior distribution.
#' @param b Rate parameter for the prior distribution.
#' @param V_list List of matrices for variance calculations.
#' @param q Degrees of freedom.
#' @param residuals Residuals.
#'
#' @return A simulated residual variance value.
#' @keywords internal
simulate_sigma2_fused = function(y,Z,lambda,mu,tauO2,W2,Omega2,a,b,V_list,q,residuals){
  G = length(lambda)
  a = a + q - G/2 + length(y)/2
  b = b + sum(residuals^2)/2 + (1/2)*Reduce("+",lapply(1:length(lambda),function(g){
    sum(lambda[[g]]^2/W2[[g]]) + sum(diff((lambda[[g]]))^2/Omega2[[g]])/tauO2
  }))
  return(1 / rgamma(1,a,b))
}


#' Simulate Global Control of Coefficients (Horseshoe)
#'
#' Simulate global control of coefficients following the conditional complete posterior distribution (fused version).
#'
#' @param lambda List of coefficients.
#' @param sigma2 Residual variance.
#' @param W2 List of local variance parameters.
#' @param xi_tau Control variable of the Horseshoe prior.
#' @param q Degrees of freedom.
#'
#' @return A simulated global control variable.
#' @keywords internal
simulate_tau2_horseshoe = function(lambda,sigma2,W2,xi_tau,q){
  a = (q+1)/2
  b = 1/xi_tau + (1/(2*sigma2)) * Reduce("+",lapply(1:length(lambda),function(g){
    sum(lambda[[g]]^2/W2[[g]])
  }))
  return(1 / rgamma(1,a,b))
}

#' Simulate Global Control of Coefficients(Fusion)
#'
#' Simulate global control of coefficients following the conditional complete posterior distribution.
#'
#' @param lambda List of coefficients.
#' @param sigma2 Residual variance.
#' @param Omega2 List of local shrinkage parameters for differences.
#' @param xi_tauO Control variable of the Horseshoe prior.
#' @param q Degrees of freedom.
#'
#' @return A simulated global control variable.
#' @keywords internal
simulate_tauO2_fusion = function(lambda,sigma2,Omega2,xi_tauO,q){
  G = length(lambda)
  a = (q + 1 - G)/2
  b = 1/xi_tauO + (1/(2*sigma2)) * Reduce("+",lapply(1:length(lambda),function(g){
    sum(diff((lambda[[g]]))^2/Omega2[[g]])
  }))
  return(1 / rgamma(1,a,b))
}


#' Simulate global shrinkage parameter of differences (Fused)
#'
#' Simulate global shrinkage parameters of difference following the conditional complete posterior distribution (fused version).
#'
#' @param lambda List of coefficients.
#' @param sigma2 Residual variance.
#' @param Omega2 List of local shrinkage parameters for differences.
#' @param xi_tauO Control variable of the Horseshoe prior.
#' @param q Degrees of freedom.
#'
#' @return A simulated global control variable.
#' @keywords internal
simulate_tauO2_fused = function(lambda,sigma2,Omega2,xi_tauO,q){
  G = length(lambda)
  a = (q+1-G)/2
  b = 1/xi_tauO + (1/(2*sigma2)) * Reduce("+",lapply(1:length(lambda),function(g){
    sum(diff((lambda[[g]]))^2/Omega2[[g]])
  }))
  return(1 / rgamma(1,a,b))
}


#' Simulate Local Control of Differences
#'
#' Simulate local control of differences following the conditional complete posterior distribution (fused version).
#'
#' @param lambda List of coefficients.
#' @param sigma2 Residual variance.
#' @param tauO2 Global shrinkage parameter of differences.
#' @param xi_Omega Control variable of the Horseshoe prior.
#' @param Omega2 List of local shrinkage parameters for differences.
#'
#' @return A list of updated local control variables for differences.
#' @keywords internal
simulate_Omega2 = function(lambda,sigma2,tauO2,xi_Omega,Omega2){
  for(g in 1:length(Omega2)) {
    n = length(Omega2[[g]])
    a = 1
    b = 1/xi_Omega[[g]] + diff((lambda[[g]]))^2/(2*sigma2*tauO2)
    Omega2[[g]] = 1 / rgamma(n,a,b)
  }
  return(Omega2)
}
