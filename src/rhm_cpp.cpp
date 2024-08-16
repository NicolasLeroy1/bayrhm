#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' @name simulate_Z_lambda_horseshoe_cpp
//' @title Simulate Random Effects and Coefficients using Rcpp and Armadillo
//' @description This function simulates random effects and coefficients following the conditional complete posterior distribution using Rcpp and Armadillo for performance improvements.
//' @param y A numeric vector representing the response variable.
//' @param Z A list of lists of matrices representing random effects.
//' @param mu A double representing the mean value.
//' @param lambda A list of vectors representing coefficients.
//' @param W2 A list of vectors representing the local control of random effects.
//' @param sigma2 A double representing the residual variance.
//' @param tau2 A double representing the global control of random effects.
//' @param V_list A list of lists of left singular matrices used for variance calculations.
//' @param D_list A list of lists of matrices for the prior distribution.
//' @param residuals A numeric vector representing the residuals.
//' @return An integer value.
// [[Rcpp::export]]
int simulate_Z_lambda_horseshoe_cpp(const arma::vec& y, Rcpp::List& Z, double mu,
                                    Rcpp::List& lambda, const Rcpp::List& W2, double sigma2,
                                    double tau2, const Rcpp::List& V_list,const Rcpp::List& D_list, arma::vec& residuals) {
  double temp_var_inv_lambda;
  double center_lambda;
  double sd_lambda;
  for (int g = 0; g < Z.size(); ++g) {
    Rcpp::List Z_g = Z[g];
    Rcpp::List V_g = V_list[g];
    Rcpp::List D_g = D_list[g];
    arma::vec lambda_g = as<arma::vec>(lambda[g]);
    arma::vec W2_g = as<arma::vec>(W2[g]);
    for (int l = 0; l < Z_g.size(); ++l) {
      //INITIALIZE
      arma::mat V_gl = as<arma::mat>(V_g[l]);
      arma::vec Z_gl = as<arma::vec>(Z_g[l]);
      arma::vec D_gl = as<arma::vec>(D_g[l]);
      arma::vec u_l = V_gl * Z_gl;
      residuals += u_l * lambda_g[l];
      //SIMULATE LAMBDA
      temp_var_inv_lambda = (arma::dot(Z_gl, Z_gl) + 1.0 / W2_g[l] / tau2);
      center_lambda = arma::dot(residuals, u_l) / temp_var_inv_lambda;
      sd_lambda = std::sqrt(sigma2 / temp_var_inv_lambda);
      lambda_g[l] = R::rnorm(center_lambda, sd_lambda);
      //SIMULATE Z
      arma::vec temp_var_inv = ((sigma2/ D_gl) + std::pow(lambda_g[l],2));
      arma::vec center = lambda_g[l]* (V_gl.t() *  residuals)/temp_var_inv;
      arma::vec sd = arma::sqrt(sigma2 / temp_var_inv);
      arma::vec new_Z_gl = (arma::randn(Z_gl.n_elem) % sd) + center;
      Z_g[l] = Rcpp::wrap(new_Z_gl);
      residuals -= V_gl * new_Z_gl * lambda_g[l];
    }
    lambda[g] = lambda_g;
    Z[g] = Z_g;
  }
  return 0;
}

//' @name simulate_Z_lambda_fusion_cpp
//' @title Simulate Random Effects and Coefficients using Rcpp and Armadillo
//' @description This function simulates random effects and coefficients following the conditional complete posterior distribution using Rcpp and Armadillo for performance improvements.
//' @param y A numeric vector representing the response variable.
//' @param Z A list of lists of matrices representing random effects.
//' @param mu A double representing the mean value.
//' @param lambda A list of vectors representing coefficients.
//' @param W2 A list of vectors representing the local control of random effects limit-conditions of each group.
//' @param Omega2 A list of vectors representing the local control of random effects differences.
//' @param sigma2 A double representing the residual variance.
//' @param tauO2 A double representing the global control of random effects differences.
//' @param V_list A list of lists of left singular matrices used for variance calculations.
//' @param D_list A list of lists of matrices for the prior distribution.
//' @param residuals A numeric vector representing the residuals.
//' @return An integer value.
// [[Rcpp::export]]
int simulate_Z_lambda_fusion_cpp(const arma::vec& y, Rcpp::List& Z, double mu,
                                    Rcpp::List& lambda,const Rcpp::List& W2, const Rcpp::List& Omega2, double sigma2,
                                    double tauO2, const Rcpp::List& V_list,const Rcpp::List& D_list, arma::vec& residuals) {
  double temp_var_inv_lambda;
  double center_lambda;
  double sd_lambda;
  double znorm;
  for (int g = 0; g < Z.size(); ++g) {
    Rcpp::List Z_g = Z[g];
    Rcpp::List V_g = V_list[g];
    Rcpp::List D_g = D_list[g];
    arma::vec lambda_g = as<arma::vec>(lambda[g]);
    arma::vec Omega2_g = as<arma::vec>(Omega2[g]);
    arma::vec W2_g = as<arma::vec>(W2[g]);
    for (int l = 0; l < Z_g.size(); ++l) {
      //INITIALIZE
      arma::mat V_gl = as<arma::mat>(V_g[l]);
      arma::vec Z_gl = as<arma::vec>(Z_g[l]);
      arma::vec D_gl = as<arma::vec>(D_g[l]);
      arma::vec u_l = V_gl * Z_gl;
      znorm = arma::dot(Z_gl, Z_gl);
      residuals += u_l * lambda_g[l];
      //SIMULATE LAMBDA
      if(l==0){temp_var_inv_lambda =  znorm + (1.0 / Omega2_g[l]/ tauO2) + (1.0/W2_g[0]);}
      else if(l==(Z_g.size()-1)) {temp_var_inv_lambda = (1.0/Omega2_g[l-1]/ tauO2) + (1.0/W2_g[1]) + znorm  ;}
      else {temp_var_inv_lambda = ((1.0 / Omega2_g[l] + 1.0/Omega2_g[l-1] )/ tauO2) + znorm;}
      center_lambda = arma::dot(residuals, u_l);
      if(l==0){center_lambda = (center_lambda+lambda_g[l+1]/Omega2_g[l]/tauO2)/temp_var_inv_lambda;}
      else if(l==(Z_g.size()-1)) {center_lambda = (center_lambda+lambda_g[l-1]/Omega2_g[l-1]/tauO2)/temp_var_inv_lambda;}
      else {center_lambda = (center_lambda+lambda_g[l-1]/Omega2_g[l-1]/tauO2 + (lambda_g[l+1])/Omega2_g[l]/tauO2)/temp_var_inv_lambda;}
      sd_lambda = std::sqrt(sigma2 / temp_var_inv_lambda);
      lambda_g[l] = R::rnorm(center_lambda, sd_lambda);
      //SIMULATE Z
      arma::vec temp_var_inv = ((sigma2/ D_gl) + std::pow(lambda_g[l],2));
      arma::vec center = lambda_g[l]* (V_gl.t() *  residuals)/temp_var_inv;
      arma::vec sd = arma::sqrt(sigma2 / temp_var_inv);
      arma::vec new_Z_gl = (arma::randn(Z_gl.n_elem) % sd) + center;
      Z_g[l] = Rcpp::wrap(new_Z_gl);
      residuals -= V_gl * new_Z_gl * lambda_g[l];
    }
    lambda[g] = lambda_g;
    Z[g] = Z_g;
  }
  return 0;
}


//' @name simulate_Z_lambda_fused_cpp
//' @title Simulate Random Effects and Coefficients using Rcpp and Armadillo
//' @description This function simulates random effects and coefficients following the conditional complete posterior distribution using Rcpp and Armadillo for performance improvements.
//' @param y A numeric vector representing the response variable.
//' @param Z A list of lists of matrices representing random effects.
//' @param mu A double representing the mean value.
//' @param lambda A list of vectors representing coefficients.
//' @param W2 A list of vectors representing the local control of random effects.
//' @param Omega2 A list of vectors representing the local control of random effects differences.
//' @param sigma2 A double representing the residual variance.
//' @param tauO2 A double representing the global control of random effects differences.
//' @param V_list A list of lists of left singular matrices used for variance calculations.
//' @param D_list A list of lists of matrices for the prior distribution.
//' @param residuals A numeric vector representing the residuals.
//' @return An integer value.
// [[Rcpp::export]]
int simulate_Z_lambda_fused_cpp(const arma::vec& y, Rcpp::List& Z, double mu,
                                 Rcpp::List& lambda,const Rcpp::List& W2, const Rcpp::List& Omega2, double sigma2,
                                 double tauO2, const Rcpp::List& V_list,const Rcpp::List& D_list, arma::vec& residuals) {
  double temp_var_inv_lambda;
  double center_lambda;
  double sd_lambda;
  double znorm;
  for (int g = 0; g < Z.size(); ++g) {
    Rcpp::List Z_g = Z[g];
    Rcpp::List V_g = V_list[g];
    Rcpp::List D_g = D_list[g];
    arma::vec lambda_g = as<arma::vec>(lambda[g]);
    arma::vec W2_g = as<arma::vec>(W2[g]);
    arma::vec Omega2_g = as<arma::vec>(Omega2[g]);
    for (int l = 0; l < Z_g.size(); ++l) {
      //INITIALIZE
      arma::mat V_gl = as<arma::mat>(V_g[l]);
      arma::vec Z_gl = as<arma::vec>(Z_g[l]);
      arma::vec D_gl = as<arma::vec>(D_g[l]);
      arma::vec u_l = V_gl * Z_gl;
      znorm = arma::dot(Z_gl, Z_gl);
      residuals += u_l * lambda_g[l];
      //SIMULATE LAMBDA
      if(l==0){temp_var_inv_lambda = znorm + (1.0 / Omega2_g[l]/tauO2 + 1.0 / W2_g[l]);}
      else if(l==(Z_g.size()-1)) {temp_var_inv_lambda = znorm + (1.0/Omega2_g[l-1]/tauO2 + 1.0 / W2_g[l]);}
      else {temp_var_inv_lambda = znorm + (1.0 / Omega2_g[l]/tauO2 + 1.0/Omega2_g[l-1]/tauO2 + 1.0 / W2_g[l] );}
      center_lambda = arma::dot(residuals, u_l);
      if(l==0){center_lambda = (center_lambda+lambda_g[l+1]/Omega2_g[l]/tauO2)/temp_var_inv_lambda;}
      else if(l==(Z_g.size()-1)) {center_lambda = (center_lambda+lambda_g[l-1]/Omega2_g[l-1]/tauO2)/temp_var_inv_lambda;}
      else {center_lambda = (center_lambda+lambda_g[l-1]/Omega2_g[l-1]/tauO2 + (lambda_g[l+1])/Omega2_g[l]/tauO2)/temp_var_inv_lambda;}
      sd_lambda = std::sqrt(sigma2 / temp_var_inv_lambda);
      lambda_g[l] = R::rnorm(center_lambda, sd_lambda);
      //SIMULATE Z
      arma::vec temp_var_inv = ((sigma2/ D_gl) + std::pow(lambda_g[l],2));
      arma::vec center = lambda_g[l]* (V_gl.t() *  residuals)/temp_var_inv;
      arma::vec sd = arma::sqrt(sigma2 / temp_var_inv);
      arma::vec new_Z_gl = (arma::randn(Z_gl.n_elem) % sd) + center;
      Z_g[l] = Rcpp::wrap(new_Z_gl);
      residuals -= V_gl * new_Z_gl * lambda_g[l];
    }
    lambda[g] = lambda_g;
    Z[g] = Z_g;
  }
  return 0;
}

