#' Generate Test Data for Simulation Studies
#'
#' This function generates simulated test data for statistical modeling, particularly useful
#' in contexts where grouped and multilevel data structures are involved. It constructs
#' a hierarchical structure of scaled matrix data (X_list), singular vectors (V_list),
#' normalized singular values (D_list), random effects coefficients (lambda), and
#' responses (y) based on the parameters specified.
#'
#' @param n Integer, the number of observations in each generated matrix. Default is 400
#' @param k Integer, the number of groups of matrices to generate. Each group will have `p` matrices. Default is 3.
#' @param p Integer, the number of matrices (features) in each group. Default is 15.
#' @param nu Integer, the number of columns in each matrix within a group, effectively controlling the dimensionality. Default is 5.
#' @param mu Numeric, the global mean used to generate the response vector `y`. Default is 0.
#' @param signal_percentage Numeric, the percentage of signal in the data. Default is 0.5.
#' @param n_signal Integer, the number of non-zero coefficients in each group. Default is 3.
#' @param contiguous Logical, whether to generate contiguous signals in the matrices. Default is FALSE.
#'
#' @return A list containing several components:
#' \itemize{
#'   \item{y}{Numeric vector, the response variable generated from the simulated data.}
#'   \item{V_list}{List of lists, where each sublist contains the left singular vectors of the matrices from X_list.}
#'   \item{D_list}{List of lists, where each sublist contains the normalized squared singular values of the matrices from X_list.}
#'   \item{sigma2}{Numeric, the variance used to generate the response vector `y`.}
#'   \item{mu}{Numeric, the global mean used in the generation of `y`.}
#'   \item{lambda}{List of numeric vectors, the random effects coefficients used for each group.}
#'   \item{X_list}{List of lists, where each sublist contains the scaled matrices generated from the data.}
#'   \item{Z_list}{List of lists, where each sublist contains the random effects used to generate the response vector `y`.}
#' }
#'
#' @examples
#' test_data = generate_test_data()
#' names(test_data)
#' @export
generate_test_data = function(n=400,k=3,p=15,nu=5,mu=0,signal_percentage=0.5,n_signal=3,contiguous=FALSE) {
  contiguous_start = sample(n_signal:(p-n_signal+1),1)
  contiguous_zone = (contiguous_start+1):(contiguous_start+n_signal-1)
  X_list = lapply(1:k,function(g){
    if(!contiguous){
      X_g = lapply(1:p,function(l){
        scale(matrix(rnorm(n*nu),n,nu))
      })
    }
    else {
      X_g = lapply(1:p,function(l){
        scale(matrix(rnorm(n*nu),n,nu))
      })
      for(i in contiguous_zone)
      X_g[[i]] = X_g[[contiguous_start]] + scale(matrix(rnorm(n*nu),n,nu))/3
    }
    return(X_g)
  })
  V_list = lapply(X_list,function(X_g){
    lapply(X_g,function(X_l){
      svd(X_l)$u
    })
  })
  D_list = lapply(X_list,function(X_g){
    lapply(X_g,function(X_l){
      temp_svd = svd(X_l)
      temp_svd$d^2*nrow(X_l)/sum(temp_svd$d^2)
    })
  })
  if(!contiguous) {
    lambda = lapply(1:k,function(g){
      result = rep(0,p)
      signal_index = round((p-1)/(n_signal-1+1)*(1:(n_signal)-1) +1 )
      result[signal_index] = sqrt(signal_percentage/(n_signal*k))
      return(result)
    })
  }
  else {
    lambda = lapply(1:k,function(g){
      result = rep(0,p)
      result[c(contiguous_start,contiguous_zone)] = sqrt(signal_percentage/(n_signal*k))
      return(result)
    })
  }
  Z = lapply(D_list,function(D_g){
    lapply(D_g,function(D_l){
      temp = sqrt(D_l) * rnorm(length(D_l))
      temp/sqrt(sum(temp^2)/n)
    })
  })
  y = Reduce("+",lapply(1:k,function(g){
    Reduce("+",lapply(1:p,function(l){
      lambda[[g]][l]*V_list[[g]][[l]]%*%(Z[[g]][[l]])
    }))
  })) + mu + (rnorm(n)*sqrt(1-signal_percentage))
  return(list(y=y,V_list=V_list,D_list=D_list,signal_percentage=signal_percentage,mu=mu,lambda=lambda,X_list=X_list,Z_list=Z))
}

#' Generate SNP Data for Simulation Studies
#'
#' This function generates simulated SNP data for statistical modeling, particularly useful
#' in contexts where genetic data are involved. It constructs a hierarchical structure of
#' chromosome data (chr_list) and position data (position_list) based on the parameters specified.
#'
#' @param n_chr Integer, the number of groups of chromosomes to generate. Default is 1.
#' @param n_rows Integer, the number of rows in each chromosome matrix. Default is 100.
#' @param n_cols Integer, the number of columns in each chromosome matrix. Default is 100.
#' @param by_chr Logical, whether to generate data by chromosome (TRUE) or by SNP (FALSE). Default is TRUE.
#'
#' @return A list containing two components:
#' \itemize{
#'   \item{chr_list}{List of matrices, where each matrix represents a chromosome.}
#'   \item{position_list}{List of numeric vectors, where each vector represents the position of each column in the corresponding chromosome matrix.}
#' }
#'
#' @examples
#' snp_data = generate_snp_data()
#' names(snp_data)
#' @export
generate_snp_data = function(n_chr=3,n_rows=100,n_cols=100,by_chr=TRUE) {
  if(by_chr){
    chr_list = lapply(1:n_chr,function(chr){
      matrix(rbinom(n_rows * n_cols, 2, 0.5), n_rows)
    })
    position_list = lapply(1:n_chr,function(chr){
      runif(n_cols, 0, n_cols)
    })
    lambda_list = lapply(1:n_chr,function(chr){
      rnorm(n_cols)
    })
    y = Reduce("+",lapply(1:length(chr_list),function(chr){
      chr_list[[chr]]%*%(lambda_list[[chr]]*rnorm(n_cols))
    }))
    return(list(chr_list=chr_list,position_list=position_list,y=y))
  }
  else {
    X = matrix(rbinom(n_chr*n_rows*n_cols,2,0.5),n_rows)
    chromosome_list = unlist(lapply(1:n_chr,function(chr){
      rep(chr,n_cols)
    }))
    position_list = unlist(lapply(1:n_chr,function(chr){
      runif(n_cols,0,n_cols)
    }))
    y = X%*%(rnorm(n_cols*n_chr)*rnorm(n_cols*n_chr))
    return(list(y=y,X=X,chromosome_list=chromosome_list,position_list=position_list))
  }
}

