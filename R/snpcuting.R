#' Cut Chromosome SNP Data into Windows
#'
#' This function divides a list of chromosome matrices into smaller windows of a specified size, which is useful for genomic analysis.
#'
#' @param chr_list A list of matrices, each matrix representing a chromosome's SNP data.
#' @param position_list A list of numeric vectors, where each vector indicates the positions of SNPs within the corresponding chromosome.
#' @param interval Numeric, the size of each window.
#' @param radius Numeric, the radius of each window. If NULL, it defaults to half of the interval size.
#' @param keep_position Logical, whether to return the position list for each window (default is FALSE).
#'
#' @return A list of lists of matrices, where each inner list contains matrices representing the windows for a specific chromosome.
#'
#' @examples
#' chr_list <- lapply(1:5, function(g) { matrix(rbinom(100 * 150, 2, 0.5), 100) })
#' position_list <- lapply(1:5, function(g) { runif(150, 0, 100) })
#' windows <- cut_snp_matrices(chr_list, position_list, interval = 10)
#' @export
cut_snp_matrices = function(chr_list,position_list,interval,radius=NULL,keep_position=FALSE) {
  if(is.null(radius)){radius = interval/2}
  new_position_list = lapply(1:length(chr_list),function(g){
    a = min(position_list[[g]]) + radius
    b = max(position_list[[g]]) - radius
    seq(a,b,interval)
  })
  X_list = lapply(1:length(new_position_list),function(g){
    lapply(new_position_list[[g]],function(center){
      a = center - radius
      b = center + radius
      window = (position_list[[g]]>= a ) & (position_list[[g]]<=b)
      result = as.matrix(chr_list[[g]][,window])
      if (length(result)>0) {return(result)}
      warning("Empty window")
      return(matrix(0,nrow(chr_list[[g]]),1))
    })
  })
  if(keep_position){return(list(X_list=X_list,position_list=new_position_list))}
  return(X_list)
}

#' Calculate Similarity Matrices for SNP Data Windows
#'
#' This function computes the similarity matrices for a list of SNP data windows using specified scaling methods, returning the SVD decomposition for reduced size and increased stability.
#'
#' @param X_list A list of lists of matrices, each representing a window of SNP data.
#' @param svd_inertia Numeric, the proportion of inertia (variance) of the SVD to retain (default is 0.95).
#' @param method Character or function, the method used to scale the data before computing the Gram matrix. Options are "simple", "van_raden", "yang", or a custom function.
#'
#' @return A list of lists, where each inner list contains pairs of matrices (V and D) representing the SVD decomposition of the similarity matrix for each window.
#'
#' @examples
#' chr_list <- lapply(1:5, function(g) { matrix(rbinom(100 * 150, 2, 0.5), 100) })
#' position_list <- lapply(1:5, function(g) { runif(150, 0, 100) })
#' windows <- cut_snp_matrices(chr_list, position_list, interval = 10)
#' similarity_matrices <- get_snp_similarity_matrices(windows, svd_inertia = 0.95, method = "simple")
#' @export
get_snp_similarity_matrices = function(X_list,svd_inertia=0.95,method="simple") {
  if (method %in% c("van_raden","simple","yang")) {
    method = switch(method,{
      "van_raden" = function(X) {
        freq = colMeans(as.matrix(X))/2
        scale(as.matrix(X),center=TRUE,scale=FALSE)/mean(sqrt(freq*(1-freq)))}
      "simple" = function(X) {scale(as.matrix(X))}
      "yang" = function(X) {
        freq = colMeans(as.matrix(X))/2
        scale(as.matrix(X),center=TRUE,scale=sqrt(freq*(1-freq)))}
    })
  }
  else if (!is.function(method)) {
    stop("Method must be a scaling function or one of 'van_raden', 'simple', or 'yang'")
  }
  lapply(X_list,function(X_g){
    lapply(X_g,function(X){
      X = method(X)
      X[which(is.na(X))] = 0
      temp = svd(X)
      if(sum(temp$d^2)!=0){
        temp$d = temp$d^2/sum(temp$d^2)
        nu = min(which(cumsum(temp$d)>=svd_inertia)[1],length(temp$d))
        return(list(V=as.matrix(temp$u[,1:nu]),D=temp$d[1:nu]*nrow(X)/sum(temp$d[1:nu])))
      }
      warning("Similarity_matrix with 0 variance")
      return(list(V=temp$u,D=temp$d))
    })
  })
}

