#' covboost
#'
#' \code{covboost} estimates covariance matrices through boosting.
#' This allows automatic sparsity detection. The biggest benefits are perhaps
#' for high dimensional \code{p>>n} problems. Imports RcppArmadillo and OpenMP.
#'
#' Important functions:
#'
#' \itemize{
#' \item \code{\link{covboost}}: function for boosting out a covariance matrix from the Identity matrix
#' \item \code{\link{covboost_cv}}: The same as \code{\link{covboost}}, but adds k-fold cross validation
#' }
#'
#' See individual function documentation for usage.
#'
#' @docType package
#' @author Berent Ånund Strømnes Lunde
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib covboost
#' @name covboost-package
#' @title Estimate Covariance Matrices Through Boosting
NULL
