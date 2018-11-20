#'
#' Matrix pseudo-inverse
#'
#' Matrix pseudo-inverse
#'
#' @param X matrix
#' @keywords pseudoinverse
#' @export
#' @return
#' X
#'
pinv <- function(X, tol = sqrt(.Machine$double.eps)){
  # Generalized Inverse of a Matrix
  # Venables, W. N. and Ripley, B. D. (1999)
  # Modern Applied Statistics with S-PLUS. Third Edition. Springer. p.100.
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
    dimnames = dnx[2:1])
}
