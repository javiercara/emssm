#'
#' Build a Hankel matrix
#'
#' Hankel matrix with i rows and j columns from y
#'
#' @param y: data matrix (m x N), i: number of row blocks, j: number of columns
#' @keywords Hankel matrix
#' @export
#' @return
#' H in R^\{m*i x j\}
#' @examples
#'
hankel_yij <- function(y,i,j){

  # y dimensions
  m <- nrow(y)
  N <- ncol(y)

  # Hankel matrix
  H <- matrix(0,m*i,j)
  for (k in 1:i){
    H[((k-1)*m+1):(k*m),] <- y[,k:(k+j-1)]
  }

  # output
  return(H)
}
