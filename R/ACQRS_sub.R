#'
#' Estimate the state space model using subspace algorithm
#'
#' Subspace algorithm for the estimation of model
#' \deqn{x_{t+1} = Ax_{t} + w_{t}}
#' \deqn{y_{t} = Cx_{t} + v_{t}}
#' where \eqn{w_{t}} and \eqn{v_{t}} have zero mean and covariance matrices
#' \deqn{Var(w_{t})=Q, Var(v_{t})=R, Cov(w_{t},v_{t})=S}
#'
#' @param \bold{y} data
#' @param \bold{nx} number of rows of matrix A
#' @param \bold{i} auxiliar parameter for building the Hankel matrix. By default, i = nx+1
#' @export
#' @return
#' A, C, Q, R, S
#' @references
#'  \emph{Subspace Identification for Linear Systems.}
#'  \emph{Theory - Implementation - Applications.}
#'  Peter Van Overschee / Bart De Moor.
#'  Kluwer Academic Publishers, 1996
#'
ACQRS_sub <- function(y,nx,ny,i=nx+1){

  # data as matrix and by rows
  y = as.matrix(y)
  if (nrow(y) != ny){
    y = t(y)
  }
  nt = ncol(y)

  # Hankel matrix
  # --------------------------------------------------
  # number of Hankel matrix columns
  j <- nt - 2*i + 1
  if (j < ny*2*i){
    # rows(H) has to be > columns(H)
    stop("Not enough data for building the Hankel matrix")
  }
  H <- hankel_yij(y/sqrt(j),2*i,j)

  # LQ factorization
  # --------------------------------------------------
  q <- qr(t(H))
  L <- t(qr.R(q))
  L21 <- L[(ny*i+1):(2*ny*i),1:(ny*i)]

  # singular values
  # --------------------------------------------------
  s <- svd(L21)
  if (nx==1){
    U1 <- matrix(s$u[,1],ncol=1) # as matrix, not vector
    S1 <- s$d[1]
    }
  else{
    U1 <- s$u[,1:nx]
    S1 <- s$d[1:nx]
    }

  # Matrices gam and gam1
  # --------------------------------------------------
  if (nx==1){
    gam  <- U1 %*% sqrt(S1)
    gam1 <- U1[1:(ny*(i-1)),] %*% sqrt(S1)
  }
  else{
    gam  <- U1 %*% diag(sqrt(S1))
    gam1 <- U1[1:(ny*(i-1)),] %*% diag(sqrt(S1))
  }
  # and pseudo-inverses
  gam_inv  <- pinv(gam)
  gam1_inv <- pinv(gam1)

  # Determine the states Xi and Xi1
  # --------------------------------------------------
  Xi  <- gam_inv  %*% L21
  Xi1 <- gam1_inv %*% L[(ny*(i+1)+1):(2*ny*i),1:(ny*(i+1))]

  # Computing the state matrices A and C
  # --------------------------------------------------
  Rhs <- cbind(Xi, matrix(0,nx,ny))	# Right hand side
  Lhs <- rbind(Xi1, L[(ny*i+1):(ny*(i+1)),1:(ny*(i+1))]) # Left hand side

  # Least squares
  sol <- Lhs %*% pinv(Rhs)

  A <- sol[1:nx,1:nx]
  C <- sol[(nx+1):(nx+ny),1:nx]

  # Computing the covariance matrices Q, R and S
  # -------------------------------------------------
  # Residuals
  res <- Lhs - sol %*% Rhs
  cov <- res %*% t(res)

	Q <- cov[1:nx,1:nx]
	S <- cov[1:nx,(nx+1):(nx+ny)]
	R <- cov[(nx+1):(nx+ny),(nx+1):(nx+ny)]

	return( list(A=A,C=C,Q=Q,R=R,S=S) )

}
