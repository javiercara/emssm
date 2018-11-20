#'
#' Kalman filter
#'
#' Kalman filter for state space model
#'
#' @param y: data. Matrix ny*nt
#' A: matrix nx*nx
#' C: matrix ny*nx
#' @export
#' @return
#' y by rows, ny = nrow(y), nt = ncol(y)
#'
ACQR_kfilter <- function(y,A,C,Q,R,x10,P10,nx,ny){
  #
  # Kalman filter for model
  #
  # x_{t+1} = A*x_{t} + w_{t}
  # y_{t}   = C*x_{t} + v_{t}
  #
  # cov(w_{t},v_{t}) = [Q 0;0 R]
  #
  # javier.cara@upm.es, 2018-11
  #

  # data as matrices
  y = as.matrix(y)
  if (nrow(y) != ny){
    y = t(y)
  }
  nt = ncol(y)

  A <- matrix(A, nrow = nx, ncol = nx) # matrix works beter than as.matrix
  C <- matrix(C, nrow=ny, ncol = nx)
  Q <- matrix(Q, nrow = nx, ncol = nx)
  R <- matrix(R, nrow = ny, ncol = ny)
  x10 <- x10
  P10 <- matrix(P10, nrow = nx, ncol = nx)

  # allocation
  xtt <- array(0,c(nx,nt))
  Ptt <- array(0,c(nx,nx,nt))
  xtt1 <- array(0,c(nx,nt+1))
  Ptt1 <- array(0,c(nx,nx,nt+1))
  et <- array(0,c(ny,nt))
  St <- array(0,c(ny,ny,nt))
  Kt <- array(0,c(nx,ny,nt))
  loglik <- 0.0

  # Filter
  xtt1[,1] <- x10
  Ptt1[,,1] <- P10
  for (t in 1:nt){

    #  innovations
    et[,t] = y[,t] - C %*% matrix(xtt1[,t],ncol=1)
    St[,,t] = C %*% Ptt1[,,t] %*% t(C) + R # et variance

    # Kalman gain
    if (ny==1){ Stinv = 1/St[,,t] }
    else{ Stinv = solve(St[,,t]) }
    Kt[,,t] = Ptt1[,,t] %*% t(C) %*% Stinv

    # filtered values
    xtt[,t] = xtt1[,t] + Kt[,,t] %*% matrix(et[,t],ncol=1)
    Ptt[,,t] = (diag(nx) - Kt[,,t] %*% C) %*% Ptt1[,,t]

    # one-step ahead prediction
    xtt1[,t+1] = A %*% xtt[,t]
    Ptt1[,,t+1] = A %*% Ptt[,,t] %*% t(A) + Q

    # likelihood
    if (ny==1){
      Stdet = St[,,t]
    }
    else{
      Stdet = det(St[,,t])
    }
    loglik = loglik + log(Stdet) + t(et[,t]) %*% Stinv %*% et[,t]
  }

  loglik =  - ny*nt/2*log(2*pi) - 0.5*loglik

  return(list(xtt=xtt,Ptt=Ptt,xtt1=xtt1,Ptt1=Ptt1,et=et,St=St,Kt=Kt,loglik=loglik))

}

