#'
#' forecasting function
#'
#' forecasting function for state space model
#'   x_{t+1} = Ax_{t} + w_{t}
#'   y_{t}   = Cx_{t} + v_{t}
#'   cov(w_{t},v_{t}) = [Q 0;0 R]
#'
#' @param y data
#' @export
#' @return
#' y by rows, m = nrow(y), nt = ncol(y)
#'
ACQR_predict <- function(y,A,C,Q,R,m1,P1,nx,ny,n_ahead,conf_level=0.95){

  # data as matrices
  y = as.matrix(y)
  if (nrow(y) != ny){
    y = t(y)
  }
  nt = ncol(y)

  A <- matrix(A, nrow = nx, ncol = nx)
  #C <-
  Q <- matrix(Q, nrow = nx, ncol = nx)
  R <- matrix(R, nrow = ny, ncol = ny)
  x10 <- m1
  P10 <- matrix(P1, nrow = nx, ncol = nx)

  # predicted values
  yp <- array(0,c(ny,n_ahead))
  # root mean square perdiction errors
  rmspe <- array(0,c(ny,n_ahead))

  # kalman filter
  kf <- ACQR_kfilter(y,A,C,Q,R,x10,P10,nx,ny)

  # forecasting
  xp <- kf$xtt1[,nt+1]
  Pp <- kf$Ptt1[,,nt+1]
  for (t in 1:n_ahead){

    # predictions
    yp[,t] <- C %*% xp

    #  innovations
    St = C %*% Pp %*% t(C) + R # variance
    rmspe[,t] <- sqrt(diag(St))

    # Kalman gain
    Kt = Pp %*% t(C) %*% solve(St)

    # values for the next step
    xp = A %*% xp
    Pp = A %*% ( Pp - Kt %*% C %*% Pp ) %*% t(A) + Q
  }
  # prediction interval
  alpha = 1-conf_level
  ypi1 = yp - qnorm(1-alpha/2)*rmspe
  ypi2 = yp + qnorm(1-alpha/2)*rmspe

  return(list(yp=yp,rmspe=rmspe,ypi1=ypi1,ypi2=ypi2))

}
