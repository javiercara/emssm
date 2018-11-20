#'
#' Kalman smoother
#'
#' Kalman smoother for state space model
#'   x_{t+1} = Ax_{t} + w_{t}
#'   y_{t}   = Cx_{t} + v_{t}
#'
# '  cov(w_{t},v_{t}) = [Q 0;0 R]
#'
#' @param y output data
#' @export
#' @return
#' y by rows, ny = nrow(y), nt = ncol(y)
#'
ACQR_ksmoother <- function(A,xtt,Ptt,xtt1,Ptt1){
  #
  # Kalman smoother for model
  #
  # x_{t+1} = Ax_{t} + w_{t}
  # y_{t}   = Cx_{t} + v_{t}
  #
  # cov(w_{t},v_{t}) = [Q 0;0 R]
  #
  # javier.cara@upm.es, 2018-11
  #

  nx = nrow(xtt)
  nt = ncol(xtt)

  # allocation
  xtN = array(0,c(nx,nt+1))
  PtN = array(0,c(nx,nx,nt+1))
  Pt1tN = array(0,c(nx,nx,nt))

  # values for t = nt+1
  xtN[,nt+1] = xtt1[,nt+1]
  PtN[,,nt+1] = Ptt1[,,nt+1]

  # values for t=nt
  xtN[,nt] = xtt[,nt]
  PtN[,,nt] = Ptt[,,nt]
  Pt1tN[,,nt] = A %*% Ptt[,,nt]

  # smother
  for (t in (nt-1):1){
    # Kalman Smoother matrix J
    Jt = ( Ptt[,,t] %*% t(A) ) %*% solve( Ptt1[,,t+1] )

		xtN[,t] = xtt[,t] + Jt %*% ( xtN[,t+1] - xtt1[,t+1] )
		PtN[,,t] = Ptt[,,t] + Jt %*% ( PtN[,,t+1] - Ptt1[,,t+1] ) %*% t(Jt)

    Pt1tN[,,t] = PtN[,,t+1] %*% t(Jt)
  }

	return(list(xtN = xtN, PtN = PtN, Pt1tN = Pt1tN))

}


