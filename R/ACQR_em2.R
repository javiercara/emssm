#'
#' em algorithm
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
ACQR_em2 <- function(y,Ai,Ci,Qi,Ri,m1i,P1i,nx,ny,
                    Ae=matrix(1,nx,nx),Ce=matrix(1,ny,nx),Qe=matrix(1,nx,nx),
                    Re=matrix(1,ny,ny),m1e=matrix(1,nx,1),P1e=matrix(1,nx,nx),
                    max_iter = 100,tol = 1e-6,txo = FALSE){
  #
  # estimate A, C, Q, R, m1, P1 using the EM algorithm for model
  #
  # x_{t+1} = A*x_{t} + w_{t}
  # y_{t}   = C*x_{t} + v_{t}
  #
  # cov(w_{t},v_{t}) = [Q 0;0 R]
  # x1 -> N(m1,P1)
  #
  # javier.cara@upm.es, 2108-11
  #

  # data as matrices
  y = as.matrix(y)
  if (nrow(y) != ny){
    y = t(y)
  }
  nt = ncol(y)

  # initial values (as matrices)
  A <- matrix(Ai, nrow = nx, ncol = nx) # matrix works better than as.matrix
  C <- matrix(Ci, nrow = ny, ncol = nx)
  Q <- matrix(Qi, nrow = nx, ncol = nx)
  R <- matrix(Ri, nrow = ny, ncol = ny)
  m1 <- m1i
  P1 <- matrix(P1i, nrow = nx, ncol = nx)

  # log-likelihood values
  loglikv = rep(0,max_iter)

  # Syy does not depend on the iterations
  Syy = array(0,c(ny,ny))
  for (t in 1:nt){
    Syy = Syy + y[,t] %*% t(y[,t])
  }

  tol1 = 1.0
  iter = 1
  while ( (iter <= max_iter) && (tol1 > tol) ){
    time1 = proc.time()

    # E-step
    # ---------------------------------------------------------------------------------
    # Kalmanfilter
    kf = ACQR_kfilter(y,A,C,Q,R,m1,P1,nx,ny)
    ks = ACQR_ksmoother(A,kf$xtt,kf$Ptt,kf$xtt1,kf$Ptt1)
    xtN = ks$xtN
    PtN = ks$PtN
    Pt1tN = ks$Pt1tN

    loglikv[iter] = kf$loglik
    if (iter > 1){
      tol1 = abs( (loglikv[iter] - loglikv[iter-1])/loglikv[iter-1] )
    }

    # initial values
    Sxx = array(0,c(nx,nx))
    Sx1x = array(0,c(nx,nx))
    Syx = array(0,c(ny,nx))

    # matrices Sxx, Sx1x, Syx, Sx1x1
    for (t in 1:nt){
      Sxx = Sxx + PtN[,,t] + xtN[,t] %*% t(xtN[,t])
      Sx1x = Sx1x + Pt1tN[,,t] + xtN[,t+1] %*% t(xtN[,t])
      Syx = Syx + y[,t] %*% t(xtN[,t])
    }

    Sx1x1 = Sxx - (PtN[,,1] + xtN[,1] %*% t(xtN[,1]))  + (PtN[,,nt+1] + xtN[,nt+1] %*% t(xtN[,nt+1]))

    # M-step
    # -------------------------------------------------------------------------------------
    # Matrices m1 y P1
    m1 = xtN[,1]
    for (i in 1:nx){
      if (m1e[i] == 0){m1[i] = m1i[i]}
    }
    #
    P1 = matrix(PtN[,,1],nrow = nx)
    for (i in 1:nx){
      for (j in 1:nx){
        if (P1e[i,j] == 0){P1[i,j] = P1i[i,j]}
      }
    }

    # Matrix A
    A = Sx1x %*% solve(Sxx)
    for (i in 1:nx){
      for (j in 1:nx){
        if (Ae[i,j] == 0){A[i,j] = Ai[i,j]}
      }
    }

    # Matrix Q
    M1 = Sx1x %*% t(A)
    Q = Sx1x1 - M1 - t(M1) + A %*% Sxx %*% t(A)
    Q = 1/nt*Q
    Q = (Q + t(Q))/2 # to make sure it's a symmetric matrix
    for (i in 1:nx){
      for (j in 1:nx){
        if (Qe[i,j] == 0){Q[i,j] = Qi[i,j]}
      }
    }

    # Matrix C
    C = Syx %*% solve(Sxx)
    for (i in 1:ny){
      for (j in 1:nx){
        if (Ce[i,j] == 0){C[i,j] = Ci[i,j]}
      }
    }

    # Matrix R
    M1 = Syx %*% t(C)
    R = Syy - t(M1) - M1 + C %*% Sxx %*% t(C)
    R = 1/nt*R
    R = (R + t(R))/2 # to make sure it's a symmetric matrix
    for (i in 1:ny){
      for (j in 1:ny){
        if (Re[i,j] == 0){R[i,j] = Ri[i,j]}
      }
    }

    etime0 = proc.time() - time1
    etime = etime0[1] + etime0[2] # user + system
    if (txo){
      cat( "Iter ", sprintf("%3d",iter), ",   @time = ", sprintf("%.2E",etime), ",   logLik = ", sprintf("%.6E",kf$loglik), ",   tol = ", sprintf("%.2E",tol1), "\n", sep = "" )
    }

    iter = iter + 1
  }

  loglikv = loglikv[1:(iter-1)]

  # Akaike Information Criterion
  #P = ny*nx + nx*nx + nx*(nx+1)/2 + ny*(ny+1)/2
  P = sum(Ae) + sum(Qe) + sum(Ce) + sum(Re) # number of estimated parameters
  aic = -2*loglikv[length(loglikv)] + 2*P

  # output
  return( list(A = A, C = C, Q = Q, R = R, m1 = m1, P1 = P1, loglikv = loglikv, aic = aic) )
}
