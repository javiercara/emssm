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
ABCDQR_em <- function(y,u,Ai,Bi,Ci,Di,Qi,Ri,m1i,P1i,nx,ny,nu,
                    Ae=TRUE,Be=TRUE,Ce=TRUE,De=TRUE,Qe=TRUE,Re=TRUE,m1e=TRUE,P1e=TRUE,
                    max_iter = 100,tol = 1e-6,txo = FALSE){
  #
  # estimate A, B, C, D, Q, R, m1, P1 using the EM algorithm for model
  #
  # x_{t+1} = A*x_{t} + B*u_{t} + w_{t}
  # y_{t}   = C*x_{t} + D*u_{t} v_{t}
  #
  # cov(w_{t},v_{t}) = [Q 0;0 R]
  # x1 -> N(m1,P1)
  #
  # javier.cara@upm.es, 2018-11
  #

  # data as matrices
  y = as.matrix(y)
  if (nrow(y) != ny){
    y = t(y)
  }
  nt = ncol(y)

  u = as.matrix(u)
  if (nrow(u) != nu){
    u = t(u)
  }

  # initial values (as matrices)
  A <- matrix(Ai, nrow = nx, ncol = nx) # matrix works beter than as.matrix
  B <- matrix(Bi, nrow = nx, ncol = nu)
  C <- matrix(Ci, nrow = ny, ncol = nx)
  D <- matrix(Di, nrow = ny, ncol = nu)
  Q <- matrix(Qi, nrow = nx, ncol = nx)
  R <- matrix(Ri, nrow = ny, ncol = ny)
  m1 <- m1i
  P1 <- matrix(P1i, nrow = nx, ncol = nx)

  # log-likelihood values
  loglikv = rep(0,max_iter)

  # Syy does not depend on the iterations
  Syy = array(0,c(ny,ny))
  Suu = array(0,c(nu,nu))
  Syu = array(0,c(ny,nu))
  for (t in 1:nt){
    Syy = Syy + y[,t] %*% t(y[,t])
    Suu = Suu + u[,t] %*% t(u[,t])
    Syu = Syu + y[,t] %*% t(u[,t])
  }

  tol1 = 1.0
  iter = 1
  while ( (iter <= max_iter) && (tol1 > tol) ){
    time1 = proc.time()

    # E-step
    # ---------------------------------------------------------------------------------
    # Kalmanfilter
    kf = ABCDQR_kfilter(y,u,A,B,C,D,Q,R,m1,P1,nx,ny,nu)
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
    Sx1u = array(0,c(nx,nu))
    Sxu = array(0,c(nx,nu))

    # matrices Sxx, Sx1x, Syx, Sx1x1, Sxu, Sx1x1
    for (t in 1:nt){
      Sxx = Sxx + PtN[,,t] + xtN[,t] %*% t(xtN[,t])
      Sx1x = Sx1x + Pt1tN[,,t] + xtN[,t+1] %*% t(xtN[,t])
      Syx = Syx + y[,t] %*% t(xtN[,t])
      Sx1u = Sx1u + xtN[,t+1] %*% t(u[,t])
      Sxu = Sxu + xtN[,t] %*% t(u[,t])
    }

    Sx1x1 = Sxx - (PtN[,,1] + xtN[,1] %*% t(xtN[,1]))  + (PtN[,,nt+1] + xtN[,nt+1] %*% t(xtN[,nt+1]))

    # M-step
    # -------------------------------------------------------------------------------------
    # Matrices m1 y P1
    if (m1e){
      m1 = xtN[,1]
    }
    if (P1e){
      P1 = matrix(PtN[,,1],nrow = nx)
    }

    # Matrix AB
    # AB = [Sx1x Sx1u]/[Sxx Sxu;Sxu' Suu]
    if (Ae || Be){
      AB1 = cbind(Sx1x, Sx1u)
      AB2a = cbind(Sxx, Sxu)
      AB2b = cbind(t(Sxu), Suu)
      AB2 = rbind(AB2a, AB2b)

      AB = AB1 %*% solve(AB2)
    }

    # Matrix A
    if (Ae){
      A = AB[,1:nx]
    }

    # Matrix B
    if (Be){
      B = matrix(AB[,(nx+1):(nx+nu)], nrow = nx)
    }

    # Matrix Q
    if (Qe){
      M1 = Sx1x %*% t(A)
      M2 = Sx1u %*% t(B)
      M3 = A %*% Sxu %*% t(B)
      Q = Sx1x1 - M1 - t(M1) - M2 - t(M2) + M3 + t(M3) + A %*% Sxx %*% t(A) + B %*% Suu %*% t(B)
      Q = 1/nt*Q
      Q = (Q + t(Q))/2 # to make sure it's a symmetric matrix
    }

    # Matrix CD
    # CD = [Syx Syu]/[Sxx Sxu;Sxu' Suu]
    if (Ce || De){
      CD1 = cbind(Syx, Syu)
      CD2a = cbind(Sxx, Sxu)
      CD2b = cbind(t(Sxu), Suu)
      CD2 = rbind(CD2a, CD2b)

      CD = CD1 %*% solve(CD2)
    }

    # Matrix C
    if (Ce){
      C = matrix(CD[,1:nx], nrow = ny)
    }

    # Matrix D
    if (De){
      D = matrix(CD[,(nx+1):(nx+nu)], nrow = ny)
    }

    # Matrix R
    if (Re){
      M1 = Syx %*% t(C)
      M2 = Syu %*% t(D)
      M3 = C %*% Sxu %*% t(D)
      R = Syy - t(M1) - M1 - t(M2) - M2 + t(M3) + M3 + C %*% Sxx %*% t(C) + D %*% Suu %*% t(D)
      R = 1/nt*R
      R = (R + t(R))/2 # to make sure it's a symmetric matrix
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
  P = ny*nx + nx*nx + nx*(nx+1)/2 + ny*(ny+1)/2 + nx*nu + ny*nu
  aic = -2*loglikv[length(loglikv)] + 2*P

  # output
  return( list(A = A, B = B, C = C, D = D, Q = Q, R = R, m1 = m1, P1 = P1,loglikv = loglikv, aic = aic) )
}
