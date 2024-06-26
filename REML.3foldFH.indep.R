
## This R code corresponds to the sections "Real data simulation experiments" and
## ''Prediction of sex occupational segregation by Spanish province'' of the paper 
## Model-based estimation of small area dissimilarity indexes: An application to sex occupational segregation in Spain.
## This script contains the Fisher-Scoring algorithm

## AUTHORS: Bugallo M., Esteban M.D., Morales D., Pagliarella M. C.,


REML.3foldFH.indep <- function(X, ydrt, D, R, mr, T, sigma2edrt, sigma1.0, sigma2.0, sigma3.0, MAXITER) {
  
  sigma1.f <- sigma1.0
  sigma2.f <- sigma2.0
  sigma3.f <- sigma3.0
  theta.f  <- c(sigma1.f, sigma2.f, sigma3.f)
  Bad <- 0
  FLAG <- 0
  p <- ncol(X)
  q <- length(theta.f)
  
  a <- list(1:mr[1])
  mrcum <- cumsum(mr); mrcum
  aa <- list(); aa
  for(d in 1:D){
    if(d>1){
      mrcum <- mrcum+(R*T)
      a[[1]] <- (mrcum[1]-T+1):mrcum[1]
    }
    for(r in 2:R){
      a[[r]] <- (mrcum[r-1]+1):mrcum[r]
    }
    aa[[d]] <- a
  }
  
  yd <- Xd <- list()
  for(d in 1:D) {
    yd[[d]] <- ydrt[unlist(aa[[d]])]
    Xd[[d]] <- X[unlist(aa[[d]]),]
  }
  
  for(ITER in 1:MAXITER){
    Vd.inv <- Vinvyd <- VinvXd <- list()
    
    VinvV1d <- XtVinvV1dVinvX <- VinvV1dVinvV1d <- XtVinvV1dVinvV1dVinvX <- list()
    tr.VinvV1d <- ytVinvV1dVinvy <- ytVinvV1dVinvX <- SumXtVinvV1dVinvX <- tr.VinvV1dVinvV1d <- 0
    
    VinvV2d <- XtVinvV2dVinvX <- VinvV2dVinvV2d <- XtVinvV2dVinvV2dVinvX <- list()
    tr.VinvV2d <- ytVinvV2dVinvy <- ytVinvV2dVinvX <- SumXtVinvV2dVinvX <- tr.VinvV2dVinvV2d <- 0
    
    VinvV3d <- XtVinvV3dVinvX <- VinvV3dVinvV3d <- XtVinvV3dVinvV3dVinvX <- list()
    tr.VinvV3d <- ytVinvV3dVinvy <- ytVinvV3dVinvX <- SumXtVinvV3dVinvX <- tr.VinvV3dVinvV3d <- 0
    
    VinvV1dVinvV2d <- XtVinvV1dVinvV2dVinvX <- list()
    VinvV1dVinvV3d <- XtVinvV1dVinvV3dVinvX <- list()
    VinvV2dVinvV3d <- XtVinvV2dVinvV3dVinvX <- list()
    
    ytVinvX <- tr.VinvV1dVinvV2d <- tr.VinvV1dVinvV3d <- tr.VinvV2dVinvV3d <-0
    
    Q.inv <- matrix(0, nrow=p, ncol=p)
    for(d in 1:D) {
      unoRT <- rep(1,R*T)
      unoT <- rep(1,T)
      V1d <- unoRT%*%t(unoRT)
      V2d <- diag(R) %x% (unoT%*%t(unoT))
      V3d <- diag(unoRT)
      Vd <- (sigma1.f*V1d + sigma2.f*V2d + sigma3.f*V3d + diag(sigma2edrt[unlist(aa[[d]])]))
      
      if(abs(det(Vd))<0.000000001 || abs(det(Vd))>100000000000) {
      #  FLAG <- 1
        Bad <- Bad+1
        #cat("\n Check occurred at domain d =", d, "due to det(Vd)=",det(Vd),"\n")
        #cat("\n ITER = ", ITER, ", Bad = ", Bad, "and FLAG = ", FLAG, "\n")
        #cat("\n thetas = ", sigma1.f, sigma2.f, sigma3.f, "\n")
      #  break
      }
      
      Vd.inv[[d]] <- solve(Vd)
      Vinvyd[[d]] <- Vd.inv[[d]]%*%yd[[d]]
      VinvXd[[d]] <- Vd.inv[[d]]%*%Xd[[d]]
      Q.inv <- Q.inv + t(Xd[[d]])%*%VinvXd[[d]]
      
      VinvV1d[[d]] <- Vd.inv[[d]]%*%V1d
      tr.VinvV1d <- tr.VinvV1d + sum(diag(VinvV1d[[d]]))
      XtVinvV1dVinvX[[d]] <- t(VinvXd[[d]])%*%V1d%*%VinvXd[[d]]
      ytVinvX <- ytVinvX + t(yd[[d]])%*%VinvXd[[d]]
      ytVinvV1dVinvy <- ytVinvV1dVinvy + t(Vinvyd[[d]])%*%V1d%*%Vinvyd[[d]]
      ytVinvV1dVinvX <- ytVinvV1dVinvX + t(Vinvyd[[d]])%*%V1d%*%VinvXd[[d]]
      SumXtVinvV1dVinvX <- SumXtVinvV1dVinvX + XtVinvV1dVinvX[[d]]
      
      VinvV2d[[d]] <- Vd.inv[[d]]%*%V2d
      tr.VinvV2d <-  tr.VinvV2d + sum(diag(VinvV2d[[d]]))
      XtVinvV2dVinvX[[d]] <- t(VinvXd[[d]])%*%V2d%*%VinvXd[[d]]
      ytVinvV2dVinvy <- ytVinvV2dVinvy + t(Vinvyd[[d]])%*%V2d%*%Vinvyd[[d]]
      ytVinvV2dVinvX <- ytVinvV2dVinvX + t(Vinvyd[[d]])%*%V2d%*%VinvXd[[d]]
      SumXtVinvV2dVinvX <- SumXtVinvV2dVinvX + XtVinvV2dVinvX[[d]]
      
      VinvV3d[[d]] <- Vd.inv[[d]]%*%V3d
      tr.VinvV3d <-  tr.VinvV3d + sum(diag(VinvV3d[[d]]))
      XtVinvV3dVinvX[[d]] <- t(VinvXd[[d]])%*%V3d%*%VinvXd[[d]]
      ytVinvV3dVinvy <- ytVinvV3dVinvy + t(Vinvyd[[d]])%*%V3d%*%Vinvyd[[d]]
      ytVinvV3dVinvX <- ytVinvV3dVinvX + t(Vinvyd[[d]])%*%V3d%*%VinvXd[[d]]
      SumXtVinvV3dVinvX <- SumXtVinvV3dVinvX + XtVinvV3dVinvX[[d]]
      
      VinvV1dVinvV1d[[d]] <- VinvV1d[[d]]%*%VinvV1d[[d]]
      tr.VinvV1dVinvV1d <- tr.VinvV1dVinvV1d + sum(diag(VinvV1dVinvV1d[[d]]))
      XtVinvV1dVinvV1dVinvX[[d]]  <- t(Xd[[d]])%*%VinvV1dVinvV1d[[d]]%*%VinvXd[[d]]
      
      VinvV1dVinvV2d[[d]] <- VinvV1d[[d]]%*%VinvV2d[[d]]
      tr.VinvV1dVinvV2d <- tr.VinvV1dVinvV2d + sum(diag(VinvV1dVinvV2d[[d]]))
      XtVinvV1dVinvV2dVinvX[[d]] <- t(Xd[[d]])%*%VinvV1dVinvV2d[[d]]%*%VinvXd[[d]]
      
      VinvV1dVinvV3d[[d]] <- VinvV1d[[d]]%*%VinvV3d[[d]]
      tr.VinvV1dVinvV3d <- tr.VinvV1dVinvV3d + sum(diag(VinvV1dVinvV3d[[d]]))
      XtVinvV1dVinvV3dVinvX[[d]] <- t(Xd[[d]])%*%VinvV1dVinvV3d[[d]]%*%VinvXd[[d]]
      
      VinvV2dVinvV2d[[d]] <- VinvV2d[[d]]%*%VinvV2d[[d]]
      tr.VinvV2dVinvV2d <- tr.VinvV2dVinvV2d + sum(diag(VinvV2dVinvV2d[[d]]))
      XtVinvV2dVinvV2dVinvX[[d]] <- t(Xd[[d]])%*%VinvV2dVinvV2d[[d]]%*%VinvXd[[d]]
      
      VinvV2dVinvV3d[[d]] <- VinvV2d[[d]]%*%VinvV3d[[d]]
      tr.VinvV2dVinvV3d <- tr.VinvV2dVinvV3d + sum(diag(VinvV2dVinvV3d[[d]]))
      XtVinvV2dVinvV3dVinvX[[d]] <- t(Xd[[d]])%*%VinvV2dVinvV3d[[d]]%*%VinvXd[[d]]
      
      VinvV3dVinvV3d[[d]] <- VinvV3d[[d]]%*%VinvV3d[[d]]
      tr.VinvV3dVinvV3d <- tr.VinvV3dVinvV3d + sum(diag(VinvV3dVinvV3d[[d]]))
      XtVinvV3dVinvV3dVinvX[[d]] <- t(Xd[[d]])%*%VinvV3dVinvV3d[[d]]%*%VinvXd[[d]]
      
      if(FLAG==1){
        FLAG <- 0
        cat("\n Flag = ", FLAG, "\n")
        cat("\n ITER = ", ITER, "\n")
        ITER <- MAXITER
        cat("\n ITER new = ", ITER, "\n")
        break
      }
    }
    Q <- solve(Q.inv)
    #Q <- ginv(Q.inv)
    
    tr.XtVinvV1dVinvXQ <- tr.XtVinvV2dVinvXQ <- tr.XtVinvV3dVinvXQ <- 0
    tr.XtVinvV1dVinvV1dVinvXQ <- tr.XtVinvV1dVinvV2dVinvXQ <- tr.XtVinvV1dVinvV3dVinvXQ <- 0
    tr.XtVinvV2dVinvV2dVinvXQ <- tr.XtVinvV2dVinvV3dVinvXQ <- 0
    tr.XtVinvV3dVinvV3dVinvXQ <- 0
    for(d in 1:D){
    
      tr.XtVinvV1dVinvXQ <- tr.XtVinvV1dVinvXQ + sum(diag(XtVinvV1dVinvX[[d]]%*%Q))
      tr.XtVinvV2dVinvXQ <- tr.XtVinvV2dVinvXQ + sum(diag(XtVinvV2dVinvX[[d]]%*%Q))
      tr.XtVinvV3dVinvXQ <- tr.XtVinvV3dVinvXQ + sum(diag(XtVinvV3dVinvX[[d]]%*%Q))
      
      tr.XtVinvV1dVinvV1dVinvXQ <- tr.XtVinvV1dVinvV1dVinvXQ + sum(diag(XtVinvV1dVinvV1dVinvX[[d]]%*%Q))
      tr.XtVinvV1dVinvV2dVinvXQ <- tr.XtVinvV1dVinvV2dVinvXQ + sum(diag(XtVinvV1dVinvV2dVinvX[[d]]%*%Q))
      tr.XtVinvV1dVinvV3dVinvXQ <- tr.XtVinvV1dVinvV3dVinvXQ + sum(diag(XtVinvV1dVinvV3dVinvX[[d]]%*%Q))
      
      tr.XtVinvV2dVinvV2dVinvXQ <- tr.XtVinvV2dVinvV2dVinvXQ + sum(diag(XtVinvV2dVinvV2dVinvX[[d]]%*%Q))
      tr.XtVinvV2dVinvV3dVinvXQ <- tr.XtVinvV2dVinvV3dVinvXQ + sum(diag(XtVinvV2dVinvV3dVinvX[[d]]%*%Q))
      
      tr.XtVinvV3dVinvV3dVinvXQ <- tr.XtVinvV3dVinvV3dVinvXQ + sum(diag(XtVinvV3dVinvV3dVinvX[[d]]%*%Q))
    }
    
    tr.PV1 <- tr.VinvV1d - tr.XtVinvV1dVinvXQ 
    tr.PV2 <- tr.VinvV2d - tr.XtVinvV2dVinvXQ
    tr.PV3 <- tr.VinvV3d - tr.XtVinvV3dVinvXQ
    
    SumXtVinvV1dVinvXQ <- SumXtVinvV1dVinvX%*%Q
    SumXtVinvV2dVinvXQ <- SumXtVinvV2dVinvX%*%Q
    SumXtVinvV3dVinvXQ <- SumXtVinvV3dVinvX%*%Q
    
    
    tr.PV1PV1  <- tr.VinvV1dVinvV1d - 2*tr.XtVinvV1dVinvV1dVinvXQ + sum(diag(SumXtVinvV1dVinvXQ%*%SumXtVinvV1dVinvXQ))
    tr.PV1PV2  <- tr.VinvV1dVinvV2d - 2*tr.XtVinvV1dVinvV2dVinvXQ + sum(diag(SumXtVinvV1dVinvXQ%*%SumXtVinvV2dVinvXQ))
    tr.PV1PV3  <- tr.VinvV1dVinvV3d - 2*tr.XtVinvV1dVinvV3dVinvXQ + sum(diag(SumXtVinvV1dVinvXQ%*%SumXtVinvV3dVinvXQ))
    
    tr.PV2PV2  <- tr.VinvV2dVinvV2d - 2*tr.XtVinvV2dVinvV2dVinvXQ + sum(diag(SumXtVinvV2dVinvXQ%*%SumXtVinvV2dVinvXQ))
    tr.PV2PV3  <- tr.VinvV2dVinvV3d - 2*tr.XtVinvV2dVinvV3dVinvXQ + sum(diag(SumXtVinvV2dVinvXQ%*%SumXtVinvV3dVinvXQ))
    
    tr.PV3PV3  <- tr.VinvV3dVinvV3d - 2*tr.XtVinvV3dVinvV3dVinvXQ + sum(diag(SumXtVinvV3dVinvXQ%*%SumXtVinvV3dVinvXQ))
    
    ytVinvXQ <- ytVinvX%*%Q
    
    ytPV1Py  <- ytVinvV1dVinvy - 2*ytVinvV1dVinvX%*%t(ytVinvXQ) + ytVinvXQ%*%SumXtVinvV1dVinvX%*%t(ytVinvXQ)
    ytPV2Py  <- ytVinvV2dVinvy - 2*ytVinvV2dVinvX%*%t(ytVinvXQ) + ytVinvXQ%*%SumXtVinvV2dVinvX%*%t(ytVinvXQ)
    ytPV3Py  <- ytVinvV3dVinvy - 2*ytVinvV3dVinvX%*%t(ytVinvXQ) + ytVinvXQ%*%SumXtVinvV3dVinvX%*%t(ytVinvXQ)
    
    S1 <- -0.5*tr.PV1 + 0.5*ytPV1Py
    S2 <- -0.5*tr.PV2 + 0.5*ytPV2Py
    S3 <- -0.5*tr.PV3 + 0.5*ytPV3Py
    
    F11 <- 0.5*tr.PV1PV1
    F12 <- 0.5*tr.PV1PV2
    F13 <- 0.5*tr.PV1PV3
    F22 <- 0.5*tr.PV2PV2
    F23 <- 0.5*tr.PV2PV3
    F33 <- 0.5*tr.PV3PV3
    
    Ssig <- c(S1,S2,S3)
    Fsig <- matrix(c(F11,F12,F13,F11,F22,F23,F13,F23,F33), ncol=q)
    
    Fsig.inv <- solve(Fsig)
    #Fsig.inv <- ginv(Fsig)
    dif <- Fsig.inv%*%Ssig
    theta.f <- theta.f + dif
    sigma1.f <- theta.f[1,1]
    sigma2.f <- theta.f[2,1]
    sigma3.f <- theta.f[3,1]
    
    if(identical(as.numeric(abs(dif)<0.000001),rep(1,q))){
      break
    }
    
    if(sigma1.f<0 || sigma2.f<0 || sigma3.f<0) {
      #cat("\n Bad occurred at check for sigma",which(theta.f<0)," < 0 at iteration ",ITER,"\n")
      Bad <- Bad+1
      
      theta.f[which(theta.f<0)] <- 0.00001
      sigma1.f <- theta.f[1,1]
      sigma2.f <- theta.f[2,1]
      sigma3.f <- theta.f[3,1]
    }
    
    
  }
  
  #if(sigma1.f<0 || sigma2.f<0 || sigma3.f<0) {
  #  ITER <- MAXITER
  #  Bad <- Bad+1
  #}
  
  
  return(list(as.vector(theta.f), Fsig, ITER, Bad, Q))
  
}



