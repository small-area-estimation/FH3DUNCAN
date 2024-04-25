
## This R code corresponds to the sections "Real data simulation experiments" and
## ''Prediction of sex occupational segregation by Spanish province'' of the paper 
## Model-based estimation of small area dissimilarity indexes: An application to sex occupational segregation in Spain.
## This script contains the estimation of beta and the prediction of u.

## AUTHORS: Bugallo M., Esteban M.D., Morales D., Pagliarella M. C.,

BETA.U.3foldFH.indep <- function(X, ydrt, D, R, mr, T, sigma2edrt, sigmau1, sigmau2, sigmau3) {
  
  p <- ncol(X)
  
  
  a <- list(1:mr[1]); a
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
  
  
  yd <- Xd <- Vd.inv <- list()
  Q.inv <- matrix(0, nrow=p, ncol=p)
  XVy <- 0
  for(d in 1:D) {
    yd[[d]] <- ydrt[unlist(aa[[d]])]
    Xd[[d]] <- X[unlist(aa[[d]]),]
    
    unoRT <- rep(1,R*T)
    unoT <- rep(1,T)
    
    V1d <- unoRT%*%t(unoRT)
    V2d <- diag(R) %x% (unoT%*%t(unoT))
    V3d <- diag(unoRT)
    Vd <- (sigmau1*V1d + sigmau2*V2d + sigmau3*V3d + diag(sigma2edrt[unlist(aa[[d]])]))
    
    Vd.inv[[d]] <- solve(Vd)                                   
    Q.inv <- Q.inv + t(Xd[[d]])%*%Vd.inv[[d]]%*%Xd[[d]]         
    XVy <- XVy + t(Xd[[d]])%*%Vd.inv[[d]]%*%yd[[d]]            
  }
  Q <- solve(Q.inv)
  
  beta <- Q%*%XVy
  
  u1 <- u2 <- u3 <- list()
  for(d in 1:D){
    
    u1[[d]] <- sigmau1*( t(unoRT)%*%Vd.inv[[d]]%*%(yd[[d]]-Xd[[d]]%*%beta) )
    u2[[d]] <- sigmau2*( diag(R) %x% t(unoT) %*% Vd.inv[[d]]%*%(yd[[d]]-Xd[[d]]%*%beta) )
    u3[[d]] <- sigmau3*( Vd.inv[[d]]%*%(yd[[d]]-Xd[[d]]%*%beta) )
  }
  u1 <- as.matrix(unlist(u1))
  u2 <- as.matrix(unlist(u2))
  u3 <- as.matrix(unlist(u3))
  
  return(list(beta,u1,u2,u3))
}

