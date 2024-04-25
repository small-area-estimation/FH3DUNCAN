

## This R code corresponds to the sections "Real data simulation experiments" and
## ''Prediction of sex occupational segregation by Spanish province'' of the paper 
## Model-based estimation of small area dissimilarity indexes: An application to sex occupational segregation in Spain.
## This script contains the estimation of the variances and the calculation of the p-values

## AUTHORS: Bugallo M., Esteban M.D., Morales D., Pagliarella M. C.,

stderr.3foldFH.indep <- function(fit) {
  Finv <- solve(fit[[2]])
  sigma1.std.err <- sqrt(Finv[1,1])
  sigma2.std.err <- sqrt(Finv[2,2])
  sigma3.std.err <- sqrt(Finv[3,3])
  beta.std.err <- sqrt(as.vector(diag(fit[[5]])))
  
  return( list(beta.std.err, c(sigma1.std.err, sigma2.std.err, sigma3.std.err)) )
  
}

pvalue <- function(beta.fitted, fit) {
  z <- abs(beta.fitted)/sqrt(as.vector(diag(fit[[5]])))
  p <- pnorm(z, lower.tail=F)
  
  return( 2*p )
  
}