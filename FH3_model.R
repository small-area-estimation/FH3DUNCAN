
## This R code corresponds to the section ''Prediction of sex occupational segregation by Spanish province'' of the paper 
## Model-based estimation of small area dissimilarity indexes: An application to sex occupational segregation in Spain.
## This script contains the main code of the case study

## AUTHORS: Bugallo M., Esteban M.D., Morales D., Pagliarella M. C.,

rm(list=ls())

source("REML.3foldFH.indep.R")
source("BETA.U.3foldFH.indep.R")
source("stderr.3foldFH.indep.R")

library(RColorBrewer)

### Three-fold Fay-Herriot model and Duncan Segregation Index

########################################
######## PART 1 
########################################

# Real data on man proportion (Y) and auxiliary variables got from EPA_2020_T4 and EPA_2021_T1 -- EPA_2021_T4
data = read.csv("YHajek_covariates.csv", sep=',', header=T) 

data$Province = factor(data$Province)
data$Ocupation = factor(data$Ocupation)
data$Period = factor(data$Period)

# Small Area Estimation Problem
summary(data$nd)
par(mar=c(4, 4.3, 4, 4), xpd=T)
hist(data$nd, cex.axis=1.8, cex.main=1.85, cex.lab=1.85, main='', xlab='Sample size', col='gray80',
     freq=F, ylab='Probability mass')
quantile(data$nd, probs=seq(0,1, by=0.1))	

summary(100*data$nd/data$hatNd)
par(mar=c(4, 4.3, 4, 4), xpd=T)
hist(100*data$nd/data$hatNd, cex.axis=1.8, cex.main=1.85, cex.lab=1.85, main='', xlab='Sampling fractions (%)', col='gray80',
     freq=F, ylab='Probability mass')
round(quantile(100*data$nd/data$hatNd, probs=seq(0,1, by=0.1)), 4)	

# Sampling variances
sigma2edrt <- data$Yvar.mean
summary(sigma2edrt)

# Auxiliary variables: delete 'edu3', 'age41', 'age43 and 'situ3'
X <- as.matrix(data[,c(6, 9:10, 12:15, 17)]) 
p <- ncol(X)

# Response variable Y (direct estimates)
ydrt <- data$Ymean

df_1<-data.frame(aggregate(ydrt, by=list(data$Ocupation, data$Period), sum)) 
df_1$Prop=df_1$x/52
df_1$Tot=data.frame(aggregate(data$nd, by=list(data$Ocupation, data$Period), sum))$x 

# Number of domains D
D <- length(unique(data[,"Province"])); D

# Number of subdomain R
R <- length(unique(data[,"Ocupation"])); R

# Number of periods R
T <- length(unique(data[,"Period"])); T

# Vector of length D, containing repeated values of R
md <- rep(R, D)

# Vector of length R, containing repeated values of T
mr <- rep(T, R)


########################################
######## PART 2 
########################################

# Fitting REML algorithm
# Setting a starting value 
sigma.0 = 0.5
fitReml <- REML.3foldFH.indep(X, ydrt, D, R, mr, T, sigma2edrt, sigma1.0 = sigma.0, 
                              sigma2.0 = sigma.0, sigma3.0 = sigma.0, MAXITER = 40)

# Estimated variance parameters
sigmau1.hat <- fitReml[[1]][1]
sigmau2.hat <- fitReml[[1]][2]
sigmau3.hat <- fitReml[[1]][3]
Fsig <- fitReml[[2]] # Fisher-scoring matrix

# Regression parameters beta
fitBeta <- BETA.U.3foldFH.indep(X, ydrt, D, R, mr, T, sigma2edrt, sigmau1.hat, sigmau2.hat, sigmau3.hat)

# Estimated betas coefficients
beta <- as.vector(fitBeta[[1]]); beta

# Random effects
u1 <- as.vector(fitBeta[[2]])
u2 <- as.vector(fitBeta[[3]])
u3 <- as.vector(fitBeta[[4]])

# Standard error for betas hat and sigmas hat
stderr <- stderr.3foldFH.indep(fitReml)

# Confidence intervals at 95% for betas
conf.level <- 0.95
alpha <- 1-conf.level
k <- 1-alpha/2
z <- qnorm(k)

lower.b <- beta-z*stderr[[1]]
upper.b <- beta+z*stderr[[1]]

# p-value for betas
pvalue.b <- pvalue(beta, fitReml)

# Summary of betas estimates
summary_betas <- cbind(beta.hat=beta, std.error=stderr[[1]], t=beta/stderr[[1]], 
                       p.value=round(pvalue.b,6), lower=lower.b, upper=upper.b)
row.names(summary_betas) <- colnames(X); summary_betas


# Confidence intervals at 95% for sigmas
lower.s <- c(sigmau1.hat-z*stderr[[2]][1], sigmau2.hat-z*stderr[[2]][2], sigmau3.hat-z*stderr[[2]][3])
upper.s <- c(sigmau1.hat+z*stderr[[2]][1], sigmau2.hat+z*stderr[[2]][2], sigmau3.hat+z*stderr[[2]][3])
IC.s <- cbind(lower.s, upper.s)
row.names(IC.s) <- c('sigma1', 'sigma2', 'sigma3'); round(IC.s, 4)

test.s <- c(lower.s[1]<0 & upper.s[1]>0, lower.s[2]<0 & upper.s[2]>0, lower.s[3]<0 & upper.s[3]>0); test.s

# Confidence intervals at 95% for the differences of sigmas
lower.s.dif <- c((sigmau1.hat - sigmau2.hat) - z*sqrt(solve(Fsig)[1,1] + solve(Fsig)[2,2] -2*solve(Fsig)[1,2]),
                 (sigmau1.hat - sigmau3.hat) - z*sqrt(solve(Fsig)[1,1] + solve(Fsig)[3,3] -2*solve(Fsig)[1,3]),
                 (sigmau2.hat - sigmau3.hat) - z*sqrt(solve(Fsig)[2,2] + solve(Fsig)[3,3] -2*solve(Fsig)[2,3]))

upper.s.dif <- c((sigmau1.hat - sigmau2.hat) + z*sqrt(solve(Fsig)[1,1] + solve(Fsig)[2,2] -2*solve(Fsig)[1,2]),
                 (sigmau1.hat - sigmau3.hat) + z*sqrt(solve(Fsig)[1,1] + solve(Fsig)[3,3] -2*solve(Fsig)[1,3]),
                 (sigmau2.hat - sigmau3.hat) + z*sqrt(solve(Fsig)[2,2] + solve(Fsig)[3,3] -2*solve(Fsig)[2,3]))
test.s.dif <- c(lower.s.dif[1]<0 & upper.s.dif[1]>0, lower.s.dif[2]<0 & upper.s.dif[2]>0, lower.s.dif[3]<0 & upper.s.dif[3]>0); test.s.dif 

# Summary of sigma estimates
summary_sigmas <- cbind(sigmas.hat=fitReml[[1]], std.error=stderr.3foldFH.indep(fitReml)[[2]], 
                        lower=lower.s, upper=upper.s, test=as.character(test.s))
row.names(summary_sigmas) <- c("sigmau1", "sigmau2", "sigmau3"); summary_sigmas

summary_sigmas.diff <- cbind(lower=lower.s.dif, upper=upper.s.dif, test=as.character(test.s.dif))
row.names(summary_sigmas) <- c("sigmau1-sigmau2", "sigmau1-sigmau3", "sigmau2-sigmau3"); summary_sigmas.diff


########################################
######## PART 3 
########################################

### An index of dissimilarity: The Duncan Segregation Index (DSI)
## Using the Hajek direct estimator for domain means

Ydrt1_ = data$Ymean
Ydrt2_ = 1 - data$Ymean
Ndrt = data$hatNd

Nd.t = data.frame(aggregate(Ndrt, by=list(data$Province, data$Period), sum))$x

Yd.t1_ = data.frame(aggregate(Ydrt1_*Ndrt, by=list(data$Province, data$Period), sum))$x/Nd.t
Yd.t2_ = data.frame(aggregate(Ydrt2_*Ndrt, by=list(data$Province, data$Period), sum))$x/Nd.t

Sdrt_dir = Ndrt /Nd.t * abs (Ydrt1_/Yd.t1_ - Ydrt2_/Yd.t2_)
Sd.t_prov_period_dir = data.frame(aggregate(Sdrt_dir, by=list(data$Province, data$Period), sum))
Sd.t_prov_period_dir$x = 1/2 *Sd.t_prov_period_dir$x
names(Sd.t_prov_period_dir)=c('Province', 'Period', 'Sd.t_dir')

# Hajek DSI for each province d and period time t
Sd.t_dir = Sd.t_prov_period_dir$Sd.t_dir

# For T4 2021
Sd.t_2021T4_dir= Sd.t_dir[209:260]


########################################
######## PART 4 
########################################
### An index of dissimilarity: The Duncan Segregation Index (DSI)

mudrt1 = X%*%beta + rep(u1, 7*5) + rep(u2, 5) + u3
mudrt2 = 1 - mudrt1 

mud.t1 = data.frame(aggregate(mudrt1*Ndrt, by=list(data$Province, data$Period), sum))$V1
mud.t2 = data.frame(aggregate(mudrt2*Ndrt, by=list(data$Province, data$Period), sum))$V1

Sdrt = abs(mudrt1*Ndrt/mud.t1 - mudrt2*Ndrt/mud.t2)

Sd.t_prov_period = data.frame(aggregate(Sdrt, by=list(data$Province, data$Period), sum))
Sd.t_prov_period$V1 = 1/2 *Sd.t_prov_period$V1
names(Sd.t_prov_period)=c('Province', 'Period', 'Sd.t')

# Model DSI for each province d and period time t
Sd.t = Sd.t_prov_period$Sd.t

# For T4 2021
Sd.t_2021T4 = Sd.t[209:260]

# For T4 2020
Sd.t_2020T4 = Sd.t[1:52]

write.table(data.frame('PROV'=1:52, 'Sd.t_2021T4'=Sd.t_2021T4), "PI_Duncan.txt", 
            sep=";",col.names = TRUE, row.names = FALSE)

write.table(data.frame('PROV'=1:52, 'Sd.t_dif'=Sd.t_2021T4 - Sd.t_2020T4), "PI_DuncanDIFER.txt", 
            sep=";",col.names = TRUE, row.names = FALSE)	

##### OCUP1

# Model
Sdrt_OCUP1.1 = Sdrt[data$Ocupation==1 & data$Period==5]
Sdrt_OCUP1.2 = Sdrt[data$Ocupation==2 & data$Period==5]
Sdrt_OCUP1.3 = Sdrt[data$Ocupation==3 & data$Period==5]
Sdrt_OCUP1.4 = Sdrt[data$Ocupation==4 & data$Period==5]
Sdrt_OCUP1.5 = Sdrt[data$Ocupation==5 & data$Period==5]
Sdrt_OCUP1.6 = Sdrt[data$Ocupation==6 & data$Period==5]
Sdrt_OCUP1.7 = Sdrt[data$Ocupation==7 & data$Period==5]

write.table(data.frame('PROV'=1:52, 'Sdrt_OCUP1.2'=Sdrt_OCUP1.2), "PI_DuncanR2.txt", 
            sep=";",col.names = TRUE, row.names = FALSE)

Sdrt_OCUP1.2DIFER = Sdrt[data$Ocupation==2 & data$Period==5] - Sdrt[data$Ocupation==2 & data$Period==1]			

write.table(data.frame('PROV'=1:52, 'Sdrt_OCUP1.2'=Sdrt_OCUP1.2DIFER), "PI_DuncanR2DIFER.txt", 
            sep=";",col.names = TRUE, row.names = FALSE)

# Main occupation
sd = 2* Sd.t[209:260]/R

write.table(data.frame('PROV'=1:52, 'sd'=sd), "PI_Duncan_sd.txt", 
            sep=";",col.names = TRUE, row.names = FALSE)

par(mfrow=c(2,4), mar=c(4, 4, 4, 4), xpd=T)
plot(Sdrt_OCUP1.1, type='l', ylim=c(0,0.41), lwd=2, main='OC1', ylab='', xlab='Province', 
     cex.lab=1.6, cex.axis=1.5, cex.main=1.6)
lines(sd,  col=brewer.pal(10,"RdYlBu")[9], lwd=2, lty=3)

plot(Sdrt_OCUP1.2, type='l', ylim=c(0,0.41), lwd=2, main='OC2', ylab='', xlab='Province', 
     cex.lab=1.6, cex.axis=1.5, cex.main=1.6)
lines(sd,  col=brewer.pal(10,"RdYlBu")[9], lwd=2, lty=3)

plot(Sdrt_OCUP1.3, type='l', ylim=c(0,0.41), lwd=2, main='OC3', ylab='', xlab='Province', 
     cex.lab=1.6, cex.axis=1.5, cex.main=1.6)
lines(sd,  col=brewer.pal(10,"RdYlBu")[9], lwd=2, lty=3)

plot(Sdrt_OCUP1.4, type='l', ylim=c(0,0.41), lwd=2, main='OC4', ylab='', xlab='Province', 
     cex.lab=1.6, cex.axis=1.5, cex.main=1.6)
lines(sd,  col=brewer.pal(10,"RdYlBu")[9], lwd=2, lty=3)

plot(Sdrt_OCUP1.5, type='l', ylim=c(0,0.41), lwd=2, main='OC5', ylab='', xlab='Province', 
     cex.lab=1.6, cex.axis=1.5, cex.main=1.6)
lines(sd,  col=brewer.pal(10,"RdYlBu")[9], lwd=2, lty=3)

plot(Sdrt_OCUP1.6, type='l', ylim=c(0,0.41), lwd=2, main='OC6', ylab='', xlab='Province', 
     cex.lab=1.6, cex.axis=1.5, cex.main=1.6)
lines(sd,  col=brewer.pal(10,"RdYlBu")[9], lwd=2, lty=3)

plot(Sdrt_OCUP1.7, type='l', ylim=c(0,0.41), lwd=2, main='OC7', ylab='', xlab='Province', 
     cex.lab=1.6, cex.axis=1.5, cex.main=1.6)
lines(sd,  col=brewer.pal(10,"RdYlBu")[9], lwd=2, lty=3)

plot(c(1,1), col='white', axes=F, xlab='', ylab='')
legend("center", legend=c("DSI estimates",' ',"Average contributions"), lty=c(1,1,3), col=c('black','white', brewer.pal(10,"RdYlBu")[9]),
       cex=1.6, lwd=3, bty = "n")	


# Proportions
# 0
par(mar=c(7, 4.3, 4, 4), mfrow=c(1,1))
plot(data$Ymean[(52*7*4+1):(52*7*5)], pch=19, ylim=c(0,1), xlab='Province', ylab='Proportion',
     cex.axis=1.85, cex.main=1.85, cex.lab=1.85, lwd=1.2,  xaxt = 'n')
abline(c(0.5,0), col='gray30', lwd=2, lty=2)
axis(1, at=c(26, 26+52, 26+52*2, 26+52*3, 26+52*4, 26+52*5, 26+52*6), labels=c('OC1', 'OC2', 'OC3', 'OC4', 'OC5', 'OC6', 'OC7'),
     cex.lab=1.55, cex.axis=1.47); axis(2, cex.lab=1.85, cex.axis=1.85)

# 1
par(mar=c(7, 4.3, 4, 4), mfrow=c(1,1))
plot(data$Ymean[(52*7*4+1):(52*7*5)]~data$Ocupation[(52*7*4+1):(52*7*5)], pch=19, ylim=c(0,1), xlab=' ', ylab='Proportion',
     cex.axis=1.85, cex.main=1.85, cex.lab=1.85, lwd=1.2,  xaxt = 'n')
abline(c(0.5,0), col='gray30', lwd=2, lty=2)
axis(1, at=seq(1,7, by=1), labels=c('OC1', 'OC2', 'OC3', 'OC4', 'OC5', 'OC6', 'OC7'),
     cex.lab=1.555, cex.axis=1.4); axis(2, cex.lab=1.85, cex.axis=1.85)

# 2
par(mar=c(7, 4.3, 4, 4), mfrow=c(1,1))
ind=sort(data$nd[(52*7*4+1):(52*7*5)], index.return=T)$ix
plot((sqrt(sigma2edrt[(52*7*4+1):(52*7*5)])/data$Ymean[(52*7*4+1):(52*7*5)])[ind], ylim=c(0,0.4), pch=19, col='black',
     xlab='Sample size', ylab='CV', cex.axis=1.85, cex.main=1.85, cex.lab=1.85, lwd=1.2, xaxt = 'n')
axis(1, at=round(seq(1,D*R,length=6),0), labels=floor(sort(data$nd[(52*7*4+1):(52*7*5)], index.return=T)$x)[round(seq(1,D*R,length=6),0)],
     cex.lab=1.85, cex.axis=1.85); axis(2, cex.lab=1.85, cex.axis=1.85)

# 3
par(mar=c(7, 4.3, 4, 4), xpd=T)
plot(data$Ymean[(52*7*4+1):(52*7*5)], type='l', ylim=c(0,1), xlab='Subdomain index', ylab='Proportion', col='red',
     cex.axis=1.85, cex.main=1.85, cex.lab=1.85, lwd=1.2,  xaxt = 'n')
points(mudrt1[(52*7*4+1):(52*7*5)], col='black', type='l')
new <- data.frame(x = seq(0, 365))
lines(new$x, rep(0.5,366), col='gray30', lwd=2, lty=2)
axis(1, at=c(26, 26+52, 26+52*2, 26+52*3, 26+52*4, 26+52*5, 26+52*6), labels=c('OC1', 'OC2', 'OC3', 'OC4', 'OC5', 'OC6', 'OC7'),
     cex.lab=1.85, cex.axis=1.25); axis(2, cex.lab=1.85, cex.axis=1.85)	
legend("top", inset=c(0,-0.14), legend=c("Hajek","EBLUP"), lty=c(1,1), col=c('red', 'gray30'),
       horiz=T, cex=1.7, lwd=3, bty = "n")	

# 4
par(mar=c(7, 4.3, 4, 4), xpd=T)
pred.df=data.frame('prop'=c(data$Ymean[(52*7*4+1):(52*7*5)], mudrt1[(52*7*4+1):(52*7*5)]), 
                   'occ'=rep(data$Ocupation[(52*7*4+1):(52*7*5)], 2), 'id'= c(rep(1, 364), rep(2, 364)))
boxplot(pred.df$prop ~ pred.df$id:pred.df$occ, col=c('#ff6961', 'gray80'), xaxt = 'n', 
        ylab='Proportion', cex.lab=1.85, cex.axis=1.85, xlab='OC:METHOD', pch=19)
new <- data.frame(x = seq(0, 15))
lines(new$x, rep(0.5,16), col='gray30', lwd=2, lty=2)
axis(1, at=c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5), labels=c('OC1', 'OC2', 'OC3', 'OC4', 'OC5', 'OC6', 'OC7'),
     cex.lab=1.85, cex.axis=1.25)	
legend("top", inset=c(0,-0.14), legend=c("Hajek","EBLUP"), lty=c(1,1), col=c('red', 'gray30'),
       horiz=T, cex=1.7, lwd=3, bty = "n")	

# 5
RRMSE_PI=read.csv(file='rrmse2000_PI.csv', header=T)
ind=sort(data$nd[(52*7*4+1):(52*7*5)], index.return=T)$ix
RRMSE_PI$RRMSE_PI_H = RRMSE_PI$RMSE_PI/data$Ymean

par(mar=c(4, 4.3, 4, 4), xpd=T)
plot((sqrt(data$Yvar.meanHajek[(52*7*4+1):(52*7*5)])/data$Ymean[(52*7*4+1):(52*7*5)])[ind], 
     (RRMSE_PI$RRMSE_PI_H[(52*7*4+1):(52*7*5)])[ind], pch=19, xlab='CV Hajek', ylab='RRMSE EBLUP',
     cex.axis=1.85, cex.main=1.85, cex.lab=1.85, xlim=c(0, 0.4), ylim=c(0, 0.4))
new <- data.frame(x = seq(0, 0.4, by=0.001))
lines(new$x, new$x, col='red', lwd=2, lty=2)


# Residuals
res.bt0 = data$Ymean-mudrt1  # 5 EPA
res.st0 = (res.bt0-mean(res.bt0))/sd(res.bt0)

par(mfrow=c(1,3))
boxplot(res.st0~data$Period, main='Time period', ylab='', xlab='', cex.main=1.8, cex.lab=1.6, cex.axis=1.5)

res.bt = data$Ymean[(52*7*4+1):(52*7*5)]-mudrt1[(52*7*4+1):(52*7*5)] # EPA T4
res.st = (res.bt-mean(res.bt))/sd(res.bt)

boxplot(res.st~data$Province[(52*7*4+1):(52*7*5)], main='Province', ylab='', xlab='', cex.main=1.8, cex.lab=1.6, cex.axis=1.5)
boxplot(res.st~data$Ocupation[(52*7*4+1):(52*7*5)], main='Main occupation', ylab='', xlab='', cex.main=1.8, cex.lab=1.6, 
        cex.axis=1.5, xaxt = 'n')
axis(1, at=c(1,2,3,4,5,6,7), labels=c('OC1', 'OC2', 'OC3', 'OC4', 'OC5', 'OC6', 'OC7'),
     cex.lab=1.6, cex.axis=1.25, las=2); axis(2, cex.lab=1.6, cex.axis=1.5)	



