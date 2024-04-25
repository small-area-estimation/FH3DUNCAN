
## This R code corresponds to the section "Prediction of sex occupational segregation by Spanish province" of the paper 
## Model-based estimation of small area dissimilarity indexes: An application to sex occupational segregation in Spain.
## This is the script to obtain the bootstrap estimation of mse for the model-based predictors

## AUTHORS: Bugallo M., Esteban M.D., Morales D., Pagliarella M. C.,

rm(list=ls())
set.seed(2509)

source("REML.3foldFH.indep.R")
source("BETA.U.3foldFH.indep.R")
source("stderr.3foldFH.indep.R")


########################################
######## PART 1 
########################################

data = read.csv("YHajek_covariates.csv", sep=',', header=T)
data.var.cov = read.csv("var_covariates.csv", sep=',', header=T) 

data$Province = factor(data$Province)
data$Ocupation = factor(data$Ocupation)
data$Period = factor(data$Period)

sigma2edrt <- data$Yvar.mean

# Auxiliary variables deleted 'edu3', 'age41', 'age43' and 'situ3'
X <- as.matrix(data[,c(6, 9:10, 12:15, 17)])
p <- ncol(X)

X.var <-cbind(1, data.var.cov[, c(6:7, 9:12,14)])

ydrt <- data$Ymean

D <- length(unique(data[,"Province"])) 
R <- length(unique(data[,"Ocupation"]))
T <- length(unique(data[,"Period"]))

md <- rep(R, D); mr <- rep(T, R)

########################################
######## PART 2 
########################################

# Fitting REML algorithm

sigma.0 = 0.5
fitReml <- REML.3foldFH.indep(X, ydrt, D, R, mr, T, sigma2edrt, sigma1.0 = sigma.0, 
	sigma2.0 = sigma.0, sigma3.0 = sigma.0, MAXITER = 40)

sigmau1.hat <- fitReml[[1]][1]
sigmau2.hat <- fitReml[[1]][2]
sigmau3.hat <- fitReml[[1]][3]

fitBeta <- BETA.U.3foldFH.indep(X, ydrt, D, R, mr, T, sigma2edrt, sigmau1.hat, sigmau2.hat, sigmau3.hat)
beta <- as.vector(fitBeta[[1]])

u1 <- as.vector(fitBeta[[2]])
u2 <- as.vector(fitBeta[[3]])
u3 <- as.vector(fitBeta[[4]])

########################################
######## PART 3 
########################################

### An index of dissimilarity: The Duncan Segregation Index (DSI)
## Hajek

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

## FH3 model for domain means

Ndrt = data$hatNd
mudrt1 = X%*%beta + rep(u1, 7*5) + rep(u2, 5) + u3
mudrt1 [mudrt1 >1] = 1
mudrt2 = 1 - mudrt1

mud.t1 = data.frame(aggregate(mudrt1*Ndrt, by=list(data$Province, data$Period), sum))$V1
mud.t2 = data.frame(aggregate(mudrt2*Ndrt, by=list(data$Province, data$Period), sum))$V1

Sdrt = abs(mudrt1*Ndrt/mud.t1 - mudrt2*Ndrt/mud.t2)

Sd.t_prov_period = data.frame(aggregate(Sdrt, by=list(data$Province, data$Period), sum))
Sd.t_prov_period$V1 = 1/2 *Sd.t_prov_period$V1
names(Sd.t_prov_period)=c('Province', 'Period', 'Sd.t')

# Model DSI for each province d and period time t
Sd.t = Sd.t_prov_period$Sd.t

Sdrt_OCUP1.2 = Sdrt[data$Ocupation==2]

########################################
######## PART 4 
########################################

# Bootstrap resampling

B = 2000 

y_boot=rep(NA,D*R*T)
mu1_boot=rep(NA,D*R*T)

mu1_boot_vec=matrix(nrow=B, ncol=D*R*T)
mu1_boot.hat_vec=matrix(nrow=B, ncol=D*R*T)

Sd.t_prov_period_boot = Sd.t_prov_period_boot.hat = Sd.t_prov_period_boot_OCUP1.2 =
	Sd.t_prov_period_boot.hat_OCUP1.2 = matrix(nrow=B, ncol=D*T)
	
Sd.t_prov_period_boot_DIF.hat=matrix(nrow=B, ncol=D)

sigmau1.hat_boot=rep(NA,B)
sigmau2.hat_boot=rep(NA,B)
sigmau3.hat_boot=rep(NA,B)
beta_boot=matrix(nrow=B, ncol=length(beta))

t <- proc.time() 
for (b in 1:B){
	print(b)
	### Step 2
	u1d<-rnorm(D,mean=0,sd=sqrt(sigmau1.hat)) # (a)
	u2dr<-rnorm(D*R,mean=0,sd=sqrt(sigmau2.hat)) # (b)
	u3drt<-rnorm(D*R*T,mean=0,sd=sqrt(sigmau3.hat)) # (c)
	edrt<-rnorm(D*R*T,mean=0,sd=sqrt(sigma2edrt)) # (d)
	
	vdrt<-matrix(0, nrow=D*R*T, ncol=p)
	for(i in 2:p){ vdrt[,i]<-sapply(X.var[,i], function(x){ rnorm(1, mean=0, sd=x)} ) }

	# (e)
	mu1_boot<-(X+vdrt)%*%beta+rep(u1d, R*T)+rep(u2dr, T)+u3drt 
	mu1_boot_vec[b,]<-mu1_boot
	y_boot<-mu1_boot+edrt
	
	# (f)
	mu2_boot=1-mu1_boot 
	mud.t1_boot = data.frame(aggregate(mu1_boot*Ndrt, by=list(data$Province, data$Period), sum))$V1
	mud.t2_boot = data.frame(aggregate(mu2_boot*Ndrt, by=list(data$Province, data$Period), sum))$V1
	
	Sdrt_boot = abs(mu1_boot*Ndrt/mud.t1_boot - mu2_boot*Ndrt/mud.t2_boot)
	Sd.t_prov_period_DFboot = data.frame(aggregate(Sdrt_boot, by=list(data$Province, data$Period), sum))
	Sd.t_prov_period_boot[b, ] = 1/2 *Sd.t_prov_period_DFboot$V1

	# (g)
	fitReml_boot <- REML.3foldFH.indep(X, y_boot, D, R, mr, T, sigma2edrt, sigma1.0 = sigma.0, 
		sigma2.0 = sigma.0, sigma3.0 = sigma.0, MAXITER = 40)	
	sigmau1.hat_boot[b] <- fitReml_boot[[1]][1]
	sigmau2.hat_boot[b] <- fitReml_boot[[1]][2]
	sigmau3.hat_boot[b] <- fitReml_boot[[1]][3]

	fitBeta_boot <- BETA.U.3foldFH.indep(X, y_boot, D, R, mr, T, sigma2edrt, sigmau1.hat_boot[b], sigmau2.hat_boot[b], sigmau3.hat_boot[b])
	beta_boot[b,] <- as.vector(fitBeta_boot[[1]])
	
	u1_boot <- as.vector(fitBeta_boot[[2]])
	u2_boot <- as.vector(fitBeta_boot[[3]])
	u3_boot <- as.vector(fitBeta_boot[[4]])
	
	mu1_boot.hat<-X%*%beta_boot[b,] + rep(u1_boot, 7*5) + rep(u2_boot, 5) + u3_boot
	mu1_boot.hat_vec[b,]<-mu1_boot.hat
	mu2_boot.hat<-1-mu1_boot.hat
	
	mud.t1_boot.hat = data.frame(aggregate(mu1_boot.hat*Ndrt, by=list(data$Province, data$Period), sum))$V1
	mud.t2_boot.hat = data.frame(aggregate(mu2_boot.hat*Ndrt, by=list(data$Province, data$Period), sum))$V1
	
	Sdrt_boot.hat = abs(mu1_boot.hat*Ndrt/mud.t1_boot.hat - mu2_boot.hat*Ndrt/mud.t2_boot.hat)
	Sd.t_prov_period_DFboot.hat = data.frame(aggregate(Sdrt_boot.hat, by=list(data$Province, data$Period), sum))
	Sd.t_prov_period_boot.hat[b, ] = 1/2 *Sd.t_prov_period_DFboot.hat$V1
	
	
	Sd.t_prov_period_boot_DIF.hat[b, ] = Sd.t_prov_period_boot.hat[b, 209:260] - Sd.t_prov_period_boot.hat[b, 1:52]
	
	Sd.t_prov_period_boot_OCUP1.2[b, ] = Sdrt_boot[data$Ocupation==2]
	Sd.t_prov_period_boot.hat_OCUP1.2[b, ] = Sdrt_boot.hat[data$Ocupation==2]
}

proc.time()-t 

### SAVE RESULTS

mse_prop <- apply(  (mu1_boot_vec-mu1_boot.hat_vec)^2, 2, mean, trim = 0.005 )
rmse_prop <- sqrt(mse_prop)

write.csv(data.frame('RMSE_PI'=sqrt(mse_prop)),"rrmse2000_PI.csv", row.names = FALSE)



mse_DSI<- apply(  (Sd.t_prov_period_boot-Sd.t_prov_period_boot.hat)^2, 2, mean, trim = 0.005 )
rmse_DSI <- sqrt(mse_DSI)
rrmse_DSI <- rmse_DSI/mudrt1

df=data.frame(rep(data$Province[1:52],5), c(rep(1, 52), rep(2, 52), rep(3, 52),
              rep(4, 52), rep(5, 52)), Sd.t, rmse_DSI, rrmse_DSI)
colnames(df)=c('PROV', 'PERIOD', 'PI', 'RMSE_PI','RRMSE_PI')

dff=df[df$PERIOD==5,]
write.csv(dff,"rrmse2000_DSI.2021T4.csv", row.names = FALSE)




mse_DSI_OCUP1.2 <- apply(  (Sd.t_prov_period_boot_OCUP1.2-Sd.t_prov_period_boot.hat_OCUP1.2)^2, 2, mean, trim = 0.005 )
rmse_DSI_OCUP1.2 <- sqrt(mse_DSI_OCUP1.2)
rrmse_DSI_OCUP1.2 <- rmse_DSI_OCUP1.2/Sdrt_OCUP1.2

df_OC2 <- data.frame(rep(data$Province[1:52],5), c(rep(1, 52), rep(2, 52), rep(3, 52),
              rep(4, 52), rep(5, 52)), rmse_DSI_OCUP1.2, rrmse_DSI_OCUP1.2)
              
colnames(df_OC2)=c('PROV', 'PERIOD', 'RMSE_PI_OCUP1.2','RRMSE_PI_OCUP1.2')
df_OC2 <- df_OC2[df_OC2$PERIOD==5,]
write.csv(df_OC2,"rrmse2000_OC2.csv", row.names = FALSE)



IC_Sd.t_DIF = apply(Sd.t_prov_period_boot_DIF.hat, 2, quantile, probs=c(0.025, 0.975))

df2=data.frame('PROV' = data$Province[1:52], 'LB' = IC_Sd.t_DIF[1,], 'UB' = IC_Sd.t_DIF[2,])
write.csv(df2,"PI_DuncanDIFER_IC.csv", row.names = FALSE)



