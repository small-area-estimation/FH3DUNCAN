
## This R code corresponds to the sections "Real data simulation experiments" and
## ''Prediction of sex occupational segregation by Spanish province'' of the paper 
## Model-based estimation of small area dissimilarity indexes: An application to sex occupational segregation in Spain.
## This script creates the files containing the aggregated area-level variables and estimates the design-based variances

## AUTHORS: Bugallo M., Esteban M.D., Morales D., Pagliarella M. C.,


rm(list=ls())
library(car)

EPA_2019T4=read.csv("EPA_2019T4.csv", sep='\t', header=T) 
names(EPA_2019T4)[names(EPA_2019T4) == 'EDAD5'] <- 'EDAD1' 
EPA_2019T4=EPA_2019T4[EPA_2019T4$EDAD1>15 & !is.na(EPA_2019T4$OCUP1),] 
EPA_2019T4$CICLO=1 
EPA_2019T4=EPA_2019T4[, c(1, 3, 7, 9, 14, 17, 20, 35, 36, 37, 48, 93)] 

EPA_2020T1=read.csv("EPA_2020T1.csv", sep='\t', header=T) 
names(EPA_2020T1)[names(EPA_2020T1) == 'EDAD5'] <- 'EDAD1' 
EPA_2020T1=EPA_2020T1[EPA_2020T1$EDAD1>15 & !is.na(EPA_2020T1$OCUP1),] 
EPA_2020T1$CICLO=2 
EPA_2020T1=EPA_2020T1[, c(1, 3, 7, 9, 14, 17, 20, 35, 36, 37, 48, 93)] 

EPA_2020T2=read.csv("EPA_2020T2.csv", sep='\t', header=T) 
names(EPA_2020T2)[names(EPA_2020T2) == 'EDAD5'] <- 'EDAD1' 
EPA_2020T2=EPA_2020T2[EPA_2020T2$EDAD1>15 & !is.na(EPA_2020T2$OCUP1),] 
EPA_2020T2$CICLO=3 
EPA_2020T2=EPA_2020T2[, c(1, 3, 7, 9, 14, 17, 20, 35, 36, 37, 48, 93)] 

EPA_2020T3=read.csv("EPA_2020T3.csv", sep='\t', header=T) 
names(EPA_2020T3)[names(EPA_2020T3) == 'EDAD5'] <- 'EDAD1' 
EPA_2020T3=EPA_2020T3[EPA_2020T3$EDAD1>15 & !is.na(EPA_2020T3$OCUP1),] 
EPA_2020T3$CICLO=4 
EPA_2020T3=EPA_2020T3[, c(1, 3, 7, 9, 14, 17, 20, 35, 36, 37, 48, 93)] 

EPA_2020T4=read.csv("EPA_2020T4.csv", sep='\t', header=T) 
names(EPA_2020T4)[names(EPA_2020T4) == 'EDAD5'] <- 'EDAD1' 
EPA_2020T4=EPA_2020T4[EPA_2020T4$EDAD1>15 & !is.na(EPA_2020T4$OCUP1),] 
EPA_2020T4$CICLO=5 
EPA_2020T4=EPA_2020T4[, c(1, 3, 7, 9, 14, 17, 20, 35, 36, 37, 48, 93)] 

EPA_2021T1=read.csv("EPA_2021T1.csv", sep='\t', header=T) 
EPA_2021T1=EPA_2021T1[EPA_2021T1$EDAD1>15 & !is.na(EPA_2021T1$OCUP1),] 
EPA_2021T1$CICLO=6 
EPA_2021T1=EPA_2021T1[, c(1, 3, 7, 9, 14, 17, 20, 34, 35, 36, 47, 91)] 

EPA_2021T2=read.csv("EPA_2021T2.csv", sep='\t', header=T) 
EPA_2021T2=EPA_2021T2[EPA_2021T2$EDAD1>15 & !is.na(EPA_2021T2$OCUP1),] 
EPA_2021T2$CICLO=7 
EPA_2021T2=EPA_2021T2[, c(1, 3, 7, 9, 14, 17, 20, 34, 35, 36, 47, 91)] 

EPA_2021T3=read.csv("EPA_2021T3.csv", sep='\t', header=T) 
EPA_2021T3=EPA_2021T3[EPA_2021T3$EDAD1>15 & !is.na(EPA_2021T3$OCUP1),] 
EPA_2021T3$CICLO=8 
EPA_2021T3=EPA_2021T3[, c(1, 3, 7, 9, 14, 17, 20, 34, 35, 36, 47, 91)] 

EPA_2021T4=read.csv("EPA_2021T4.csv", sep='\t', header=T) 
EPA_2021T4=EPA_2021T4[EPA_2021T4$EDAD1>15 & !is.na(EPA_2021T4$OCUP1),] 
EPA_2021T4$CICLO=9 
EPA_2021T4=EPA_2021T4[, c(1, 3, 7, 9, 14, 17, 20, 34, 35, 36, 47, 91)] 

EPA_agreggatted=data.frame(rbind(EPA_2019T4, EPA_2020T1, EPA_2020T2, EPA_2020T3, 
			EPA_2020T4, EPA_2021T1, EPA_2021T2, EPA_2021T3, EPA_2021T4)) 


# 'Military occupations' (0) & 'Technicians and Support Professionals' (3)
EPA_agreggatted$OCUP1[EPA_agreggatted$OCUP1==0]=3

# 'Artisans and skilled workers in the manufacturing, construction, mining industries' (7)
# & 'Plant and machinery operators and assemblers' (8)
EPA_agreggatted$OCUP1[EPA_agreggatted$OCUP1==8]=7

# 'Unskilled workers' (9) & 'Skilled workers in the agricultural, livestock, 
# forestry and fishing sectors' (6)
EPA_agreggatted$OCUP1[EPA_agreggatted$OCUP1==9]=6

sex_man<-as.numeric(EPA_agreggatted$SEXO1==1) 


################################################################################

dir2 <- function(data, w, domain, Nd) { if(is.vector(data)){
	last <- length(domain) + 1
	Nd.hat <- aggregate(w, by=domain, sum)[,last]
	nd <- aggregate(rep(1, length(data)), by=domain, sum)[,last] 
	Sum <- aggregate(w*data, by=domain, sum)
	mean <- Sum[,last]/Nd.hat
	dom <- as.numeric(Reduce("paste0", domain)) 
	if(length(domain)==1){
		domain.unique <- sort(unique(dom)) }
	else{
		domain.unique <- as.numeric(Reduce("paste0", Sum[,1:length(domain)]))
	}
	difference <- list() 
	for(d in 1:length(mean)){
		condition <- dom==domain.unique[d]
		difference[[d]] <- w[condition]*(w[condition]-1)*(data[condition]-mean[d])^2 }
		var.mean <- unlist(lapply(difference, sum))/Nd.hat^2 
		if(missing(Nd)){
			return(data.frame(Sum[,-last], mean, var.mean, Nd.hat, nd)) }
		else{
			tot <- mean*Nd
			var.tot <- var.mean*Nd^2
			return(data.frame(Sum[,-last], tot, var.tot, mean, var.mean, Nd.hat, Nd, nd))
			} }
		else{
			warning("Only a numeric or integer vector must be called as data", call. = FALSE)
} }


dir_SAE=dir2(data=sex_man, EPA_agreggatted$FACTOREL, domain=list(EPA_agreggatted$PROV, 
	EPA_agreggatted$OCUP1, EPA_agreggatted$CICLO))
names(dir_SAE)=c("Province", "Ocupation", "Period", "mean", "var.mean", "Nd.hat", "nd")

dir_SAE<-dir_SAE[dir_SAE$Period %in% 5:9, ]


################################################################################

# 1. Age4
age4=cut(EPA_agreggatted$EDAD1, breaks=c(0,30,50,Inf), include.lowest=TRUE,labels=c("1","2","3"))
age41<-as.numeric(age4==1) # '16-30'
age42<-as.numeric(age4==2) # '30-50'
age43<-as.numeric(age4==3) # '>50'

# 2. Citizenship
cit1<-as.numeric(EPA_agreggatted$NAC1==1) 

# 3. Education level
edu <- recode(EPA_agreggatted$NFORMA, " c('AN','P1','P2') = 1; c('S1') = 2; c('SG', 'SP') = 3; c('SU') = 4 ", 
	as.factor=TRUE, levels = c(1,2,3,4))
edu1<-as.numeric(edu==1) # 'primary'
edu2<-as.numeric(edu==2) # 'secondary'
edu3<-as.numeric(edu==3) # 'high secondary'
edu4<-as.numeric(edu==4) # 'superior'

# 4. Working hours
jor1<-as.numeric(EPA_agreggatted$PARCO1==1) 

# 5. Professional status 
situ <- recode(EPA_agreggatted$SITU, " c(1,3) = 1; c(5,6) = 2; 7=3; 8=4; 9=5", 
	as.factor=TRUE, levels = c(1,2,3,4,5))
situ1<-as.numeric(situ==1) 
situ2<-as.numeric(situ==2)
situ3<-as.numeric(situ==3) 
situ4<-as.numeric(situ==4)
situ5<-as.numeric(situ==5)
	

features_dir<-features_var<-hatNd<-list()

for(i in 5:9){
	ciclo.aux<-EPA_agreggatted$CICLO %in% c(i-4, i)
	EPA_aux<-EPA_agreggatted[ciclo.aux, ]
	age41.dir2<-dir2(data=age41[ciclo.aux], EPA_aux$FACTOREL, domain=list(EPA_aux$PROV, EPA_aux$OCUP1))[,1:5]
	age43.dir2<-dir2(data=age43[ciclo.aux], EPA_aux$FACTOREL, domain=list(EPA_aux$PROV, EPA_aux$OCUP1))[,3:4]
	cit1.dir2<-dir2(data=cit1[ciclo.aux], EPA_aux$FACTOREL, domain=list(EPA_aux$PROV, EPA_aux$OCUP1))[,3:4]
	edu1.dir2<-dir2(data=edu1[ciclo.aux], EPA_aux$FACTOREL, domain=list(EPA_aux$PROV, EPA_aux$OCUP1))[,3:4]
	edu3.dir2<-dir2(data=edu3[ciclo.aux], EPA_aux$FACTOREL, domain=list(EPA_aux$PROV, EPA_aux$OCUP1))[,3:4]
	edu4.dir2<-dir2(data=edu4[ciclo.aux], EPA_aux$FACTOREL, domain=list(EPA_aux$PROV, EPA_aux$OCUP1))[,3:4]
	jor1.dir2<-dir2(data=jor1[ciclo.aux], EPA_aux$FACTOREL, domain=list(EPA_aux$PROV, EPA_aux$OCUP1))[,3:4]
	situ1.dir2<-dir2(data=situ1[ciclo.aux], EPA_aux$FACTOREL, domain=list(EPA_aux$PROV, EPA_aux$OCUP1))[,3:4]
	situ2.dir2<-dir2(data=situ2[ciclo.aux], EPA_aux$FACTOREL, domain=list(EPA_aux$PROV, EPA_aux$OCUP1))[,3:4]
	situ3.dir2<-dir2(data=situ3[ciclo.aux], EPA_aux$FACTOREL, domain=list(EPA_aux$PROV, EPA_aux$OCUP1))[,3:4]
	situ4.dir2<-dir2(data=situ4[ciclo.aux], EPA_aux$FACTOREL, domain=list(EPA_aux$PROV, EPA_aux$OCUP1))[,3:4]
	
	features_dir[[i]]<-data.frame('Province'=age41.dir2[,1], "Ocupation" = age41.dir2[,2], 
							"Period" = rep(i, length(age41.dir2[,1])), 'age41' = age41.dir2[,3], 'age43' = age43.dir2[,1], 
							"cit1" = cit1.dir2[,1], "edu1" = edu1.dir2[,1], "edu3" = edu3.dir2[,1], "edu4" = edu4.dir2[,1], 
							"jor1" = jor1.dir2[,1], "situ1" = situ1.dir2[,1], "situ2" = situ2.dir2[,1], 
							"situ3" = situ3.dir2[,1], "situ4" = situ4.dir2[,1])
							 
	features_var[[i]]<-data.frame('Province'=age41.dir2[,1], "Ocupation" = age41.dir2[,2], 
							"Period" = rep(i, length(age41.dir2[,1])), 'age41.var' = age41.dir2[,4], 'age43.var' = age43.dir2[,2], 
							"cit1.var" = cit1.dir2[,2], "edu1.var" = edu1.dir2[,2],  "edu3.var" = edu3.dir2[,2], 
							"edu4.var" = edu4.dir2[,2], "jor1.var" = jor1.dir2[,2], "situ1.var" = situ1.dir2[,2],
							 "situ2.var" = situ2.dir2[,2], "situ3.var" = situ3.dir2[,2], "situ4.var" = situ4.dir2[,2] )	
	hatNd[[i]]<-	age41.dir2[,5]				 					 
}

features_dir <- do.call(rbind, features_dir)
features_var <- do.call(rbind, features_var)
hatNd <- unlist(hatNd)

features_var <- data.frame(dir_SAE[, 1:3], features_var[,4:14])	
features_var$Period <- features_var$Period - 4
	
write.csv(features_var, file="var_covariates.csv", row.names = F)

features_dir <- data.frame(dir_SAE[, 1:3], 'nd' = dir_SAE$nd, 'hatNd' = hatNd, 'intercept' = rep(1, length(dir_SAE$Nd.hat)), 
				features_dir[,4:14], 'Ymean' = dir_SAE$mean, 'Yvar.mean' = dir_SAE$var.mean)				
features_dir$Period <- features_dir$Period - 4		
features_dir$Yvar.mean[features_dir$Yvar.mean==0] = min(features_dir$Yvar.mean[features_dir$Yvar.mean!=0])			

attach(features_dir)

# GVF method
GVF.mod <- lm(log(Yvar.mean) ~ 1 + Ymean + nd)
v <- deviance(GVF.mod)/df.residual(GVF.mod)
var.mean2<-exp(v/2)*exp(predict(GVF.mod))

features_dir$Yvar.mean<-var.mean2
features_dir$Yvar.meanHajek<-dir_SAE$var.mean

write.csv(features_dir, file="YHajek_covariates.csv", row.names = F)
	





