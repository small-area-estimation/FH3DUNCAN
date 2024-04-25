
## This R code corresponds to the section ''Prediction of sex occupational segregation by Spanish province'' of the paper 
## Model-based estimation of small area dissimilarity indexes: An application to sex occupational segregation in Spain.
## This script maps DSI predictions and RRMSEs

## AUTHORS: Bugallo M., Esteban M.D., Morales D., Pagliarella M. C.,

           
library(maptools)
library(RColorBrewer)

GroupClassification <- function(data,datacompare,intervals)
{
   n = length(data)
   group = matrix(0,nrow=n) 
   ninterv = length(intervals)
   for (i in 1:n)
   {
      for (j in 1:ninterv)
         if (datacompare[i]<intervals[j])
         {
            group[i]= intervals[j]
            break
         }
   }
   result = split(data,group)
   return (result)
}
 
PrintSpainMap <- function(pathmap,datos,colors,titlemap,textlegend,eliminarprov)
{

   m <- matrix(c(1,1,1,2),2,2)
   layout(m, widths=c(1.5, 1), heights=c(1.5, 1), respect=F)

   xName <- readShapePoly(pathmap, IDvar="NAME", proj4string=CRS("+proj=longlat +ellps=clrk66"))     

   xName$datos <- NA
   for (i in 1:length(colors))
      xName$datos[datos[[i]]] <- colors[i]
 
   xSC <- xName[xName$ESP_PROV_I < 35 | xName$ESP_PROV_I >38 | xName$ESP_PROV_I==36 | xName$ESP_PROV_I ==37,]
   plot(xSC,  xlab="",  col=xSC$datos, axes=F)

  title(titlemap, line=-0.5, cex.main=1.8)
  legend( "topright", textlegend, pt.bg = colors, pch=21, bty = "n", cex=1.4, pt.cex = 2.2)  #  cex=1.3

   # box()

   xC <- xName[xName$ESP_PROV_I == 35 | xName$ESP_PROV_I == 38,]
   plot(xC,  xlab="",  col=xC$datos)
   box()
}

pathmap    <- "spainmap/esp_prov.shp"

################################################################################
############################## INDICE DE DUNCAN ################################
################################################################################

# 2021 T4
datos=read.table(file="PI_Duncan.txt", header=TRUE, sep=';', dec='.')
dom=datos$PROV[1:50]

estML=datos$Sd.t_2021T4[1:50]
 
intervals_prop <- c(0, 0.15, 0.20, 0.25, 0.30, Inf)  

colorsprop <- c(brewer.pal(10,"RdYlBu")[c(6,7,8,9, 10)])

legend_prop <- expression("< 15 %", "15-20 %", "20-25 %", "25-30 %", "30-35 %") 

result = GroupClassification(dom,estML,intervals_prop)
PrintSpainMap(pathmap,result,colorsprop,"DSI estimates",legend_prop,eliminarprov)


################################################################################
# 2021 T4 - PI_DuncanR2
datos=read.table(file="PI_DuncanR2.txt", header=TRUE, sep=';', dec='.')
dom=datos$PROV[1:50]

estML=datos$Sdrt_OCUP1.2[1:50]

intervals_prop <- c(0.05, 0.10, 0.15, 0.20, 0.25)  

colorsprop <- c(brewer.pal(10,"RdYlBu")[c(7,8,9, 10)])

legend_prop <- expression("5-10 %", "10-15 %", "15-20 %", "15-25 %") 

result = GroupClassification(dom,estML,intervals_prop)
PrintSpainMap(pathmap,result,colorsprop,"DSI estimates: OC2",legend_prop,eliminarprov)


################################################################################
# 2021 T4 - PI_Duncan_sd
datos=read.table(file="PI_Duncan_sd.txt", header=TRUE, sep=';', dec='.')
dom=datos$PROV[1:50]

estML=datos$sd[1:50]
 
intervals_prop <- c(0.04, 0.05, 0.07, 0.10)  

colorsprop <- c('white', brewer.pal(10,"RdYlBu")[c(6,7)])

legend_prop <- expression("4-5 %", '5-7 %', "7-10 %") 

result = GroupClassification(dom,estML,intervals_prop)
PrintSpainMap(pathmap,result,colorsprop,"Average contribution",legend_prop,eliminarprov)

##############################################################

# 2021 T4 RRMSE
datos=read.csv(file='rrmse2000_DSI.2021T4.csv', header=T)
dom=datos$PROV[1:50]

estML=datos$RRMSE_PI[1:50]

intervals_prop <- c(0, 0.10, 0.15, 0.20, 0.25, Inf)  

colorsprop <- c(brewer.pal(10,"RdGy")[c(5, 4, 3, 2, 1)])

legend_prop <- expression("< 10 %", "10-15 %", "15-20 %", "20-25 %", "25-30 %") 

result = GroupClassification(dom,estML,intervals_prop)
PrintSpainMap(pathmap,result,colorsprop,"RRMSE",legend_prop,eliminarprov)


##############################################################

# 2021 T4 RRMSE OC2
datos=read.csv(file='rrmse2000_OC2.csv', header=T)
dom=datos$PROV[1:50]

estML=datos$RRMSE_PI_OCUP1.2[1:50]

intervals_prop <- c(0.10, 0.15, 0.20, 0.25, 0.30, Inf)  

colorsprop <- c(brewer.pal(10,"RdGy")[c(4, 3, 2, 1)], '#640b0b')

legend_prop <- expression("10-15 %", "15-20 %", "20-25 %", "25-30 %", ">30 %") 

result = GroupClassification(dom,estML,intervals_prop)
PrintSpainMap(pathmap,result,colorsprop,"RRMSE: OC2",legend_prop,eliminarprov)


################################################################################

# 2021 T4 - 2020 T4 : PREDICTIONS
datos=read.table(file="PI_DuncanDIFER.txt", header=TRUE, sep=';', dec='.')
dom=datos$PROV[1:50]

estML=datos$Sd.t_dif[1:50]
 
intervals_prop <- c(-0.10, -0.03, -0.01, 0.01, 0.03, 0.11)  

colorsprop <- c(brewer.pal(10,"RdYlBu")[c(6,7,8,9, 10)])

legend_prop <- expression("(-0.07, -0.03]", "(-0.03, -0.01]", "(-0.01, 0.01]", "(0.01, 0.03]", "(0.03, 0.10]") 

result = GroupClassification(dom,estML,intervals_prop)
PrintSpainMap(pathmap,result,colorsprop,"DSI evolution",legend_prop,eliminarprov)


################################################################################

# 2021 T4 - 2020 T4 : PREDICTIONS OC2
datos=read.table(file="PI_DuncanR2DIFER.txt", header=TRUE, sep=';', dec='.')
dom=datos$PROV[1:50]

estML=datos$Sdrt_OCUP1.2[1:50]
 
intervals_prop <- c(-0.10, -0.03, -0.01, 0.01, 0.03, 0.11)  

colorsprop <- c(brewer.pal(10,"RdYlBu")[c(6,7,8,9, 10)])

legend_prop <- expression("(-0.07, -0.03]", "(-0.03, -0.01]", "(-0.01, 0.01]", "(0.01, 0.03]", "(0.03, 0.10]") 

result = GroupClassification(dom,estML,intervals_prop)
PrintSpainMap(pathmap,result,colorsprop,"DSI evolution: OC2",legend_prop,eliminarprov)


