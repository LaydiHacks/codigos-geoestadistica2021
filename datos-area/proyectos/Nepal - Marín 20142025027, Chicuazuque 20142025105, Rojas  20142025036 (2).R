library(stargazer)
library(spgwr)
library(adespatial)
library(raster)
library(tmap)
library(lmtest)
library(classInt)
library(spdep)
library(sphet)
library(pgirmess)
library(maptools)
library(ggplot2)
library(xtable)
library(RcmdrPlugin.epack)
library(olsrr)
library(RColorBrewer)
source("C:/Geoestadistica/Funciones/spcorrelogram.bi.R")
source("C:/Geoestadistica/Funciones/moranbi1.test.R")
source("C:/Geoestadistica/Funciones/moran1.bi.R")
source("C:/Geoestadistica/Funciones/moranbi.plot.R")
source("C:/Geoestadistica/Funciones/moran.bi.R")
source("C:/Geoestadistica/Funciones/moran.cluster.R")
source("C:/Geoestadistica/Funciones/getis.cluster.R")
nepal<-readShapePoly("C:/Geoestadistica/Proyecto Final/nepal/nepal/Nepal.shp")

##ESTADISTICAS DESCRIPTIVAS##
#Histogramas
par(family = "arial")
DEPECPROVhist<-hist(nepal@data$DEPECPROV,  main="" , breaks="Sturges" , col= "purple", xlab="DEPECPROV" , ylab="Frecuencia")
POVINDEXhist<-hist(nepal@data$POVINDEX,  main="" , breaks="Sturges" , col= rgb(0.2,0.8,0.5,0.5), xlab="POVINDEX" , ylab="Frecuencia")
PCINChist<-hist(nepal@data$PCINC,  main="" , breaks="Sturges" , col= rgb(1,0,0,0.5), xlab="PCINC" , ylab="Frecuencia")
PCINCPPPhist<-hist(nepal@data$PCINCPPP,  main="" , breaks="Sturges" , col= rgb(0,0,1,0.5), xlab="PCINCPPP" , ylab="Frecuencia")
PCINCMPhist<-hist(nepal@data$PCINCMP,  main="" , breaks="Sturges" , col= rgb(0.2,0.8,0.5,0.5), xlab="PCINCMP" , ylab="Frecuencia")
MALKIDShist<-hist(nepal@data$MALKIDS,  main="" , breaks="Sturges" , col= rgb(1,0,0,0.5), xlab="MALKIDS" , ylab="Frecuencia")
LIF40hist<-hist(nepal@data$LIF40,  main="" , breaks="Sturges" , col= rgb(0,0,1,0.5), xlab="LIF40" , ylab="Frecuencia")
NOSAFH20hist<-hist(nepal@data$NOSAFH20,  main="" , breaks="Sturges" , col= rgb(0.2,0.8,0.5,0.5), xlab="NOSAFH20" , ylab="Frecuencia")
AD_ILLIThist<-hist(nepal@data$AD_ILLIT,  main="" , breaks="Sturges" , col= rgb(1,0,0,0.5), xlab="AD_ILLIT" , ylab="Frecuencia")

#Boxplot
boxplot(nepal@data$DEPECPROV,col="seagreen3",ylab="DEPECPROV")
boxplot(nepal@data$PCINC,col="paleturquoise1",ylab="PCINC")

#Medidas de Tendencia central
S1<-matrix(c(SDEPECPROV,SPOVINDEX,SPCINCPPP,SPCINCMP,SMALKIDS,SLIF40,SNOSAFH20,SAD_ILLIT),5,8)
rownames(S1)<-c("Mediana","Media","Mínimo","Máximo","Desviación Estándar")
colnames(S1)<-c("DEPECPROV","POVINDEX","PCINCPPP","PCINCMP","MALKIDS","LIF40","NOSAFH20","AD_ILLIT")
xtable(t(S1))
t(S1)
SDEPECPROV<-c(median(nepal@data$DEPECPROV),mean(nepal@data$DEPECPROV),min(nepal@data$DEPECPROV),max(nepal@data$DEPECPROV), sd(nepal@data$DEPECPROV))
SPOVINDEX<-c(median(nepal@data$POVINDEX),mean(nepal@data$POVINDEX),min(nepal@data$POVINDEX),max(nepal@data$POVINDEX), sd(nepal@data$POVINDEX))
SPCINCPPP<-c(median(nepal@data$PCINCPPP),mean(nepal@data$PCINCPPP),min(nepal@data$PCINCPPP),max(nepal@data$PCINCPPP), sd(nepal@data$PCINCPPP))
SPCINCMP<-c(median(nepal@data$PCINCMP),mean(nepal@data$PCINCMP),min(nepal@data$PCINCMP),max(nepal@data$PCINCMP), sd(nepal@data$PCINCMP))
SMALKIDS<-c(median(nepal@data$MALKIDS),mean(nepal@data$MALKIDS),min(nepal@data$MALKIDS),max(nepal@data$MALKIDS), sd(nepal@data$MALKIDS))
SLIF40<-c(median(nepal@data$LIF40),mean(nepal@data$LIF40),min(nepal@data$LIF40),max(nepal@data$LIF40), sd(nepal@data$LIF40))
SNOSAFH20<-c(median(nepal@data$NOSAFH20),mean(nepal@data$NOSAFH20),min(nepal@data$NOSAFH20),max(nepal@data$NOSAFH20), sd(nepal@data$NOSAFH20))
SAD_ILLIT<-c(median(nepal@data$AD_ILLIT),mean(nepal@data$AD_ILLIT),min(nepal@data$AD_ILLIT),max(nepal@data$AD_ILLIT), sd(nepal@data$AD_ILLIT))

#QQ-Plot
qqnorm(nepal@data$DEPECPROV)
qqline(nepal@data$DEPECPROV)

##MAPA DE DEPECPROV##
tm_shape(nepal)+
  tm_polygons(("DEPECPROV"), 
      style="kmeans",
      palette="Reds",
      auto.palette.mapping=TRUE,
      title="DEPECPROV")+
  tm_borders(col= "grey40", lwd = 1, lty = "solid", alpha = 0)+
tm_format_Europe_wide() + 
  tm_style()+
  tm_compass(north = 0, type = "arrow", fontsize = 0.8, size = NA,
             show.labels = 1, cardinal.directions = c("N", "E", "S", "W"),
             text.color = NA, color.dark = NA, color.light = NA, lwd = 1,
             position =c("RIGHT" ,"TOP") , just = NA)

##ELECCIÓN MATRIZ DE PESOS ESPACIALES##
#Reina
WQueen1<-read.gal("C:/Geoestadistica/Proyecto Final/nepal/Pesos/WQueen1.gal")
WQueen2<-read.gal("C:/Geoestadistica/Proyecto Final/nepal/Pesos/WQueen2.gal")
WQueen3<-read.gal("C:/Geoestadistica/Proyecto Final/nepal/Pesos/WQueen3.gal")
WQueen4<-read.gal("C:/Geoestadistica/Proyecto Final/nepal/Pesos/WQueen4.gal")
WQueen5<-read.gal("C:/Geoestadistica/Proyecto Final/nepal/Pesos/WQueen5.gal")
WQueen6<-read.gal("C:/Geoestadistica/Proyecto Final/nepal/Pesos/WQueen6.gal")

#Torre
WRook1<-read.gal("C:/Geoestadistica/Proyecto Final/nepal/Pesos/WRook1.gal")
WRook2<-read.gal("C:/Geoestadistica/Proyecto Final/nepal/Pesos/WRook2.gal")
WRook3<-read.gal("C:/Geoestadistica/Proyecto Final/nepal/Pesos/WRook3.gal")
WRook4<-read.gal("C:/Geoestadistica/Proyecto Final/nepal/Pesos/WRook4.gal")
WRook5<-read.gal("C:/Geoestadistica/Proyecto Final/nepal/Pesos/WRook5.gal")

#K-Vecinos
coordenadas<-coordinates(nepal)
k1 <- knn2nb(knearneigh(coordenadas))
k2 <- knn2nb(knearneigh(coordenadas,2))
k3 <- knn2nb(knearneigh(coordenadas,3))
k4 <- knn2nb(knearneigh(coordenadas,4))
k5 <- knn2nb(knearneigh(coordenadas,5))

all.linked<- max(unlist(nbdists(k4, coordenadas)))
kk4 <- dnearneigh(coordenadas, 0, all.linked)
plot(nepal, border="gray")
plot(kk4, coordenadas, add=TRUE)
title(main=paste("4-vecinos más cercanos "))

#Distancia de Gabriel 
x11()
gabrielnb=graph2nb(gabrielneigh(coordenadas),sym=TRUE)
plot(nepal,border="gray")
plot(gabrielnb,coordenadas,add=T,col="red")
title(main="Gráfica de Gabriel")

#Triangulación Delaunay
trinb=tri2nb(coordenadas)
plot(nepal,border="gray")
plot(trinb,coordenadas,add=T,col="blue")
title(main="Triangulación Delaunay")

##MATRIZ DE PESOS SELECCIONADA Y VECINOS##
access.nepal<- poly2nb(nepal)
VecinosN<- nblag(access.nepal,8)
orden<-nb2listw(WQueen2, style="W", zero.policy =T)
VecinosN
orden1<-nb2listw(VecinosN[[1]], style="W", zero.policy =T)
orden2<-nb2listw(VecinosN[[2]], style="W", zero.policy =T)
orden5<-nb2listw(VecinosN[[5]], style="W", zero.policy =T)
orden7<-nb2listw(VecinosN[[7]], style="W", zero.policy =T)

##CORRELOGRAMA MORAN A PARTIR DE LA MATRIZ DE CONTIGUIDAD ESPACIAL##

cor.nepal<- sp.correlogram(WQueen6,nepal@data$DEPECPROV,order=5,method="I",style="W", zero.policy=T)
plot(cor.nepal)

##Correlogramas Bivariados##
#Matriz de Pesos Orden 6
cor.nepalbi1<- spcorrelogram.bi(WQueen6,nepal@data$DEPECPROV,nepal@data$POVINDEX,order=5,method="I",style="W", zero.policy=T)
cor.nepalbi2<- spcorrelogram.bi(WQueen6,nepal@data$DEPECPROV,nepal@data$PCINC,order=5,method="I",style="W", zero.policy=T)
cor.nepalbi3<- spcorrelogram.bi(WQueen6,nepal@data$DEPECPROV,nepal@data$NOSAFH20,order=5,method="I",style="W", zero.policy=T)
cor.nepalbi4<- spcorrelogram.bi(WQueen6,nepal@data$DEPECPROV,nepal@data$AD_ILLIT,order=5,method="I",style="W", zero.policy=T)
plot(cor.nepalbi1)
plot(cor.nepalbi2)
plot(cor.nepalbi3)
plot(cor.nepalbi4)

##I de Moran Global Variable Respuesta##
set.seed(123)
moran.test(x=nepal@data$DEPECPROV,orden,zero.policy =T,randomisation =T,alternative = "two.sided")
moran.plot(x=x1,quiet=F,zero.policy =T,listw=orden,xlab = "DEPECPROV", ylab = "W DEPECPROV")

#I de Moran Global Variables Independientes##
moran.test(x=nepal@data$POVINDEX,orden,zero.policy =T,randomisation =T,alternative = "two.sided")
moran.test(x=nepal@data$AD_ILLIT,orden,zero.policy =T,randomisation =T,alternative = "two.sided")
moran.test(x=nepal@data$NOSAFH20,orden,zero.policy =T,randomisation =T,alternative = "two.sided")
moran.test(x=nepal@data$PCINC,orden,zero.policy =T,randomisation =T,alternative = "two.sided")
moran.test(x=nepal@data$SCHOOLCNT,orden,zero.policy =T,randomisation =T,alternative = "two.sided")
moran.test(x=nepal@data$GIRLG1_5,orden,zero.policy =T,randomisation =T,alternative = "two.sided")
moran.test(x=nepal@data$BOYG1_5,orden,zero.policy =T,randomisation =T,alternative = "two.sided")
moran.test(x=nepal@data$MALKIDS,orden,zero.policy =T,randomisation =T,alternative = "two.sided")

##I de Moran Bivariado##
IB1<-moranbi1.test(x=nepal@data$DEPECPROV,y=nepal@data$POVINDEX,orden,zero.policy =T,randomisation =T,alternative = "two.sided")
IB2<-moranbi1.test(x=nepal@data$DEPECPROV,y=nepal@data$PCINC,orden,zero.policy =T,randomisation =T,alternative = "two.sided")
IB3<-moranbi1.test(x=nepal@data$DEPECPROV,y=nepal@data$NOSAFH20,orden,zero.policy =T,randomisation =T,alternative = "two.sided")
IB4<-moranbi1.test(x=nepal@data$DEPECPROV,y=nepal@data$AD_ILLIT,orden,zero.policy =T,randomisation =T,alternative = "two.sided")
IB4

##Dispersogramas de I de Moran Bivariado##
DEPECPROV <- as.data.frame(nepal)$DEPECPROV
POVINDEX<- as.data.frame(nepal)$POVINDEX
AD_ILLIT<- as.data.frame(nepal)$AD_ILLIT
NOSAFH20<- as.data.frame(nepal)$NOSAFH20
PCINC<- as.data.frame(nepal)$PCINC
x<-c(DEPECPROV)
y<-c(POVINDEX)
w<-c(NOSAFH20)
z<-c(AD_ILLIT)
u<-c(PCINC)
x1<-(x-mean(x))/sd(x)
y1<-(y-mean(y))/sd(y)
w1<-(w-mean(w))/sd(w)
z1<-(z-mean(z))/sd(z)
u1<-(u-mean(u))/sd(u)

moranbi.plot(x=x1,y=y1,quiet=F,zero.policy =T,listw=orden,xlab = "DEPECPROV", ylab = "W POVINDEX" )
moranbi.plot(x=x1,y=w1,quiet=F,zero.policy =T,listw=orden,xlab = "DEPECPROV", ylab = "W NOSAFH20" )
moranbi.plot(x=x1,y=z1,quiet=F,zero.policy =T,listw=orden,xlab = "DEPECPROV", ylab = "W AD_ILLIT" )
moranbi.plot(x=x1,y=u1,quiet=F,zero.policy =T,listw=orden,xlab = "DEPECPROV", ylab = "W PCINC" )

##I de Moran Local##
x11()
moran.cluster(nepal@data$DEPECPROV,orden,zero.policy = T, nepal, significant=T)
moran.cluster(nepal@data$POVINDEX,orden,zero.policy = T, nepal, significant=T)
moran.cluster(nepal@data$NOSAFH20,orden,zero.policy = T, nepal, significant=T)
moran.cluster(nepal@data$PCINC,orden,zero.policy = T, nepal, significant=T)
moran.cluster(nepal@data$GIRLG1_5,orden,zero.policy = T, nepal, significant=T)
moran.cluster(nepal@data$MALKIDS,orden,zero.policy = T, nepal, significant=T)

localmoran<-localmoran(nepal@data$DEPECPROV,orden,zero.policy =T)

## G de Getis Local##

localGestis<-(localG(nepal@data$DEPECPROV,orden))
localGestis
plot(localGestis)
x11()
getis.cluster(nepal@data$DEPECPROV,orden, zero.policy = T, nepal, significant=T)

class(localGestis)

##MODELOS##
nepal.data<-as.data.frame(nepal)

#MODELO CLÁSICO
mclasico<-lm(DEPECPROV~POVINDEX+PCINC+MALKIDS+LIF40+NOSAFH20+AD_ILLIT+POPULATION+KIDS1_5,data=nepal.data)
summary(mclasico)
mejorclasico<-lm(DEPECPROV~ POVINDEX+PCINC+LIF40+NOSAFH20+AD_ILLIT+lat+lon, data = nepal.data)
summary(mejorclasico)
residuosclasico<-residuals(mejorclasico)
ols_test_normality(residuosclasico)
ols_test_breusch_pagan(mejorclasico)
ols_vif_tol(mejorclasico)
lm.morantest(mejorclasico,orden)
stepwise(mclasico,"backward/forward","AIC")

#MULTIPLICADORES DE LAGRANGE
LM<-lm.LMtests(mejorclasico, orden, test="all")
print(LM)

##MODELO SPATIAL LAG
msl<-lagsarlm(DEPECPROV~POVINDEX+LIF40+MALKIDS+NOSAFH20+AD_ILLIT+lat+lon,data=nepal.data,listw = orden,method="eigen")
summary(msl, Nagelkerke=T, correlation=TRUE)
rsml<-residuals.sarlm(msl)
ols_test_normality(residuals.sarlm(msl))
moran.test(x=rsml,orden,zero.policy =T,randomisation =T,alternative = "two.sided")

#MEJOR SPATIAL LAG
msl<-lagsarlm(DEPECPROV~ POVINDEX+PCINC+LIF40+NOSAFH20+AD_ILLIT+lat+lon,data=nepal.data,listw = orden,method="eigen")
summary(msl, Nagelkerke=T, correlation=TRUE)

  #Supuestos
rsml<-residuals.sarlm(msl)
  #-Normalidad
ols_test_normality(residuals.sarlm(msl))

qqnorm(rsml, ylab="Residuos", xlab="Cuantiles teóricos",main="",pch=20,col=rgb(0.3,0.5,1,0.4))
qqline(rsml,lwd=1.9)
boxplot(rsml, data=nepal.data,col=rgb(0.3,0.5,1,0.4),alpha=0.1,main=" ",horizontal = T, xlab="Residuos",pch=20)
hist(rsml,  main="" , breaks="Sturges" , col= rgb(0.3,0.5,1,0.4), xlab="Residuos" , ylab="Frecuencia")

  #- Heterocedasticidad
bptest.sarlm(msl)
hetero.plot <- function(model) {
  plot(residuals(model) ~ fitted(model),pch=20,xlab="Predicciones", ylab="Residuales",col="blue")
  abline(h=0, lty="dotted")
  lines(lowess(fitted(model), residuals(model)), col="red")
}
hetero.plot(msl)

  #-Autocorrelación 
moran.test(x=rsml,orden,zero.policy =T,randomisation =T,alternative = "two.sided")

##MODELO SPATIAL ERROR
mse<-errorsarlm(DEPECPROV~POVINDEX+LIF40+AD_ILLIT+KIDS1_5+NOSAFH20+lat+lon, data=nepal1, listw=orden)
summary(mse, correlation=TRUE, Nagelkerke=T)
rmse<-residuals.sarlm(mse)
moran.test(x=rmse,orden,zero.policy =T,randomisation =T,alternative = "two.sided")

## SARAR - AUTOREGRESIVO EN EL ERROR
nepal.sarar <- sacsarlm(DEPECPROV~ POVINDEX+PCINC+LIF40+NOSAFH20+SCHOOLCNT+AD_ILLIT+lat+lon, data=nepal.data, listw=orden, method="eigen")    
summary(nepal.sarar, correlation=TRUE, Nagelkerke=T)
residuals.sarlm(nepal.sarar)
coef.sarlm(nepal.sarar)
bptest.sarlm(nepal.sarar)
shapiro.test(nepal.sarar$residuals)
ols_test_normality(nepal.sarar$residuals)
moran.test(x=nepal.sarar$residuals,orden,zero.policy =T,randomisation =T,alternative = "two.sided")

#MODELO DURBIN SPATIAL
nepallagsd<-lagsarlm(DEPECPROV~ POVINDEX+PCINC+LIF40+NOSAFH20+AD_ILLIT+lat+lon, listw=orden,data=nepal.data, type="mixed") 
summary(nepallagsd, correlation=TRUE)

##ANÁLISIS DE IMPACTOS##
nepal.new<-nepal
nepal.poly<-nepal
nepal_nbq1 <- poly2nb(nepal.poly) 

#Cambiando el POVINDEX
nepal.new@data[nepal.new@data$ID_3 == "74","POVINDEX"] <- 51

# Los valores de las predicciones originales
orig.pred <- as.data.frame(predict(msl))

# Los valores predichos con el nuevo índice de pobreza humana en el distrito de Udayapur
new.pred <- as.data.frame(predict(msl, newdata = nepal.new, listw = nb2listw(nepal_nbq1,style="W")))

# Las diferencias entre las predicciones
effect.51 <- new.pred$fit - orig.pred$fit                  
el <- data.frame(name = nepal.new$ID_3, dif_pred_DEPECPROV = effect.51)
nepal.new$ef51 <- el$dif_pred_DEPECPROV


# Ordenando los barrios por el valor absoluto del cambio en la predicción del DEPECPROV
el <- el[rev(order(abs(el$dif_pred_DEPECPROV))), ]
impact<-el[1:10, ]


# Mapeo de los cambios
pal5 <- brewer.pal(3, "Set1")
breaks <- c(min(nepal.new$ef51), -0.05, 0.05, max(nepal.new$ef51))
labels <- c("Efecto negativo (< -.05)", "Sin efecto (-.05 a .05)", 
            "Efecto positivo (> .05)")

np <- findInterval(nepal.new$ef51, breaks,all.inside =T) 
x11()
plot(nepal.new, col = pal5)
legend(locator(1),legend = labels, fill = pal5, bty = "n")

#Utilizando la función impacts
impacts(msl, listw=orden)

