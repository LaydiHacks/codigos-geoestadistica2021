setwd("GEOESTADISTICA/codigos-geoestadistica2021/")  
library(spgwr)
library(adespatial)
library(raster)
library(tmap)
library(lmtest)
library(RColorBrewer)
library(classInt)
library(spdep)
library(sphet)
library(pgirmess)
library(rgdal)
library(spatialreg)
library(car)
source("funciones/northarrow.r")
source("funciones/scalebar.R")
source("funciones/moran.R")
source("funciones/geary.R")
source("funciones/moranbi.test.R")
source("funciones/randomize_vector.R")
source("funciones/moran.cluster.R")
source("funciones/moran.bi.R")
source("funciones/moran.cluster.R")
source("funciones/getis.cluster.R")
source("funciones/localmoran.bi.R")
source("funciones/moranbi.plot.R") 
source("funciones/quantile.e.R")
source("funciones/sp.na.omit.R")
source("funciones/correlogram.d.R")
source("funciones/sp.correlogram.R")
source("funciones/spcorrelogram.bi.R")
source("funciones/moran.bi1.R")
source("funciones/moranbi1.test.R")
source("funciones/geary.bi.R")
source("funciones/test.w.R")

#########################################################################
#################### PARCIAL 1 #########################################
##########  Laydi Viviana Baustista - 20151025069
##########  Alejandro Dimate - 20152025827

mapa <-readOGR("datos-area/ejercicios-parcial/parcial-20211/Spain.N2.shp")
proj4string(mapa) <- CRS("+proj=longlat +datum=WGS84")
plot(mapa)


###############################################
#######      ESTADISTICOS        ###############
###############################################

boxplot(mapa$GDP,col="seagreen3",ylab="TEXTODELEJE")
hist(mapa$GDP)
summary(mapa$GDP)
shapiro.test(mapa$GDP)
qqnorm(mapa@data$GDP)
qqline(mapa@data$GDP)


coords <- coordinates(mapa)
X = coords[,1]
Y = coords[,2]


#######################################################
#######      ANALISIS EXPLORATORIO       ###############
####################################################

mapa@data

## Matriz de Pesos 
k4 <-nb2listw(knn2nb(knearneigh(coords,4)), style="W")
plot(mapa,border="gray")
plot(k4,coords,add=T,col="red")
title(main="Gráfica de K4")

orden<-nb2listw(k4$neighbours, style="W", zero.policy =T)


# I DE MORAN
mp<-moran.plot(mapa$GDP,orden, main="Gráfico de Dispersión de Moran GDP")
mp<-moran.plot(mapa$RPHP,orden, main="Gráfico de Dispersión de Moran RHP")


# BIVARIADO I DE MORAN
moran.bi(mapa$GDP, mapa$RPHP,k4,zero.policy =T)
set.seed(123)
mbi_prodg <- moranbi.test(X=mapa$RPHP,Y=mapa$GDP, k4,zero.policy =T,adjust.n=TRUE,N=999,graph=T,print.results=T)


# Correlograma Bivariado

col_nbq1 <- poly2nb(mapa)
corbi.ci <- spcorrelogram.bi(col_nbq1, mapa$GDP, mapa$RPHP, a.lw ,order=3, method="I", style="W", randomisation =T, zero.policy=T,alternative="two.sided")
corbi.ci
plot(corbi.ci)
plot.spcorbi(corbi.ci,main="Bivariado GDP  -  RPHP")


# LISA Cluster Map
moran.cluster(mapa$GDP, orden, zero.policy = T, mapa, significant=T)

###################################
### MODELOS

#REINA

o.nb <- read.gal("datos-area/ejercicios-parcial/parcial-20211/Spain.N2.gal", override = TRUE)
a.lw <- nb2listw(o.nb, style="W")
plot(mapa,border="gray")
plot(a.lw,coords,add=T,col="green")

#=============================================================
############### MODELO CLASICO - LINEAL  ##############
#=============================================================

GDP <- mapa$GDP
RPHP <- mapa$RPHP

formula_modelo = GDP~RPHP
mclasico<-lm(formula_modelo,data = mapa)
mejorclasico<-mclasico
summary(mejorclasico)

# Validación supuestos
bptest(mejorclasico)
resettest(mejorclasico)
raintest(mejorclasico)
shapiro.test(residuals(mejorclasico))

residuos <- mapa
residuos$mc <-residuals(mejorclasico)
spplot(residuos, "mc", col.regions = rev(terrain.colors(20))) #Mapa de Residuos
shapiro.test(residuos$mc)  
bptest(mejorclasico)
lm.morantest(mejorclasico,a.lw) 
#Corremos todas las pruebas
pruebas <- lm.LMtests(mejorclasico, listw = a.lw, test = "all")
summary(pruebas)


######################
# Modelo Spatial Lag:                       
mls <- lagsarlm(formula_modelo, data=as.data.frame(mapa), listw=a.lw) 
summary(mls, Nagelkerke=T, correlation=TRUE)
AIC(mls)
shapiro.test(mls$residuals)
plot(mls$residuals)                                
dens.sl <- density(mls$residuals)                   
plot(dens.sl, main="Gráfico de la función de densidad de los residuos (Spatial Lag)") 
rsml<-residuals.sarlm(mls)
moran.test(x=rsml,a.lw,zero.policy =T,randomisation =T,alternative = "two.sided")


# Modelo Spatial Error:
sem <- errorsarlm(formula_modelo, data=as.data.frame(mapa), listw=a.lw)
summary(sem, Nagelkerke=T)
AIC(sem)
residuals.sarlm(sem)
coef.sarlm(sem)
bptest.sarlm(sem)
hetero.plot(sem)
shapiro.test(sem$residuals)

