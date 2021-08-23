
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


mapa <-readOGR("datos-area/ejercicios-parcial/ejercicio-ukn1/UK.N1.shp")
proj4string(mapa) <- CRS("+proj=longlat +datum=WGS84")
plot(mapa)
x11()

o.nb <- read.gal ("datos-area/ejercicios-parcial/ejercicio-ukn1/UK.N1.gal", override = TRUE)
a.lw <- nb2listw(o.nb, style="W")
?read.gal

mapa@data

# Gráfico de dispersión del índice de Moran

mp<-moran.plot(mapa$GVA,a.lw, main="Gráfico de Dispersión de Moran")

################ MODELO
GVA <- mapa$GVA
PROD <- mapa$Lbr_Prd
NAT <- mapa$Bsn_B_R

formula_modelo = GVA~PROD+NAT
mclasico<-lm(formula_modelo)
summary(mclasico)
#Sacamos las variables no significativas
formula_mc <-  GVA~PROD+NAT
mejorclasico<-lm(formula_mc)
summary(mejorclasico)

# Validación supuestos
bptest(mejorclasico)
resettest(mejorclasico)
raintest(mejorclasico)
shapiro.test(residuals(mejorclasico))
vif(mejorclasico)

residuos <- mapa
residuos$mc <-residuals(mejorclasico)
spplot(residuos, "mc", col.regions = rev(terrain.colors(20))) #Mapa de Residuos
shapiro.test(residuos$mc)  
bptest(mejorclasico)
lm.morantest(mejorclasico,a.lw) 
#Corremos todas las pruebas
pruebas <- lm.LMtests(mejorclasico, listw=a.lw, test = "all")
summary(pruebas)


