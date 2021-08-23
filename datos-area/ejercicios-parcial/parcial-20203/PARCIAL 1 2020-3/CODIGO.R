source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/northarrow.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/scalebar.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/moran.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/geary.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/moranbi.test.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/randomize_vector.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/moran.cluster.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/moran.bi.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/moran.cluster.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/getis.cluster.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/localmoran.bi.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/moranbi.plot.R") 
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/quantile.e.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/sp.na.omit.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/correlogram.d.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/sp.correlogram.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/spcorrelogram.bi.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/moran.bi1.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/moranbi1.test.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/geary.bi.R")
source("C:/Users/leonh/Desktop/Geoestadistica/Parcial 2020-1 Geoestadistica/Funciones/test.w.R")

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

setwd("C:/Users/leonh/Desktop/PARCIALGEO20201/primerparcial2020geoestadstica")  
EU <-readOGR(dsn="C:/Users/leonh/Desktop/PARCIALGEO20201/primerparcial2020geoestadstica",layer="mapa.UE")
proj4string(EU) <- CRS("+proj=longlat +datum=WGS84")
plot(EU)
x11()

reina.nb <- read.gal ("C:/Users/leonh/Desktop/PARCIALGEO20201/primerparcial2020geoestadstica/mapa.UE.gal", override = TRUE)
?read.gal

EU@data

# BIVARIADO MORAN

EU@data

moran.bi(EU$Cr10_11, EU$PGE_2009,a.lw,zero.policy =T)
set.seed(123)

mbi0 <- moranbi.test(X=EU$Cr10_11,Y=EU$PGE_2009,a.lw,zero.policy =T,adjust.n=TRUE,N=999,graph=T,print.results=T)

mbi <- moranbi1.test(x=EU$Cr10_11,y=EU$PGE_2009,a.lw, zero.policy =T,randomisation =T,
                     alternative="two.sided",adjust.n=TRUE)
mbi

# Correlograma Moran a partir de matriz contiguidad espacial

sp.cr <- sp.correlogram(reina.nb, EU$Cr10_11, order=9, method="I", style="W", zero.policy=T)
?sp.correlogram

x11()
plot(sp.cr, main="Correlograma Cr10_11")

sp.cr1 <- sp.correlogram(reina.nb, EU$PGE_2009, order=9, method="I", style="W", zero.policy=T)
?sp.correlogram

x11()
plot(sp.cr1, main="Correlograma PGE_2009")

# Correlograma Bivariado

corbi.ci <- spcorrelogram.bi(col_nbq1, EU$Cr10_11, EU$PGE_2009, order=9, 
                             method="I", style="W", randomisation =T, zero.policy=T,alternative="two.sided")
?spcorrelogram.bi
corbi.ci
plot(corbi.ci)
plot.spcorbi(corbi.ci,main="Bivariate GDP_10_11   -   EDUCATION % GDP_09")

# ráfico de dispersión del índice de Moran

x11()
mp<-moran.plot(EU$Cr10_11, a.lw, main="Gráfico de Dispersión de Moran")      # Cambiar por "x1", para la estandarización
mp<-moran.plot(EU$PGE_2009, a.lw, main="Gráfico de Dispersión de Moran")      # Cambiar por "x1", para la estandarización

# I moran Local

lisa_obs1<-localmoran(EU$Cr10_11, a.lw, zero.policy=T)
lisa_obs1

lisa_obs3<-localmoran(EU$PGE_2009, a.lw, zero.policy=T)
lisa_obs3

# Modelo general GNS (A partir de lo general ...?)
EU@data

mod.GNS1 <- sacsarlm(EU$Cr10_11 ~ EU$PGE_2009, listw = a.lw, data = EU@data, type = "sacmixed")
summary(mod.GNS1,Nagelkerke=T)
?sacsarlm

mod.GN01 <- sacsarlm(EU$PGE_2009 ~ EU$Cr10_11, listw = a.lw, data = EU@data, type = "sacmixed")
summary(mod.GN01,Nagelkerke=T)



