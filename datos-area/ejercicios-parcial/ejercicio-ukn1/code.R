
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

#################
# BIVARIADO MORAN

mapa@data

moran.bi(mapa$GVA, mapa$Lbr_Prd,a.lw,zero.policy =T)
set.seed(123)

mbi_prodg <- moranbi.test(X=mapa$Lbr_Prd,Y=mapa$GVA,a.lw,zero.policy =T,adjust.n=TRUE,N=999,graph=T,print.results=T)
mbi_prod <- moranbi1.test(x=mapa$GVA,y=mapa$Lbr_Prd,a.lw, zero.policy =T,randomisation =T, alternative="two.sided",adjust.n=TRUE)

mbi_natg <- moranbi.test(X=mapa$Bsn_B_R,Y=mapa$GVA,a.lw,zero.policy =T,adjust.n=TRUE,N=999,graph=T,print.results=T)
mbi_nat <- moranbi1.test(x=mapa$Bsn_B_R,y=mapa$GVA,a.lw, zero.policy =T,randomisation =T, alternative="two.sided",adjust.n=TRUE)

##########################33
# Correlograma Bivariado

col_nbq1 <- poly2nb(mapa)
corbi.ci <- spcorrelogram.bi(col_nbq1, mapa$GVA, mapa$Lbr_Prd, a.lw ,order=3, method="I", style="W", randomisation =T, zero.policy=T,alternative="two.sided")
corbi.ci <- spcorrelogram.bi(col_nbq1, mapa$GVA, mapa$Bsn_B_R, a.lw ,order=3, method="I", style="W", randomisation =T, zero.policy=T,alternative="two.sided")
corbi.ci
plot(corbi.ci)
plot.spcorbi(corbi.ci,main="Bivariate GVA  -  NAT")

# Gráfico de dispersión del índice de Moran

x11()
mp<-moran.plot(mapa$Lbr_Prd, a.lw, main="Gráfico de Dispersión de Moran")      # Cambiar por "x1", para la estandarización
mp<-moran.plot(mapa$Bsn_B_R, a.lw, main="Gráfico de Dispersión de Moran")      # Cambiar por "x1", para la estandarización


###########################3
# LISA - BILISA

# LISA Cluster Map
moran.cluster(mapa$GVA, a.lw, zero.policy = T, mapa, significant=T)

LMCI <- localmoran.bi(mapa$GVA, mapa$Lbr_Prd, a.lw, zero.policy =T)
LMCH <- localmoran.bi(col.poly@data$CRIME, col.poly@data$HOVAL, a.lwq1, zero.policy =T)
