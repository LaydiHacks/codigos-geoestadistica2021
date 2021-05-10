#----------------------------------PROYECTO ECONOMETRIA ESPACIAL---------------------------------#
#---------------------------Viviana Bautista - Alejandro Dimate Rodriguez------------------------#
#------------AN?LISIS GEOESTAD?STICO DE INDICE DE VENTAJA SOCIECONOMICA NEW YORK 2010-2014--------#

library(spgwr)
library(raster)
library(tmap)
library(lmtest)
library(RColorBrewer)
library(classInt)
library(spdep)
library(sphet)
library(pgirmess)
library(spatialreg)
library(dbscan)
library(adespatial)
library(rgdal)
library(RGeostats)
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


shp <- "db/ventajasocioeconomica/condados/newyork.shp" 
mapa <- readOGR(shp)



###############################################
#######      VARIABLES        #################
###############################################


IVS <- as.data.frame(mapa)$IVS
DISC<- as.data.frame(mapa)$DISC
AGE17<- as.data.frame(mapa)$AGE17
AGE65<- as.data.frame(mapa)$AGE65
PMIN<- as.data.frame(mapa)$PMIN
NOSEGURO<- as.data.frame(mapa)$NOSEGURO
RENT30<- as.data.frame(mapa)$RENT30
RENT<- as.data.frame(mapa)$RENT
RENT <- as.numeric(RENT)
#NOAUTO<- as.data.frame(mapa)$NOAUTO
AH<- as.data.frame(mapa)$AH
EN5<- as.data.frame(mapa)$EN5
MPAREN<- as.data.frame(mapa)$MPAREN
NDS25<- as.data.frame(mapa)$NDS25
IPC<- as.data.frame(mapa)$IPC
CIVILES16<- as.data.frame(mapa)$CIVILES16
POBREZA<- as.data.frame(mapa)$POBREZA

coords <- coordinates(mapa)
X = coords[,1]
Y = coords[,2]
data  <- data.frame(mapa$IVS, X,Y)


###############################################
#######      ESTADISTICOS        ###############
###############################################

boxplot(mapa$IVS,col="seagreen3",ylab="TEXTODELEJE")
hist(mapa$IVS)
summary(mapa$IVS)
shapiro.test(mapa$IVS)
qqnorm(mapa@data$IVS)
qqline(mapa@data$IVS)



###############################################
#######     Amorfosis Gaussiana       #########
###############################################


IVSselect = data[,c("X","Y","mapa.IVS")] 
IVS.rgdb <- db.create(IVSselect,ndim=2,autoname=F)
IVS.herm <- anam.fit(IVS.rgdb,name="mapa.IVS",type="gaus")
IVS.hermtrans <- anam.z2y(IVS.rgdb,names="mapa.IVS",anam=IVS.herm)
IVS.trans <- IVS.hermtrans@items$Gaussian.mapa.IVS
shapiro.test(IVS.trans) #Normalizo



#=====================================================
############### SELECCIONAR MATRIZ W ##############
#====================================================

#REINA
  ny <- poly2nb(mapa) 
  ny.lags <- nblag(ny, 9)
  a.lwq1 <- nb2listw(ny.lags[[1]], style="W", zero.policy =T)
  a.lwq2 <- nb2listw(ny.lags[[2]], style="W", zero.policy =T)
  a.lwq3 <- nb2listw(ny.lags[[3]], style="W", zero.policy =T)
  a.lwq4 <- nb2listw(ny.lags[[4]], style="W", zero.policy =T)
  a.lwq5 <- nb2listw(ny.lags[[5]], style="W", zero.policy =T)
  a.lwq6 <- nb2listw(ny.lags[[6]], style="W", zero.policy =T)

#TORRE
  NYT <- poly2nb(mapa,queen=FALSE)
  nyt.lags <- nblag(NYT , 10)
  r.lw1  <- nb2listw(nyt.lags[[1]], style="W", zero.policy =T)
  r.lw2  <- nb2listw(nyt.lags[[2]], style="W", zero.policy =T)
  r.lw3  <- nb2listw(nyt.lags[[3]], style="W", zero.policy =T)
  r.lw4  <- nb2listw(nyt.lags[[4]], style="W", zero.policy =T)
  r.lw5  <- nb2listw(nyt.lags[[5]], style="W", zero.policy =T)
  r.lw6  <- nb2listw(nyt.lags[[6]], style="W", zero.policy =T)

#K-Vecinos
  IDs <- row.names(as(mapa, "data.frame"))
  
  ny.k1 <-nb2listw(knn2nb(knearneigh(coords)) , style="W")
  ny.k2 <-nb2listw(knn2nb(knearneigh(coords,2)), style="W")
  ny.k3 <-nb2listw(knn2nb(knearneigh(coords,3)), style="W")
  ny.k4 <-nb2listw(knn2nb(knearneigh(coords,4)), style="W")
  ny.k5 <-nb2listw(knn2nb(knearneigh(coords,5)), style="W")

#Criterios basados en Graficas

  ## Ddelaunay
    trinb=tri2nb(coords)
    delaunay <-nb2listw(trinb, style="W")
    #plot(mapa,border="gray")
    #plot(trinb,coords,add=T,col="blue")
    #title(main="Triangulaci?n Delaunay")
  
  ## Esfera de Influencia
    soinb=graph2nb(soi.graph(trinb,coords))
    esf.influencia <-nb2listw(soinb, style="W")
    #plot(mapa,border="gray")
    #plot(soinb,coords,add=T,col="green")
    #title(main="Esfera de influencia")
  
  ## Gabriel
    gabrielnb=graph2nb(gabrielneigh(coords),sym=TRUE)
    gabriel <-nb2listw(gabrielnb, style="W")
    #plot(mapa,border="gray")
    #plot(gabrielnb,coords,add=T,col="red")
    #title(main="Gr?fica de Gabriel")
  
  ## Vecinos Relativos
    relativenb=graph2nb(relativeneigh(coords),sym=TRUE)
    vec.relative <-nb2listw(relativenb, style="W")
    #plot(mapa,border="gray")
    #plot(relativenb,coords,add=T,col="orange")
    #title(main="Vecinos relativos")
    #op=par(mfrow=c(2,2)) 
    #par(op)



Pesos.list<-list(reina1=print.listw(a.lwq1,zero.policy=TRUE),
                 reina2=print.listw(a.lwq2,zero.policy=TRUE),
                 reina3=print.listw(a.lwq3,zero.policy=TRUE),
                 reina4=print.listw(a.lwq4,zero.policy=TRUE),
                 reina5=print.listw(a.lwq6,zero.policy=TRUE),
                 reina6=print.listw(a.lwq6,zero.policy=TRUE),
                 torre1=print.listw(r.lw1,zero.policy=TRUE),
                 torre2=print.listw(r.lw2,zero.policy=TRUE),
                 torre3=print.listw(r.lw3,zero.policy=TRUE),
                 torre4=print.listw(r.lw4,zero.policy=TRUE),
                 torre5=print.listw(r.lw5,zero.policy=TRUE),
                 torre6=print.listw(r.lw6,zero.policy=TRUE), 
                 kvecinos1=ny.k1,
                 kvecinos2=ny.k2,
                 kvecinos3=ny.k3,
                 kvecinos4=ny.k4,
                 kvecinos5=ny.k5,
                 gabriel=gabriel,
                 delaunay=delaunay,
                 esfera.inf=esf.influencia,
                 vec.relativos=vec.relative)



class(Pesos.list)
nbw <- length(Pesos.list)
1 - (1 - 0.05)^(nbw)
W_sel <- listw.select(gabriel$neighbours, Pesos.list, MEM.autocor = "all", p.adjust = TRUE, nperm = 50)
W_sel$candidates
W_sel$best.id


#============================================
############### CORRELOGRAMA ##############
#===========================================


correlograma_mapa <- sp.correlogram(gabrielnb, mapa@data$IVS, order=7,method="I",style="W", zero.policy=T)
plot(correlograma_mapa)

##Correlogramas Bivariados##
#Matriz de Pesos 
cor.mapabi1<- spcorrelogram.bi(gabrielnb,mapa$IVS,mapa$DISC,order=7,method="I",style="W", zero.policy=T)
plot(cor.mapabi1)


#=============================================================
############### AUTOCORRELACION ##############
#=============================================================
orden<-nb2listw(gabrielnb, style="W", zero.policy =T)

################################
####### I de Moran Global ######
set.seed(123)
moran.test(y,orden,zero.policy =T,randomisation =T,alternative = "two.sided") 
moran.plot(z_y,quiet=F,zero.policy =T,listw=orden,xlab = "IVS", ylab = "W IVS")

################################
####### I de Moran Bivariado ######

x <- mapa$POBREZA #Variables
#x <- (x-mean(x))/sd(x) #Normalizar variable
bivariado_yx<-moranbi1.test(mapa$IVS, x ,orden,zero.policy =T,randomisation =T,alternative = "two.sided")
bivariado_yx

####Dispersograma
moranbi.plot(mapa$IVS,x,quiet=F,zero.policy =T,listw=orden,xlab = "DEPECPROV", ylab = "W POVINDEX" ) 

################################
####### LISA      ##############
x11()
moran.cluster(mapa$IVS,orden,zero.policy = T, mapa, significant=T)

################################
####### GETIS     ############
localGestis<-(localG(mapa$IVS,orden))
localGestis
plot(localGestis)
x11()
getis.cluster(mapa$IVS,orden, zero.policy = T, mapa, significant=T)

class(localGestis)


#=============================================================
############### MODELOS  ##############
#=============================================================

mapa.data<-as.data.frame(mapa)
orden<-nb2listw(gabrielnb, style="W", zero.policy =T)

formula_modelo = IVS~AGE17+AGE65+AH+CIVILES16+DISC+EN5+log(IPC)+MPAREN+NDS25+NOSEGURO+PMIN+POBREZA+RENT30

#MODELO GNSS
mod.GNS <- sacsarlm(formula_modelo, listw = orden, data = mapa.data, type = "sacmixed")
summary(mod.GNS,Nagelkerke=T)
mod.GNS2 <- sacsarlm(l_salario ~ desempleo, listw = W2_list, data = spatial_data, type = "sacmixed")
summary(mod.GNS2,Nagelkerke=T)

#MODELO DURBIN SPATIAL
nepallagsd<-lagsarlm(formula_modelo, listw=orden,data=mapa.data, type="mixed") 
summary(nepallagsd, correlation=TRUE)


#MODELO CL?SICO
  mclasico<-lm(IVS~AGE17+AGE65+AH+CIVILES16+DISC+EN5+log(IPC)+MPAREN+NDS25+NOSEGURO+PMIN+POBREZA+RENT30,data = mapa.data)
  summary(mclasico)
  #Sacamos las variables no significativas
  mejorclasico<-lm(IVS~AGE65+AH+EN5+log(IPC)+MPAREN+NDS25+NOSEGURO+PMIN+POBREZA+RENT30,mapa.data)
  summary(mejorclasico)
  
  residuosclasico<-residuals(mejorclasico)
  shapiro.test(residuosclasico)  
  bptest(mejorclasico)
  lm.morantest(mejorclasico,orden)
  stepwise(mclasico,"backward/forward","AIC")

#MULTIPLICADORES DE LAGRANGE
  LM<-lm.LMtests(mejorclasico, orden, test="all")
  print(LM)

##MODELO SPATIAL LAG
  msl<-lagsarlm(IVS~AGE17+AGE65+AH+CIVILES16+DISC+EN5+log(IPC)+MPAREN+NDS25+NOSEGURO+PMIN+POBREZA9+RENT+RENT30,mapa.data,listw =orden,method="eigen")
  summary(msl, Nagelkerke=T, correlation=TRUE)
  rsml<-residuals.sarlm(msl)
  ols_test_normality(residuals.sarlm(msl))
  moran.test(x=rsml,orden,zero.policy =T,randomisation =T,alternative = "two.sided")

#MEJOR SPATIAL LAG
  msl<-lagsarlm(DEPECPROV~ POVINDEX+PCINC+LIF40+NOSAFH20+AD_ILLIT+lat+lon,data=mapa.data,listw = orden,method="eigen")
  summary(msl, Nagelkerke=T, correlation=TRUE)

##MODELO SPATIAL ERROR
  mse<-errorsarlm(DEPECPROV~POVINDEX+LIF40+AD_ILLIT+KIDS1_5+NOSAFH20+lat+lon, data=nepal1, listw=orden)
  summary(mse, correlation=TRUE, Nagelkerke=T)
  rmse<-residuals.sarlm(mse)
  moran.test(x=rmse,orden,zero.policy =T,randomisation =T,alternative = "two.sided")
  
## SARAR - AUTOREGRESIVO EN EL ERROR
  nepal.sarar <- sacsarlm(DEPECPROV~ POVINDEX+PCINC+LIF40+NOSAFH20+SCHOOLCNT+AD_ILLIT+lat+lon, data=mapa.data, listw=orden, method="eigen")    
  summary(nepal.sarar, correlation=TRUE, Nagelkerke=T)
  residuals.sarlm(nepal.sarar)
  coef.sarlm(nepal.sarar)
  bptest.sarlm(nepal.sarar)
  shapiro.test(nepal.sarar$residuals)
  ols_test_normality(nepal.sarar$residuals)
  moran.test(x=nepal.sarar$residuals,orden,zero.policy =T,randomisation =T,alternative = "two.sided")
  

#=============================================================
############### MODELOS - SUPUESTOS  ##############
#=============================================================  
  
  
#Supuestos
  rsml<-residuals.sarlm(msl)
#-Normalidad
ols_test_normality(residuals.sarlm(msl))

qqnorm(rsml, ylab="Residuos", xlab="Cuantiles te?ricos",main="",pch=20,col=rgb(0.3,0.5,1,0.4))
qqline(rsml,lwd=1.9)
boxplot(rsml, data=mapa.data,col=rgb(0.3,0.5,1,0.4),alpha=0.1,main=" ",horizontal = T, xlab="Residuos",pch=20)
hist(rsml,  main="" , breaks="Sturges" , col= rgb(0.3,0.5,1,0.4), xlab="Residuos" , ylab="Frecuencia")

#- Heterocedasticidad
bptest.sarlm(msl)
hetero.plot <- function(model) {
  plot(residuals(model) ~ fitted(model),pch=20,xlab="Predicciones", ylab="Residuales",col="blue")
  abline(h=0, lty="dotted")
  lines(lowess(fitted(model), residuals(model)), col="red")
}
hetero.plot(msl)

#-Autocorrelaci?n 
moran.test(x=rsml,orden,zero.policy =T,randomisation =T,alternative = "two.sided")




