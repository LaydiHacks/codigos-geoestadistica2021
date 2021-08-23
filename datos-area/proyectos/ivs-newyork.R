#----------------------------------PROYECTO ECONOMETRIA ESPACIAL---------------------------------#
#---------------------------Viviana Bautista - Alejandro Dimate Rodriguez------------------------#
#------------AN?LISIS GEOESTAD?STICO DE INDICE DE VENTAJA SOCIECONOMICA NEW YORK 2010-2014--------#

url_local = "GEOESTADISTICA/codigos-geoestadistica2021/"

#Instalar Paquetes
install.packages("Rcpp")
install.packages("spgwr")
install.packages("raster")
install.packages("tmap")
install.packages("lmtest")
install.packages("RColorBrewer")
install.packages("classInt")
install.packages("spdep")
install.packages("sphet")
install.packages("pgirmess")
install.packages("spatialreg")
install.packages("dbscan")
install.packages("adespatial")
install.packages("rgdal")
#Cargar Paquetes y Funciones
library(sp)
library(Rcpp) 
library(spData)
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
library(RGeostats) #Instalar poe Zip
source(paste(url_local,"funciones/northarrow.r", sep=""))
source(paste(url_local,"funciones/scalebar.R", sep=""))
source(paste(url_local,"funciones/moran.R", sep=""))
source(paste(url_local,"funciones/geary.R", sep=""))
source(paste(url_local,"funciones/moranbi.test.R", sep=""))
source(paste(url_local,"funciones/randomize_vector.R", sep=""))
source(paste(url_local,"funciones/moran.cluster.R", sep=""))
source(paste(url_local,"funciones/moran.bi.R", sep=""))
source(paste(url_local,"funciones/moran.cluster.R", sep=""))
source(paste(url_local,"funciones/getis.cluster.R", sep=""))
source(paste(url_local,"funciones/localmoran.bi.R", sep=""))
source(paste(url_local,"funciones/moranbi.plot.R", sep="")) 
source(paste(url_local,"funciones/quantile.e.R", sep=""))
source(paste(url_local,"funciones/sp.na.omit.R", sep=""))
source(paste(url_local,"funciones/correlogram.d.R", sep=""))
source(paste(url_local,"funciones/sp.correlogram.R", sep=""))
source(paste(url_local,"funciones/spcorrelogram.bi.R", sep=""))
source(paste(url_local,"funciones/moran.bi1.R", sep=""))
source(paste(url_local,"funciones/moranbi1.test.R", sep=""))
source(paste(url_local,"funciones/geary.bi.R", sep=""))
source(paste(url_local,"funciones/test.W.R", sep=""))

shp_url = "db/ventajasocioeconomica/condados/newyork.shp" 
shp <-paste(url_local,shp_url, sep="")
mapa <- readOGR(shp)



###############################################
#######      VARIABLES        #################
###############################################


IVS <- as.data.frame(mapa)$IVS_TP
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
data  <- data.frame(mapa$IVS_TP, X,Y)


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
#######     Preparación de Datos      #########
###############################################

#Se realizó una traslación de datos al cuadrante positivo
mapa$IVS_TP <- mapa$IVS + 10
boxplot(mapa$IVS_TP,col="seagreen3",ylab="TEXTODELEJE")
hist(mapa$IVS_TP)
summary(mapa$IVS_TP)
shapiro.test(mapa$IVS_TP)
qqnorm(mapa@data$IVS_TP)
qqline(mapa@data$IVS_TP)
IVS <- as.data.frame(mapa)$IVS_TP
mapa$IVS <- mapa$IVS_TP




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
##########################################
############### MODELOS  ##############
##########################################
#=============================================================

mapa.data<-as.data.frame(mapa)
orden<-nb2listw(gabrielnb, style="W", zero.policy =T)
WX <- lag.listw(orden, X)
WY <- lag.listw(orden, Y)

#=============================================================
############### MODELO CLASICO - LINEAL  ##############
#=============================================================

formula_modelo = IVS~AGE65+AH+EN5+log(IPC)+MPAREN+NDS25+NOSEGURO+PMIN+POBREZA+DISC+AGE17+CIVILES16
mclasico<-lm(formula_modelo,data = mapa.data)
summary(mclasico)
#Sacamos las variables no significativas
formula_mc <- IVS~AGE65+AH+EN5+log(IPC)+MPAREN+NDS25+NOSEGURO+PMIN+POBREZA
mejorclasico<-lm(formula_mc, data=mapa.data)
summary(mejorclasico)

# Validación supuestos
bptest(mejorclasico)
resettest(mejorclasico)
raintest(mejorclasico)
shapiro.test(residuals(mejorclasico))
library(car)
vif(mejorclasico)

BCresiduos <- mapa
residuos$mc <-residuals(mejorclasico)
spplot(residuos, "mc", col.regions = rev(terrain.colors(20))) #Mapa de Residuos
shapiro.test(residuos$mc)  
bptest(mejorclasico)
lm.morantest(mejorclasico,orden) 
#Corremos todas las pruebas
pruebas <- lm.LMtests(mejorclasico, listw = orden, test = "all")
summary(pruebas)

#=============================================================
############### MODELO SPATIAL LAG  ##############
#=============================================================
formula_mc <- IVS~AGE65+AH+EN5+log(IPC)+MPAREN+NDS25+NOSEGURO+PMIN+POBREZA
msl<-lagsarlm(formula_modelo,data=mapa,listw = orden,method="eigen")
summary(msl, Nagelkerke=T, correlation=TRUE)
rsml<-residuals.sarlm(msl)
moran.test(x=rsml,orden,zero.policy =T,randomisation =T,alternative = "two.sided")

#MODELO GNSS (Da Lambda = 0, probamos con Burbin Spatial)
mod.GNS <- sacsarlm(formula_modelo, listw = orden, data = mapa.data, type = "sacmixed")
summary(mod.GNS,Nagelkerke=T)

#SDM - MODELO DURBIN SPATIAL 
  lagsd<-lagsarlm(formula_modelo, listw=orden,data=mapa.data, type="Durbin") 
  summary(lagsd, correlation=TRUE)
  #Obtener las variables rezagadas
  W_DISC <- lag.listw(orden, mapa.data$DISC)
  W_AGE65 <- lag.listw(orden, mapa.data$AGE65)
  W_RENT30 <- lag.listw(orden, mapa.data$RENT30)
  W_MPAREN <- lag.listw(orden, mapa.data$MPAREN)
  LOGIPC =log(IPC)
  W_IPC <- lag.listw(orden, LOGIPC)
  W_POBREZA <- lag.listw(orden, mapa.data$POBREZA)
  nfm_SDM = IVS_T~AGE65+AH+EN5+log(IPC)+MPAREN+NDS25+NOSEGURO+PMIN+POBREZA+RENT30
  #Nuevo modelo
  lagsd<-lagsarlm(nfm_SDM, listw=orden,data=mapa.data, type="Durbin")
  summary(lagsd, correlation=TRUE)
  #Otra opciÃ³n del surbin spatial
  mod.SDEM1 <- errorsarlm(l_salario ~ desempleo, listw = W1_list, data = spatial_data, etype = "emixed")
  summary(mod.SDEM1,Nagelkerke=T)
  #Una forma de limpiar este modelo es optienedo las variables rezagadas e incluirlas en el modelo
  #col.poly1$WX <- lag.listw(a.lw, col.poly1$X) #Varibale rezagada
  col.lagx <- lagsarlm(CRIME ~ INC + HOVAL + X + WX, data=mapa.data, listw=orden) 
  
  




## SARAR - AUTOREGRESIVO EN EL ERROR
  m.sarar <- sacsarlm(formula_modelo, data=mapa.data, listw=orden, method="eigen")    
  summary(m.sarar, correlation=TRUE, Nagelkerke=T)
  residuals.sarlm(nepal.sarar)
  coef.sarlm(nepal.sarar)
  bptest.sarlm(nepal.sarar)
  shapiro.test(nepal.sarar$residuals)
  ols_test_normality(nepal.sarar$residuals)
  moran.test(x=nepal.sarar$residuals,orden,zero.policy =T,randomisation =T,alternative = "two.sided")
  
###GWR
  FM_XY = IVS~AGE17+AGE65+AH+CIVILES16+DISC+EN5+log(IPC)+MPAREN+NDS25+NOSEGURO+PMIN+POBREZA+RENT30+X+Y+WX+WY
  lm.global <- lm(FM_XY, data=mapa.data)
  summary(lm.global)
  N_FM_XY  = IVS~AGE65+AH+EN5+log(IPC)+MPAREN+NDS25+NOSEGURO+PMIN+POBREZA+RENT30+X
  lm.global <- lm(N_FM_XY, data=mapa.data)
  summary(lm.global)
  
  #GW Model
  # Crossvalidation of bandwidth for geographically weighted regression
  adapt <- gwr.sel(FM_XY, data=mapa.data, coords=cbind(X,Y)) 
  #Define Variables DeV: Var Dep  e  InV: Var Indep
  DeV<-"IVS"
  IPC = log(IPC)
  InV<-c("AGE17","AGE65","PMIN","NOSEGURO","NDS25","IPC","MPAREN","CIVILES16","DISC","EN5","X","Y","WX","WY")

  model.sel <- model.selection.gwr(DeV,InV, data =mapa, kernel = "bisquare", adaptive = TRUE, bw = 10)
  sorted.models <- model.sort.gwr(model.sel, numVars = length(InV), ruler.vector = model.sel[[2]][,2])
  model.list <- sorted.models[[1]]
  x11()
  
  str(model.view.gwr(DeV, InV, model.list = model.list))
  model.view.gwr(DeV, InV, model.list = model.list)
  mejor_bw_gwr <- bw.gwr(FM_XY, data = mapa, approach = "AICc", kernel = "bisquare", adaptive = TRUE)
  gwr.result<-gwr.basic(FM_XY, data=mapa, kernel="bisquare", adaptive=TRUE, bw=mejor_bw_gwr)
    
#=============================================================
############### AJUSTE DEL MODELO  ##############
#=============================================================  

  
  # primero se seleccionan las variables que son utiles para el modelo:
  datos.col <- as.data.frame(col.poly)
  d.col <- as.data.frame(cbind(datos.col$POLYID, datos.col$HOVAL, datos.col$INC, datos.col$CRIME))
  colnames(d.col) <- c("POLYID", "HOVAL", "INC", "CRIME")
  attach(d.col)
  
  # La variable dependiente de nuestro modelo serï¿½ la variable CRIME, y las variables 
  # explicativas son INC y HOVAL. Se procede a separarlas:
  Y <- CRIME
  INC <- datos.col$INC
  HOVAL <- datos.col$HOVAL
  X <- cbind(INC, HOVAL)
  
  # Adicionamos una columna de unos (1) al conjunto de las variables independientes en el objeto x 
  # previamente creado. Con base en los objetos x y y se harï¿½ posteriormente la estimaciï¿½n del 
  # parï¿½metro (rho) del modelo. Observamos las primeras 5 filas del objeto x:
  
  X <- cbind(1, X)
  
  # Modelo Spatial Lag:                                                                   col.k5
  col.lag.sm <- lagsarlm(CRIME ~ INC + HOVAL + X + Y, data=as.data.frame(col.poly), listw=col.k5) # CRIME ~1
  summary(col.lag.sm, Nagelkerke=T)
  predict.sarlm(col.lag.sm)
  AIC(col.lag.sm)
  deviance.sarlm(col.lag.sm)
  residuals.sarlm(col.lag.sm)
  coef.sarlm(col.lag.sm)
  fitted.sarlm(col.lag.sm)
  bptest.sarlm(col.lag.sm)
  hetero.plot <- function(model) {
    plot(residuals(model) ~ fitted(model))
    abline(h=0, lty="dotted")
    lines(lowess(fitted(model), residuals(model)), col="red")
  }
  hetero.plot(col.lag.sm)
  



  # Transformación Box-Cox
  BC <-  boxCox(formula_mc, data = mapa.data, lambda = seq(-3, 3, len = 20))
  lambda <- BC$x[which.max(bc$y)]
  IVS_BC = (mapa$IVS^lambda - 1)/lambda
  formula_mc <- IVS_BC ~ AH+EN5+log(IPC)+MPAREN+NDS25
  mejorclasico<-lm(formula_mc)
  summary(mejorclasico)
  
 
  
  #######     Amorfosis Gaussiana       #########
  IVSselect = data[,c("X","Y","mapa.IVS")] 
  IVS.rgdb <- db.create(IVSselect,ndim=2,autoname=F)
  IVS.herm <- anam.fit(IVS.rgdb,name="mapa.IVS",type="gaus")
  IVS.hermtrans <- anam.z2y(IVS.rgdb,names="mapa.IVS",anam=IVS.herm)
  IVS_T <- IVS.hermtrans@items$Gaussian.mapa.IVS
  shapiro.test(IVS_T) #Normalizo

