library(SpatialExtremes)
library(sp)
library(spgwr)
library(raster)
library(ggplot2)
library(sf)
library(adespatial)
library(ggspatial)
library(foreign)
library(spatial)
library(spData)
library(rgdal)
library(spatialreg)
library(spdep)
library(tmap)
library(RColorBrewer)
library(intervals)
library(classInt)
library(MASS)
library(gstat)
library(geoR)
library(sgeostat)
library(geospt)
library(scatterplot3d)
library(ggplot2)
library(car)
library(xtable)
library(stargazer)
library(RGeostats)
library(nortest)
library(intamap)
library(ggmap)
library(maptools)
library(maps)

# ---------------------------------------------------------------------------------------------------------
#                                  AN?LISIS ESTADISTICO Y GEOESTADISTICO DE LA TEMPERATURA
#                                          DE CROACIA
# ---------------------------------------------------------------------------------------------------------

# ------------------------------------------ CARGAR LA BASE DE DATOS ------------------------------------------------
load("/Geoestadística- Croacia/Geoestadística- Croacia/Datos/IDSTA.ov.rda")
load("/Geoestadística- Croacia/Geoestadística- Croacia/Datos/croatia.grilla 25s.2008.rda")

temp.geor<-as.data.frame(IDSTA.ov)#Geodata (coordenadas y datos)
datos_temp<-temp.geor[c("HRdem","HRdsea","HRtwi","Lat","Lon","LST2008_03_05")]#Lon, Lat, Coordendas geogr?ficas
names(datos_temp)<-c("HRdem","HRdsea","HRtwi","Lat","Lon","Temperatura")
dim(datos_temp)#DIMENSI?N DE LA BASE DE DATOS
names(datos_temp)#Nombre de los atributos
class(IDSTA.ov)#"SpatialPointsDataFrame"
croacia<-data.frame(datos_temp,IDSTA.ov@coords)

#Hay valores nulos en  las columnas
sapply(croacia, function(x) sum(is.na(x)))# Temperatura tiene 24 valores nulos
Bd_datostemp<-croacia[!is.na(datos_temp$Temperatura),]# Base de datos sin valores nulos
dim(Bd_datostemp)#DIMENSI?N DE LA BASE DE DATOS SIN VALORES NULOS
Bd_datostemp<-na.omit(Bd_datostemp,cols=LST2008_03_29)
class(Bd_datostemp)
temp.geodata<-as.geodata(Bd_datostemp,coords.col = 7:8,data.col =6,covar.col = 1:3)
attach(temp.geodata)


# -------------------------------------------AN?LISIS ESTAD?STICO-------------------------------------------------------
#Resumen de estad?sticas
resta<-summary(Bd_datostemp)
# Desviaci?n Estandar
varianza<-(sd(Bd_datostemp$Temperatura)^2)#sd=8.71022
names(Bd_datostemp)<-c("Altura","Dist_mar","Humedad","Lat","Lon","Temperatura","Este","Norte")
options(scipen = 999) #Desactivar notación científica, en este caso porque generaba los mapas con dicha notación
#Scatter plot Temperatura vs Altura (DEM)
X11()
ggplot(Bd_datostemp, aes(x=Temperatura, y=Altura,color=Temperatura,size=Temperatura)) +
  geom_point(alpha=0.7)
#Scatter plot Temperatura vs Humedad
X11()
ggplot(Bd_datostemp, aes(x=Temperatura, y=Humedad,color=Temperatura,size=Temperatura)) +
  geom_point(alpha=0.7)
#Scatter plot Temperatura vs Distancia a la costa 
X11()
ggplot(Bd_datostemp, aes(x=Temperatura, y=Dist_mar,color=Temperatura,size=Temperatura)) +
  geom_point(alpha=0.7)

x11()
points.geodata(ylab="Norte",xlab="Este",temp.geodata, x.leg=1, y.leg=4,  col.main=1, pt.div="quintile",main="Temperatura de Croacia para el 05 de marzo de 2008")
#Dispersogramas de los valores de temperatura y el histograma.
x11()
plot(temp.geodata) 
#Adici?n de linea de tendencia a los dispersogramas y el histograma 3D
X11()
plot(temp.geodata,lowess = TRUE) 
#Gr?fico de temperatura 3D
X11()
s3d<-scatterplot3d(main=c("Temperatura de Croacia","Marzo 05 de 2008"),Bd_datostemp$Este, Bd_datostemp$Norte, Bd_datostemp$Temperatura, angle=60,
                   col.main=2, xlab="Este", ylab="Norte", zlab="Temperatura",color="blue")
#Histograma de Temperatura
x11()
theme_set(theme_update())
hist(Bd_datostemp$Temperatura,  main="" , breaks="Sturges" , col= "cyan", xlab="Temperatura (ºC)" , ylab="Frecuencia")
x11()
qplot(Bd_datostemp$Temperatura,geom="density",col=I("cyan"),fill=I("cyan"), alpha=I(0.2), main="", xlab="Temperatura (°C)",ylab="Densidad")
x11()
qplot(Bd_datostemp$Temperatura, geom = "histogram",col=I("cYan"),fill=I("cyan"), alpha=I(0.4), main="", xlab="Temperatura (?C)",ylab="Frecuencia",bins=100)
par(mfrow=c(nrow=1, ncol=3))
#Q-Q Plot
X11()
qqnorm(Bd_datostemp$Temperatura, ylab="Temperatura (ºC)", xlab="Cuantiles teóricos",main="")
#Q-Q plot con linea media
qqline(Bd_datostemp$Temperatura)
#Box Plot
X11()
boxplot(Bd_datostemp$Temperatura, main="", notch=F, horizontal=T, xlab="Temperatura (ºC)")
points(mean(Bd_datostemp$Temperatura), y=1, pch=1, cex=2)
par(mfrow=c(nrow=1, ncol=1))

# -------------------------------------------NORMALIDAD-------------------------------------------------------
#Pruebas de normalidad de la variable temperatura
shapiro.test(Bd_datostemp$Temperatura)#p-value = 0.007724
#Transformacion BOX COX
TempBoxCox<-boxcox(Bd_datostemp$Temperatura~1)
lambdaTemp <- TempBoxCox$x[which.max(TempBoxCox$y)]
BoxTemp<- (((Bd_datostemp$Temperatura)^lambdaTemp)-1)/lambdaTemp 
#Pruebas transfomando por box cox
shapiro.test(BoxTemp)#p-value = 0.2181
#Grafica Variable de temperatura transformada por Box Cox
theme_set(theme_update())
x11()
qplot(BoxTemp,geom="histogram",main = "", xlab = "Temperatura Transformada", ylab = "Frecuencia",fill=I("blue"),col=I("black"), alpha=I(.3))
x11()
qplot(BoxTemp,geom="density",col=I("blue"),fill=I("blue"), alpha=I(0.2), main="", xlab="Temperatura Transformada",ylab="Densidad")

#Anamorfosis gaussiana
TempCro = Bd_datostemp[,c("Este","Norte","Temperatura")] #must be 3 columns: x,y,z
croacia_db = db.create(TempCro,ndim=2,autoname=F)
croacia.anam.fit = anam.fit(croacia_db,name="Temperatura",type="gaus")
croacia.anam.trans = anam.z2y(croacia_db,names="Temperatura",anam=croacia.anam.fit)
temp.anam = croacia.anam.trans@items$Gaussian.Temperatura
#Pruebas transformando anamorfosis gaussiana
shapiro.test(temp.anam)#p-value=1
#Grafica Variable de temperatura transformada
theme_set(theme_update())
x11()
qplot(temp.anam,geom="histogram",main = "", xlab = "Temperatura Transformada", ylab = "Frecuencia",fill=I("blue"),col=I("black"), alpha=I(.3))
x11()
qplot(temp.anam,geo m="density",col=I("blue"),fill=I("blue"), alpha=I(0.2), main="", xlab="Temperatura Transformada",ylab="Densidad")

#Creacion de Geodata adicionando transformada
tempnst<-data.frame(Bd_datostemp,temp.anam)
temp.geof<-as.geodata(tempnst,coords.col = 7:8,data.col =9) 
summary(temp.geof)
# -------------------------------------------ANISOTROPÍA-------------------------------------------------------
Dmax<-sqrt((max(Bd_datostemp[, 7])-min(Bd_datostemp[,7]))^2+(max(Bd_datostemp[,8])-min(Bd_datostemp[,8]))^2)
Dusar<-Dmax/2
coordinates(tempnst)=~Este+Norte
anisotropia<-estimateAnisotropy(tempnst,depVar = "temp.anam")#False: Los datos no presentan anisotropía, no se hace rotación
#Corrección de la anisotropía
#Semivariogramas direccionales

x11()
plot(variog4(temp.geof,max.dist = Dusar), xlab = "Distancia (m)", ylab = "Semivariograma estimado", legend = F)
legend(locator(1), legend = c(expression(0*degree), expression(45*degree), expression(90*degree), expression(135*degree)), col = 1:4, lty = 1:4)
title("Semivariograma en distintas direcciones")


#------------------------------------------ AN?LISIS DE TENDENCIA -----------------------------------------------------
Modelt1<-lm(temp.anam~Altura+Dist_mar+Humedad+Este+Norte,data=Bd_datostemp)
summary(Modelt1)
MM<-stepAIC(Modelt1,direction = "both")
summary(MM)
#Residuales y Test de Normalidad 
residuales<-residuals(MM)
shapiro.test(residuales)
ks.test(as.numeric(scale(sort(residuales))),pnorm)
#Creacion de Geodata adicionando transformada
tempnst<-data.frame(Bd_datostemp,temp.anam)
temp.geof<-as.geodata(tempnst,coords.col = 7:8,data.col =9) 
summary(temp.geof)
Base<-data.frame(temp.geof) 
Point <- point(Base, x = "Este", y = "Norte") 
class(Point)
Pair <- pair(Point,num.lags=50,maxdist=Dusar)
##Rezago: Usaremos 50 rezagos para la estimación del semivariograma puesto que esta cantidad es suficiente para
##definir claramente los parámetros de pepita, meseta y rango.
EstimaSV<- est.variograms(Point,Pair,"data",trim=0.1) 
#ModeloLST experimental Cl?sico
x11()
plot(EstimaSV$bins,EstimaSV$classic,lty=1,ylim=c(0,3), col =1,pch= 16,xlab="Distancia", ylab="Semivarianza",main="Modelo experimental clásico") 
#ModeloLST experimental mediana
x11()
plot(EstimaSV$bins,EstimaSV$med,lty=1,ylim=c(0,3), col =1,pch= 16,xlab="Distancia", ylab="Semivarianza") 
#ModeloLST experimental media recortada
x11()
plot(EstimaSV$bins,EstimaSV$trimmed.mean,lty=1,ylim=c(0,3), col =1,pch= 16,xlab="Distancia", ylab="Semivarianza") 
#ModeloLST experimental robusto
x11()
plot(EstimaSV$bins,EstimaSV$robust,lty=1,ylim=c(0,3), col =1,pch= 16,xlab="Distancia", ylab="Semivarianza") 


x11()
plot(EstimaSV$bins,EstimaSV$classic,lty=1,ylim=c(0,3), col =1,pch= 16,xlab="Distancia", ylab="Semivarianza",main="Modelo experimental clásico") 

Exp.ml<-likfit(geodata = temp.geof,nugget = 0.15,ini = c(1.7,125000),fix.nug=T) 
Sph.ml<-likfit(geodata = temp.geof,nugget = 0.15,ini = c(1.8,100000),cov.model="sph",fix.nug=T) 
Matern.ml<-likfit(geodata = temp.geof,nugget = 0.15,ini = c(1.7,125000),cov.model="mat",kappa=0.5,fix.nug=T) 
Cir.ml<-likfit(geodata = temp.geof,nugget = 0.15,ini = c(1.7,125000),cov.model="cir",fix.nug=T) 
Pow.ml<-likfit(geodata = temp.geof, nugget = 0.15,ini = c(1.7,125000),cov.model="powered.exponential",kappa=1.75,fix.nug=T) 
lines(Exp.ml,max.dist=Dusar,lwd=2,col="lawngreen") 
lines(Sph.ml,max.dist=Dusar,lwd=3,col="tan1") 
lines(Matern.ml,max.dist=Dusar,lwd=3,col="tomato") 
lines(Cir.ml,max.dist=Dusar,lwd=3,col="mediumorchid") 
lines(Pow.ml,max.dist=Dusar,lwd=3,col="yellow") 
legend(locator(1),c('Exponencial','Esferico','Matern','Circular','Pow.Exponencial'),col=c("lawngreen","tan1","tomato","mediumorchid", "yellow"), lty=c(1,1,1,1,1,1,1))
Exp.ml$AIC
Sph.ml$AIC#Modelo con mejor ajuste al modelo clásico AIC=377.1541
Matern.ml$AIC#
Cir.ml$AIC
Pow.ml$AIC


###Ajuste del variograma- Se selecciona el modelo clásico con el modelo esférico debido a que en este 
###modelo se presenta el mejor ajuste a los datos y el menor valor para el Criterio de Akaike
variogramaclasico<-variog(temp.geof,trend="cte",max.dist=Dusar, option = "cloud")
x11()
plot(EstimaSV$bins,EstimaSV$classic,lty=1,ylim=c(0,3), col =1,pch= 16,main = "ModeloLST experimental Clasico",xlab="Distancia", ylab="Semivarianza") 
sph.ml<-likfit(geodata = temp.geof,nugget = 0.15, ini = c(1.8,100000),cov.model="sph",fix.nug=T) 
lines(sph.ml,max.dist=Dusar,lwd=3,col="blue") 
sph.rml<-likfit(geodata = temp.geof,nugget = 0.15,ini = c(1.8,100000),fix.nugget= T,method='RML',cov.model = "sph")
lines(sph.rml,max.dist=Dusar,lwd=3,lty=1, col="green3")
sph.wls<-variofit(vario = variogramaclasico, nugget = 0.15,ini = c(1.8,100000),cov.model = "spheric",fix.nugget= T,weights="npairs")
lines(sph.wls,max.dist=Dusar,lwd=2,lty=1,col="coral")
legend(locator(1),c('ML','RML','WLS'),lty=c(1,1,1,1),col=c("blue","green3","coral"))
AIC(sph.ml)
AIC(sph.rml)#AIC=372.8422


#De acuerdo al AIC, el mejor ajuste está dado por Mínimos Cuadrados Restringidos
sph.ml
#Pepita=0.1758
#Meseta=1.7564
#Rango=100000

#---------------------------------------------Predicción espacial Kriging Simple------------------------------------------
#Predicción en un punto del espacio
mu<-mean(temp.anam)

#Norte min = 4705463;Norte Max=5137646
#Este min=387331; Este max= 842443
#Altura min=2m; Altura max=1453 m
#Estos datos pemiten definir el espacio donde podemos buscar predicciones puntuales
library(geospt)
library(gstat)
library(sp)
x<-tempnst$Este
y<-tempnst$Norte
h<-tempnst$Altura
z<-tempnst$temp.anam

xy<-tempnst[,7:8]

pts = data.frame(xy, z=z)
coordinates(pts) <- c("Este", "Norte")
proj4string(pts)<-CRS("+init=epsg:3116")


vgm.sph<-vgm(psill=1.8, model="Sph", range=100000,nugget=0.1)
#Mapa de predicción
detach("package:sgeostat")
proj4string(IDSTA.OV) <- CRS("+proj=utm +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824")
proj4string(pts) <- CRS("+proj=utm +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824")

gridded(IDSTA.OV)<-TRUE
temp.ks<-krige(z~1,pts,IDSTA.OV,vgm.sph,beta=mu)
#Transformación inversa de la anamorfosis para generación de mapa de predicción
Tempks.inv.anam = cbind(coordinates(temp.ks),temp.ks$var1.pred)
Tempks.inv.anam.db = db.create(Tempks.inv.anam,autoname = F)
tempks_pred = anam.y2z(Tempks.inv.anam.db,names="V3",anam = croacia.anam.fit )


temp.ks$pred_temp <- tempks_pred@items$Raw.V3
temp.ks$var_pred<-temp.ks$var1.var
#Graficación de mapas
x11()
spplot(temp.ks, "pred_temp", main="Temperatura media 05 de Marzo de 2008 \nPredicciones Kriging Simple", col.regions=bpy.colors(100), 
               cuts=50, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=1))

x11()

spplot(temp.ks, "var1.var", main="Temperatura media 05 de Marzo de 2008 \nVarianza Kriging Simple", col.regions=bpy.colors(100), 
       cuts=50, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
detach("package:sgeostat")
#----------------------------------------Predicción espacial Kriging Ordinario--------------------------------------------
temp.ko<-krige(z~1,pts , IDSTA.OV,vgm.sph)
Tempko.inv.anam = cbind(coordinates(temp.ko),temp.ko$var1.pred)
Tempko.inv.anam.db = db.create(Tempko.inv.anam,autoname = F)
tempko_pred = anam.y2z(Tempko.inv.anam.db,names="V3",anam = croacia.anam.fit )
temp.ko$pred_temp <- tempko_pred@items$Raw.V3
temp.ko$var_pred<-temp.ko$var1.var
x11()
spplot(temp.ko, "pred_temp", main="Temperatura media 05 de Marzo de 2008 \nPredicciones Kriging Simple", col.regions=bpy.colors(100), 
       cuts=50, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=1))
x11()
spplot(temp.ko, "var1.var", main="Temperatura media 05 de Marzo de 2008 \nPredicciones Kriging Simple", col.regions=bpy.colors(100), 
       cuts=50, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=1))
#-------------------------------------Validación cruzada-------------------------------------------------------------------
temp.ks.cv<-krige.cv(z~1,pts,IDSTA.OV,model=vgm.sph,beta=mu)
temp.ko.cv<-krige.cv(z~1,pts,IDSTA.OV,model=vgm.sph)
resultados.cv.z <- rbind(criterio.cv(temp.ks.cv), criterio.cv(temp.ko.cv))
rownames(resultados.cv.z) <- c("Kriging simple", "Kriging ordinario")
resultados.cv.z
xtable(resultados.cv.z)



##########################################
#### PREDICCION COKRIGING ################

################################################
#####MATRIZ DE CORRELACION #####################
Correlacion<-cor(tempnst,method="pearson")
round(Correlacion,digits = 2)

round(cor(var.covar),2)
library(corrplot)
x11()
corrplot(Correlacion,type="upper",method="shade",
         tl.col = "black", tl.srt = 45,addCoef.col = "Black", diag=F, addshade = "all") 

#Las mayores correlaciones están dadas con los SDT (Sólidos Disueltos Totales), dureza,
#temperatura y pH
####ANAMORFOSIS GAUSSIANA PARA COVARIABLES

##Altura
#Pruebas de normalidad
x11()
qplot(tempnst$Altura,geom="density",col=I("red"),fill=I("red"), alpha=I(0.2), 
      main="", xlab="Altura",ylab="Densidad")
shapiro.test(tempnst$Altura)#p-value = 0.00000000000001439. No hay normalidad
#Transformación
Alt.sel = tempnst[,c("Este","Norte","Altura")]
Alt_db = db.create(Alt.sel,ndim=2,autoname=F)
Alt.anam.fit = anam.fit(Alt_db,name="Altura",type="gaus")
Alt.anam.trans = anam.z2y(Alt_db,names="Altura",anam=Alt.anam.fit)
Alt.anam = Alt.anam.trans@items$Gaussian.Altura
#Pruebas de normalidad para la transformación
shapiro.test(Alt.anam)#p-value = 0.0204 

x11()
qplot(Alt.anam,geom="density",col=I("red"),fill=I("red"), alpha=I(0.2), main="", 
      xlab="Dureza transformada con anamorfosis gaussiana",ylab="Densidad")
#Transformacion BOX COX
AltBoxCox<-boxcox(tempnst$Altura~1)
lambdaAlt <- AltBoxCox$x[which.max(AltBoxCox$y)]
BoxAlt<- (((tempnst$Altura)^lambdaAlt)-1)/lambdaAlt
shapiro.test(BoxAlt)#p-value = 0.0055

##Humedad
#Pruebas de normalidad
x11()
qplot(tempnst$Humedad,geom="density",col=I("red"),fill=I("red"), alpha=I(0.2), 
      main="", xlab="Humedad",ylab="Densidad")
shapiro.test(tempnst$Humedad)#P-valor=0.0009631. No hay normalidad
#Transformación
Hum.sel = tempnst[,c("Este","Norte","Humedad")]
Hum_db = db.create(Hum.sel,ndim=2,autoname=F)
Hum.anam.fit = anam.fit(Hum_db,name="Humedad",type="gaus")
Hum.anam.trans = anam.z2y(Hum_db,names="Humedad",anam=Hum.anam.fit)
Hum.anam = Hum.anam.trans@items$Gaussian.Humedad
#Pruebas de normalidad para la transformación
shapiro.test(Hum.anam)#p-value = 0.9922 - Normalizado

#Transformacion BOX COX
HumBoxCox<-boxcox(tempnst$Humedad~1)
lambdaHum <- HumBoxCox$x[which.max(HumBoxCox$y)]
BoxHum<- (((tempnst$Humedad)^lambdaHum)-1)/lambdaHum
shapiro.test(BoxHum)#p-value = 0.0055

x11()
qplot(Hum.anam,geom="density",col=I("red"),fill=I("red"), alpha=I(0.2), main="", 
      xlab="Temperatura transformada con anamorfosis gaussiana",ylab="Densidad")

##Distancia a costa
#Pruebas de normalidad

tempnst$Dist_mar[tempnst$Dist_mar==0.0] <- 0.1

x11()
qplot(tempnst$Dist_mar,geom="density",col=I("red"),fill=I("red"), alpha=I(0.2), 
      main="", xlab="Distancia a costa",ylab="Densidad")
shapiro.test(tempnst$Dist_mar)#P-valor=1.05e-10 No hay normalidad
#Transformación

Dist.sel = tempnst[,c("Este","Norte","Dist_mar")]
Dist_db = db.create(Dist.sel,ndim=2,autoname=F)
Dist.anam.fit = anam.fit(Dist_db,name="Dist_mar",type="gaus")
Dist.anam.trans = anam.z2y(Dist_db,names="Dist_mar",anam=Dist.anam.fit)
Dist.anam = Dist.anam.trans@items$Gaussian.Dist_mar
#Pruebas de normalidad para la transformación
shapiro.test(Dist.anam)#p-value = 5.622e-10 - No Normalizado
theme_set(theme_update())
x11()
qplot(Dist.anam,geom="density",col=I("red"),fill=I("red"), alpha=I(0.2), main="",
      xlab="SDT transformada por anamorfosis gaussiana",ylab="Densidad")

#Transformacion BOX COX
DistBoxCox<-boxcox(tempnst$Dist_mar~1)
lambdaDist <- DistBoxCox$x[which.max(DistBoxCox$y)]
BoxDist<- (((tempnst$Dist_mar)^lambdaDist)-1)/lambdaDist
shapiro.test(BoxDist)#p-value = 8.378e-13

#Debido a que no hubo buena normalización de la distancia a costa y la altura, se procede a trabajar con humedad como covariable

Dframe.cok<-data.frame(tempnst$Este,tempnst$Norte,tempnst$temp.anam,Hum.anam)
names(Dframe.cok)<-c("Este","Norte","Temperatura","Humedad")
coordinates(Dframe.cok)<-~Este+Norte
library(sp)
proj4string(Dframe.cok)<-CRS("+init=epsg:3116")
library(gstat)
g2<-gstat(id="Temperatura",formula=Temperatura~1,data=Dframe.cok)
g2<-gstat(g2,id="Humedad",formula=Humedad~1,data=Dframe.cok)
g2<-gstat(g2,id="Temperatura",model=vgm.sph)
g2<-gstat(g2,id="Humedad",model=vgm.sph)
vcross<-variogram(g2,cutoff=(Dusar))
g2<-gstat(g2,id="Temperatura",model=vgm.sph,fill.all=T)
#Ajuste del modelo
fit.g2<-fit.lmc(vcross,g2,correct.diagonal=1.01)
x11()
plot(variogram(g2), model=fit.g2$model, col="red",pch=20)
#Validación cruzada
cv<-gstat.cv(fit.g2)
criterio.cv(cv)

######GRILLA PARA LA PREDICCION#######
proj4string(IDSTA.OV) <- CRS("+proj=tmerc +lat_0=4.59620041666667 +lon_0=-74.0775079166667 +k=1 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +units=m +no_defs")

Temp.cok <- predict(fit.g2, newdata = IDSTA.OV)

#Transformación inversa de anamorfosis gaussiana para co-kriging ordinario

Tempcok.inv.anam = cbind(coordinates(Temp.cok),Temp.cok$Temperatura.pred)
Tempcok.inv.anam.db = db.create(Tempcok.inv.anam,autoname = T)
Tempcok_pred = anam.y2z(Tempcok.inv.anam.db,names="V3",anam = croacia.anam.fit )

Temp.cok$pred_Temp <- Tempcok_pred@items$Raw.V3
Temp.cok$var_pred<-Temp.cok$var1.var

pred_cok<-st_as_sf(Temp.cok)
summary(pred_cok)
x11()
spplot(Temp.cok, "pred_Temp", main="Temperatura media 05 de Marzo de 2008 \nPredicciones Kriging Simple", col.regions=bpy.colors(100), 
       cuts=50, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=1))

x11()
spplot(Temp.cok, "Temperatura.var", main="Temperatura media 05 de Marzo de 2008 \nPredicciones Kriging Simple", col.regions=bpy.colors(100), 
       cuts=50, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=1))

