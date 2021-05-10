#----------------------------------PROYECTO GEOESTADÍSTICA-----------------------------------------#
#---------------------------Geraldine Xiomar Buitrago Mendieta-------------------------------------#
#--------------------------------Isabella Gonzalez Mendez------------------------------------------#
#------------------------------Paula Nicole Becerra Quintero---------------------------------------#
#------------ANÁLISIS GEOESTADÍSTICO DE LA VARIABLE TEMPERATURA EN CROACIA EL 2008_02_26-----------#

#-1. Cargar Datos para ajuste de variograma

library(sp)
load("G:/Mi unidad/11. Geoestadistica/Geoestadistica/Proyecto Geoestadística/1. Datos/IDSTA.ov.rda")
load("G:/Mi unidad/11. Geoestadistica/Geoestadistica/Proyecto Geoestadística/1. Datos/croatia.grilla 25s.2008.rda")
class(IDSTA.ov)
class(IDSTA.OV)

#-2. Construcción de Base de datos a usar correspondiente a datos capturados el 26 de febrero del 2008

#--2.1. Selección de variables a usar

Temp_260208<-data.frame(IDSTA.ov$HRdem, IDSTA.ov$HRdsea, IDSTA.ov$HRtwi, IDSTA.ov$LST2008_02_26, 
                        IDSTA.ov@coords)

#--2.2. Renombrar variables

names(Temp_260208)<-c("Altura","D_costera","I_humedad", "Temperatura", "Coo_X",
                      "Coo_Y")

#--2.3. Eliminar registros con valores nulos en variable "Temperatura"

#---2.3.1. Revición de cuales valores son nulos en cada campo

sapply(Temp_260208, function(x) sum(is.na(x))) # R: 5 NA en Temperatura

#---2.3.2. Omisión en la BDD registros con valores nulos en Temperatura

Temp_260208<-na.omit(Temp_260208,cols=Temperatura)
sapply(Temp_260208, function(x) sum(is.na(x))) # R: 0 NA en Temperatura

#-3. ANÁLISIS EXPLORATORIO DE LOS DATOS

#--3.1. variable Respuesta

#---3.1.1. Análisis de atípicos

X11()
Bp<-boxplot(Temp_260208$Temperatura, main="Identificación de atípicos", ylab="Temperatura")
Bp$out  # R : El dato extremo corresponde a una temperatura de 6.51, la menor observada
Temp_260208<-Temp_260208[-c(40),] #Se elimina el dato extremo, dado que corresponde a una sola
#observación. Ello garantiza menor dispersión en los datos.
Bp<-boxplot(Temp_260208$Temperatura) #Sin precencia de atípicos

#---3.1.2. Resumen estadísticas

summary(Temp_260208$Temperatura)


#---3.1.3. Análisis de normalidad

#----3.1.3.1. Gráficas

library(ggplot2)
X11()
qplot(Temp_260208$Temperatura,geom="histogram", xlab = "Temperatura", 
      ylab = "Frecuencia",fill=I("darkviolet"),col=I("black"), alpha=I(.5))

X11()
qplot(Temp_260208$Temperatura,geom="density",col=I("black"),fill=I("darkviolet"), alpha=I(0.5),
      xlab="Temperatura",ylab="Densidad")

X11()
qqnorm(Temp_260208$Temperatura, ylab="Temperatura", xlab="Cuantiles teoricos",main="",
       pch=1,col=I("black"))
qqline(Temp_260208$Temperatura,col=I("red"))

#----3.1.3.2. Test normalidad

shapiro.test(Temp_260208$Temperatura) # R: W=0.98752, p-value=0.188, siguen una distribución normal
library(nortest)
lillie.test(Temp_260208$Temperatura) # R: D=0.036846, p-value=0.8792, siguen una distribución normal

#--3.2.Variables de posicionamiento

#---3.2.1. Medidas de tendencia central

summary(Temp_260208)

#---3.2.2. Creación de geodata para análisis gráfico vs variable Temperatura

library(geoR)
Temp_260208.geodata<-as.geodata(Temp_260208,coords.col = 5:6, data.col=4)

X11()
options(scipen = 999)
plot(Temp_260208.geodata, lowess=TRUE)

X11()
points.geodata(ylab="Coordenada Y",xlab="Coordenada X", Temp_260208.geodata, x.leg=380000, 
               y.leg=4805000, col.main=1, dig.leg = 1, pt.div="quintile",
               main="Temperatura de Croacia, Febrero 26 de 2008")

#--3.3. Covariables

#---3.3.1. Medidas de tendencia central

summary(Temp_260208)

#---3.3.2. Gráficas Vs Temperatura

x11()
par(mfrow=c(1,3))
plot(Temp_260208$Altura,Temp_260208$Temperatura,pch=20,xlab="Altura",ylab="Temperatura", 
     xlim=c(1,1453))
plot(Temp_260208$D_costera,Temp_260208$Temperatura,pch=20,xlab="Distancia costera",
     ylab="Temperatura", xlim=c(0,253)) 
plot(Temp_260208$I_humedad,Temp_260208$Temperatura,pch=20,xlab="Indice de humedad topografico",
     ylab="Temperatura", xlim=c(12,22))

#-4. Análisis de tendencia

#--4.1. Regresiones para análisis

Modelo1<-lm(Temperatura~Altura+D_costera+I_humedad+Coo_X+Coo_Y, data=Temp_260208)
summary(Modelo1) #Significativas al 5%: I_humedad; al 10%: Intercepto y D_Costera

library(MASS)
Mejor_modelo<-stepAIC(Modelo1,direction = "both")
summary(Mejor_modelo) #Significativas al 5%: Intercepto, D_Costera, I_humedad y Coo_Y
                      #R-cuadrado 0.1. Sin tendencia orden 1. 

         #Teniendo en cuenta la ausencia de tendencia,se eliminan las variables que no se van a usar
         #y se establece la posibilidad para uso de Kriging ordinario y kriging simple. Por ello, 
         #se calcula un atributo de errores en la base de datos para la construcción del 
         #semivariograma para KS

Temp_260208$Altura<-NULL
Temp_260208$I_humedad<-NULL
Temp_260208$D_costera<-NULL
mu<-mean(Temp_260208$Temperatura)
Temp_260208$E_Temperatura<-Temp_260208$Temperatura-mu
         
#-5. Análisis de Anisotropia

#--5.1. Creación de objeto tipo SpatialPointsDataFrame

library(sp)
Anis_Temp<-Temp_260208
coordinates(Anis_Temp)<-~Coo_X+Coo_Y
SpatialPoints(coordinates(Anis_Temp))
class(Anis_Temp) #class: SpatialPointsDataFrame

#--5.2. Evaluación de anisotropía 

#---5.2.1. Uso de estadístico para estimación de parámetros de anisotropía

library(intamap)
Est_anis<-estimateAnisotropy(Anis_Temp, depVar = "Temperatura")#Presencia de anisotropía.
Est_anis #Presencia de anisopia, $doRotation = TRUE

#---5.2.2. Análisis de semivariogramas direccionales

Dist_max<-sqrt((max(Temp_260208$Coo_X)-min(Temp_260208$Coo_X))^2+
                  (max(Temp_260208$Coo_Y)-min(Temp_260208$Coo_Y))^2)
Dist_var<-Dist_max/2 #R: 313810.6
VarDireccionales<-variog4(Temp_260208.geodata,option='bin',max.dist =Dist_var) 
X11()
plot(VarDireccionales, omni=TRUE, lwd=3, ylab="Semivarianza", xlab="Distancia") 

         # R: Los semivariogramas direccionales obtenidos muestran presencia de ruido. Por ello, 
         #    se considera necesaria una rotación.

#--5.3. Remoción de anisotropía en los datos 

#---5.3.1. Aplicación de rotación

Coors<-cbind(Temp_260208$Coo_X,Temp_260208$Coo_Y)
Anis_par<-c((Est_anis$direction)*(pi/180),Est_anis$ratio)
Coors_Rotadas<-data.frame(coords.aniso(Coors,Anis_par,reverse=FALSE))
Temp_260208_Rot<-data.frame(Coors_Rotadas$X1,Coors_Rotadas$X2,Temp_260208$Temperatura, 
                        Temp_260208$E_Temperatura)

#---5.3.2. Análisis de anisotropía sobre datos rotados

coordinates(Temp_260208_Rot)<-~Coors_Rotadas.X1+Coors_Rotadas.X2
class(Temp_260208_Rot)
estimateAnisotropy(Temp_260208_Rot,depVar = "Temp_260208.Temperatura") #Rotación FALSA

#--5.3.3. Creación de geodata y dataframe finales.

Temp_260208<-data.frame(Temp_260208_Rot@coords, Temp_260208_Rot$Temp_260208.Temperatura,
                        Temp_260208_Rot$Temp_260208.E_Temperatura)
names(Temp_260208)<-c("Coo_X","Coo_Y", "Temperatura", "E_Temperatura")

Temp_260208.geodata<-as.geodata(Temp_260208,coords.col = 1:2, data.col=3)

X11()
plot(Temp_260208.geodata)



#-6. Estimación del Semivariograma

#--6.1. Determinación de distancia máxima para cálculo del semivariograma

Dist_max<-sqrt((max(Temp_260208$Coo_X)-min(Temp_260208$Coo_X))^2+
                  (max(Temp_260208$Coo_Y)-min(Temp_260208$Coo_Y))^2)
Dist_var<-Dist_max/2 #R: 263839.78

#--6.2. Estimación semivariogramas experimentales

library(sgeostat)
Point<-point(Temp_260208, x = "Coo_X", y = "Coo_Y")
class(Point)  # R: "point" "data.frame"
Pair<-pair(Point,num.lags=40,maxdist=Dist_var)
library(geospt)
EstimaSV<-est.variograms(Point,Pair,"Temperatura",trim=0.1) 
X11()
par(mfrow=c(2,2))
plot(EstimaSV$bins,EstimaSV$robust,lty=1,ylim=c(0,13), col =1,pch= 16, 
     main = "SV experimental Robusto",xlab="Distancia", ylab="Semivarianza") 
plot(EstimaSV$bins,EstimaSV$med,lty=1,ylim=c(0,13), col =1,pch= 16, 
     main = "SV experimental Mediana",xlab="Distancia", ylab="Semivarianza") 
plot(EstimaSV$bins,EstimaSV$classic,lty=1,ylim=c(0,13), col =1,pch= 16,
     main = "SV experimental Clásico",xlab="Distancia", ylab="Semivarianza") 
plot(EstimaSV$bins,EstimaSV$trimmed.mean,lty=1,ylim=c(0,13), col =1,pch= 16,
     main = "SV experimental Media recortada",xlab="Distancia", ylab="Semivarianza") 

#sE OPTA POR UTILIZAR EL SEMIVARIOGRAMA CLÁSICO, DADO QUE MUESTRA MENOS DISPERSIÓN Y QUE
#INICIALMENTE FUE EXTRAÍDO EL ÚNICO DATO EXTREMO Y ELLO POSIBILITA SU USO. 

#---6.3. Cálculo de modelos Teóricos

#----6.3.1. Modelo Teórico exponencial

Exp.ml<-likfit(coords = cbind(Point$x,Point$y), data=cbind(Point$Temperatura), 
               nugget = 0.1,ini.cov.pars = c(3,40000),cov.model= "exp", fix.nug=T)

#----6.3.2. Modelo Teórico esférico

Sph.ml<-likfit(coords = cbind(Point$x,Point$y), data=cbind(Point$Temperatura),
               nugget = 0.75, ini = c(8,90000),cov.model="sph",fix.nug=T) 

#----6.3.3. Modelo Teórico Matern

Matern.ml<-likfit(coords = cbind(Point$x,Point$y), data=cbind(Point$Temperatura),
                  nugget = 0.1, ini = c(6,55000),cov.model="matern", kappa = 0.4, fix.nugget = T)

#----6.3.4. Modelo Teórico circular

Cir.ml<-likfit(coords = cbind(Point$x,Point$y), data=cbind(Point$Temperatura),
               nugget = 0.1, ini = c(5,50000),cov.model="cir",fix.nug=T) 

#----6.3.5. Modelo Teórico power exponencial

Pow.ml<-likfit(coords = cbind(Point$x,Point$y), data=cbind(Point$Temperatura), 
               nugget = 0.1, ini = c(8,90000),cov.model="powered.exponential",kappa=0.85,fix.nug=T) 

#---6.3.6. Gráfica de superposición con semivariograma teórico

x11()
plot(EstimaSV$bins,EstimaSV$classic,lty=1,ylim=c(0,12), col =1,pch= 16, 
     main = "SV experimental Clásico VS Teóricos",xlab="Distancia", ylab="Semivarianza") 
lines(Exp.ml,max.dist=Dist_var,lwd=3,col="lawngreen") 
lines(Sph.ml,max.dist=Dist_var,lwd=3,col="tan1") 
lines(Matern.ml,max.dist=Dist_var,lwd=3,col="tomato")
lines(Cir.ml,max.dist=Dist_var,lwd=3,col="mediumorchid")
lines(Pow.ml,max.dist=Dist_var,lwd=3,col="yellow") 
legend(x=150000, y= 5 ,c('Exponencial','Esférico','Matern','Circular','Pow. Exponencial'),
       col=c("lawngreen","tan1","tomato","mediumorchid", "yellow"), lty=c(1,1,1,1,1), 
       lwd =c(3,3,3,3,3))

#---6.3.1. Criterio AIC para elección de variograma teórico que mejor ajusta los datos

Exp.ml$AIC       #R : 673.2514
Sph.ml$AIC       #R : 672.2681
Matern.ml$AIC    #R : 667.0441   <--- modelo teórico con mejor ajuste
Cir.ml$AIC       #R : 690.4032
Pow.ml$AIC       #R : 674.1363



#---6.3.2. SCE de cada modelo

sum(Exp.ml$model.components[,3]^2)
sum(Sph.ml$model.components[,3]^2)
sum(Matern.ml$model.components[,3]^2)
sum(Cir.ml$model.components[,3]^2)
sum(Pow.ml$model.components[,3]^2)

#--6.4. Estimación parámetros modelo teórico

#---6.4.1. ML

Matern.ml<-likfit(coords = cbind(Point$x,Point$y), data=cbind(Point$Temperatura),
                  nugget = 0.1, ini = c(5,55000),cov.model="matern", kappa = 0.4, fix.nugget = T)
Matern.ml$parameters.summary

#----6.4.1.2. cambio de parámetros por estimados

Matern.ml<-likfit(coords = cbind(Point$x,Point$y), data=cbind(Point$Temperatura),
                  nugget = 0.1, ini = c(8,40000),cov.model="matern", kappa = 0.45, fix.nugget = T)

#--6.4.2. RML

Matern.rml<-likfit(coords = cbind(Point$x,Point$y), data=cbind(Point$Temperatura),
                   nugget = 0.1, ini = c(10,55000),cov.model="matern", 
                   fix.nugget = F, fix.kappa = F, lik.method= "RML")
Matern.rml$parameters.summary

#----6.4.1.2. cambio de parámetros por estimados

Matern.rml<-likfit(coords = cbind(Point$x,Point$y), data=cbind(Point$Temperatura),
                   nugget = 0.1, ini = c(7,50000),cov.model="matern", 
                   kappa = 0.39, fix.nugget = T, fix.kappa = T)
Matern.rml$AIC

#--6.4.3. OLS (Etimado en R)

Matern.ols<-likfit(coords = cbind(Point$x,Point$y), data=cbind(Point$Temperatura),
                   nugget = 0.005, ini = c(7,47000),
                   cov.model="matern",kappa = 0.43, fix.nugget = T, fix.kappa = T)

#--6.4.3. WLS (Estimado en R)

Matern.wls<-likfit(coords = cbind(Point$x,Point$y), data=cbind(Point$Temperatura),
                   nugget = 0.05, ini = c(7,80000),
                   cov.model="matern", kappa = 0.34, fix.nugget = T, fix.kappa = T)

#--6.5. Gráficas y valores de AIC

x11()
plot(EstimaSV$bins,EstimaSV$classic,lty=1,ylim=c(0,12), col =1,pch= 16, 
     main = "SV experimental Clásico VS Teóricos",xlab="Distancia", ylab="Semivarianza") 
lines(Matern.ml,max.dist=Dist_var,lwd=3,col="lawngreen") 
lines(Matern.rml, max.dist=Dist_var, lwd=3, col="tomato")
lines(Matern.ols,max.dist=Dist_var, lwd=3, col="tan1")
lines(Matern.wls,max.dist=Dist_var, lwd=3, col="Yellow")
legend(x=200000, y= 3 ,c('ML','RML','OLS','WLS'),
       col=c("lawngreen","tomato","tan1", "yellow"), lty=c(1,1,1,1), 
       lwd =c(3,3,3,3))


Matern.ml$AIC     # R: 667.5007
Matern.rml$AIC    # R: 664.5727
Matern.ols$AIC    # R: 669.5281
Matern.wls$AIC    # R: 664.8025

#-7. Métodos de interpolación

#--7.1. Kriging simple

#---7.1.1. Predicción en un punto

#----7.1.1.1. Definición de espacio en el que se pueden buscar predicciones

#summary(Temp_260208) 

      #Coo_X              Coo_Y            
      #Min.   :387331   Min.   :4705465     
      #Max.   :842443   Max.   :5137646     
      #Se seleccionan coordenadas (390000, 4800001) para predicción.

#----7.1.1.2. Predicción 

#library(gstat)
#mu<-mean(Temp_260208$Temperatura)
#x<-Temp_260208$Coo_X
#y<-Temp_260208$Coo_Y
#z<-Temp_260208$Temperatura
#df<-data.frame(x,y,z)
#so<-data.frame(390000, 4800001)
#names(so)<-c("x","y")
#coordinates(df)=~x+y
#coordinates(so)=~x+y
#vgm.matern<-vgm(15.62, "Mat", 55000,0)
#Temp_KS<-krige(z~1, df,so, vgm.matern, beta=mu)

#---7.1.2. Mapa de predicciones
#proj4string(IDSTA.OV) <- CRS("+proj=utm +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824")
#proj4string(df) <- CRS("+proj=utm +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824")
#Temp_KS<-krige(z~1,df , IDSTA.OV, vgm.matern,beta=mu)
#gridded(Temp_KS)<-TRUE

#x11()
#spplot(Temp_KS, "var1.pred", 
#       main="Temperatura media 29 de Marzo de 2008 \nPredicciones Kriging Simple", 
#       col.regions=bpy.colors(100), cuts=50, cex.main=0.2, xlab="Este (m)", ylab = "Norte (m)",
#       key.space=list(space="right", cex=1))

#x11()
#spplot(Temp_KS, "var1.var", main="Varianza Temperatura media 29 de Marzo de 2008 \nVarianza Kriging Simple", 
#       col.regions=bpy.colors(100),cuts=50, cex.main=0.2, scales = list(draw =T), 
#       xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))


#--7.2. Kriging ordinario

library(gstat)
x<-Temp_260208$Coo_X
y<-Temp_260208$Coo_Y
z<-Temp_260208$Temperatura
df<-data.frame(x,y,z)
coordinates(df)=~x+y
vgm.matern<-vgm(psill= 7, model = "Mat", range= 50000,nugget=0.1, kappa = 0.39)
proj4string(IDSTA.OV) <- CRS("+proj=utm +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824")
proj4string(df) <- CRS("+proj=utm +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824")
Temp_KO<-krige(z~1,df , IDSTA.OV, vgm.matern)
gridded(Temp_KO)<-TRUE

x11()
spplot(Temp_KO, "var1.pred", 
       main="Temperatura media 29 de Marzo de 2008 \nPredicciones Kriging Ordinario", 
       col.regions=bpy.colors(100), cuts=50, cex.main=0.2, xlab="Este (m)", ylab = "Norte (m)", 
       key.space=list(space="right", cex=1))

x11()
spplot(Temp_KO, "var1.var", 
       main="Varianza Temperatura media 29 de Marzo de 2008 \nPredicciones Kriging Ordinario", 
       col.regions=bpy.colors(100), cuts=50, cex.main=0.2, xlab="Este (m)", ylab = "Norte (m)", 
       key.space=list(space="right", cex=1))

#--7.3. Validación cruzada métodos de interpolación

#Valcruzc_KS<-krige.cv(z~1,df,IDSTA.OV,model=vgm.matern,beta=mu)
Valcruz_KO<-krige.cv(z~1,df,IDSTA.OV,model=vgm.matern)
Calcruz_Temp<- criterio.cv(Valcruz_KO)
Calcruz_Temp

#-8. Interpolación por Métodos deterministicos

#--8.0. Cargar datos con observación extrema

load("G:/Mi unidad/11. Geoestadistica/Geoestadistica/Proyecto Geoestadística/1. Datos/IDSTA.ov.rda")
Temp_260208<-data.frame(IDSTA.ov$HRdem, IDSTA.ov$HRdsea, IDSTA.ov$HRtwi, IDSTA.ov$LST2008_02_26, 
                        IDSTA.ov@coords)
names(Temp_260208)<-c("Altura","D_costera","I_humedad", "Temperatura", "Coo_X",
                      "Coo_Y")
Temp_260208<-na.omit(Temp_260208,cols=Temperatura)
sapply(Temp_260208, function(x) sum(is.na(x)))

#--8.1. IDW

source("G:/Mi unidad/11. Geoestadistica/Geoestadistica/Funciones/idw.cv.R")          
source("G:/Mi unidad/11. Geoestadistica/Geoestadistica/Funciones/graph.idw.R")

#---8.1.1. Optimización de p

p.optimo <- function(p, formula, locations, data, newdata, nmax, nmin, maxdist, var.reg){
   idw.pred <- as.data.frame(matrix(NA,nrow= nrow(data), ncol=4))
   colnames(idw.pred) <- c("x","y","var1.pred","var1.var")
   for(i in 1:(nrow(data))){
      idw.pred[i,] <- idw(formula, locations, da[-i,], newdata[i,], nmax, nmin, maxdist, idp=p)
   } 
   RMSPE <-  sqrt(sum((idw.pred$var1.pred-var.reg)^2)/nrow(da))
   RMSPE
}
x<-Temp_260208$Coo_X
y<-Temp_260208$Coo_Y
xy <- Temp_260208[5:6] #Se toman las coordenadas xy
z <- Temp_260208$Temperatura #Variable respuesta
da = data.frame(xy, z)

#---8.1.2. P óptimo

P <- optimize(p.optimo, c(0,10), formula=z~Coo_X+Coo_Y, locations=~Coo_X+Coo_Y, 
              data=da, newdata=da, nmax=15, nmin=15, maxdist=Inf, var.reg=z)
P
RMSPE_IDW<-P$objective

#--8.1.3. Crear borde para representación de IDW

load("G:/Mi unidad/11. Geoestadistica/Geoestadistica/Proyecto Geoestadística/1. Datos/croatia.grilla 25s.2008.rda")
library(geosptdb)
gridded(IDSTA.OV) <- TRUE
borde<-spsample(IDSTA.OV,n=600000,type="regular")
x11()
plot(borde)
ptos<-as.data.frame(borde)
names(ptos)<-c("x","y")
coordinates(ptos) = c("x", "y")
gridded(ptos)<-TRUE

#---8.1.4. Se calcula IDW con el p óptimo mínimo // Mapa de salida

pron.idw <- idw(z~ x+y, ~ Coo_X+Coo_Y, da, ptos, nmax=15, nmin=15, idp=P$minimum) 
max(pron.idw$var1.pred)
min(pron.idw$var1.pred)
l2 = list("sp.points", ptos, pch = 3, col="grey")
p1<-spplot(pron.idw, "var1.pred", main="Temperatura de Croacia Predicciones IDW", 
           col.regions=bpy.colors(100), cuts=70, cex.main=0.2, scales = list(draw =T), 
           xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
x11()
print(p1, split = c(1, 1, 1, 1), more = T)

#--8.2. Funciones de base radial

#---8.2.1. Multicuadrática (MQ)

datosdf<-as.data.frame(Temp_260208)
class(datosdf)
x<-Temp_260208$Coo_X
y<-Temp_260208$Coo_Y
z<-Temp_260208$Temperatura
df<-data.frame(x,y,z)
coordinates(df)=~x+y
proj4string(IDSTA.OV) <- CRS("+proj=utm +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824")
proj4string(df) <- CRS("+proj=utm +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824")
graph.rbf(z~x+y,df,eta.opt = TRUE,rho.opt = TRUE, func = "M",iter = 20,n.neigh = 15,eta.dmax = 3,
          rho.dmax = 3)
pred.rbf_MQ <- rbf(z~x+y,df, eta=0.00001, rho=0, newdata=IDSTA.OV,n.neigh=15, func="M")
coordinates(pred.rbf_MQ) = c("x","y")
x11()
gridded(pred.rbf_MQ) <- TRUE
spplot(pred.rbf_MQ,"var1.pred",cuts=30,col.regions=bpy.colors(100),
       main="Funcion de Base Radial Multicuadratica",scales = list(draw =T), xlab="Este (m)", 
       ylab = "Norte (m)", sp.layout=list(l2),key.space=list(space="rigth",cex=0.8))


#---8.2.2. Spline capa delgada (TPS)

datosdf<-as.data.frame(Temp_260208)
class(datosdf)
x<-Temp_260208$Coo_X
y<-Temp_260208$Coo_Y
z<-Temp_260208$Temperatura
data(IDSTA.OV)
df<-data.frame(x,y,z)
coordinates(df)=~x+y
proj4string(IDSTA.OV) <- CRS("+proj=utm +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824")
proj4string(df) <- CRS("+proj=utm +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824")
graph.rbf(z~x+y,df,eta.opt = TRUE,rho.opt = TRUE,func = "TPS",iter = 20,n.neigh = 15,eta.dmax = 3,
          rho.dmax = 3)
pred.rbf_TPS<- rbf(z~x+y,df, eta=1.649208, rho=0.6781803, newdata=IDSTA.OV,n.neigh=15, func="TPS")
coordinates(pred.rbf_TPS) = c("x","y")
x11()
gridded(pred.rbf_TPS) <- TRUE
spplot(pred.rbf_TPS,"var1.pred",cuts=20,col.regions=bpy.colors(100),main="Funcion de Base Radial TPS",
       scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", sp.layout=list(l2),
       key.space=list(space="rigth",cex=0.8))

#---8.2.3. Multicuadrática inversa (IM)

datosdf<-as.data.frame(Temp_260208)
class(datosdf)
x<-Temp_260208$Coo_X
y<-Temp_260208$Coo_Y
z<-Temp_260208$Temperatura
data(IDSTA.OV)
df<-data.frame(x,y,z)
coordinates(df)=~x+y
proj4string(IDSTA.OV) <- CRS("+proj=utm +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824")
proj4string(df) <- CRS("+proj=utm +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824")
graph.rbf(z~x+y,df,eta.opt = TRUE,rho.opt = TRUE,func = "IM",iter = 30, n.neigh = 15,eta.dmax = 3,
          rho.dmax = 3)
pred.rbf_IM <- rbf(z~x+y,df, eta=3, rho=0, newdata=IDSTA.OV,n.neigh=15, func="IM")
coordinates(pred.rbf_IM) = c("x","y")
x11()
gridded(pred.rbf_IM) <- TRUE
spplot(pred.rbf_IM,"var1.pred",cuts=30,col.regions=bpy.colors(100),
       main="Funcion de Base Radial Multicuadratica Inversa",scales = list(draw =T), 
       xlab="Este (m)",ylab = "Norte (m)",sp.layout=list(l2),key.space=list(space="rigth",cex=0.8))

#---8.2.4. Sline con tensión (ST)

datosdf<-as.data.frame(Temp_260208)
class(datosdf)
x<-Temp_260208$Coo_X
y<-Temp_260208$Coo_Y
z<-Temp_260208$Temperatura
data(IDSTA.OV)
df<-data.frame(x,y,z)
coordinates(df)=~x+y
proj4string(IDSTA.OV) <- CRS("+proj=utm +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824")
proj4string(df) <- CRS("+proj=utm +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824")
graph.rbf(z~x+y,df,eta.opt = TRUE,rho.opt = TRUE,func = "ST",iter = 40,n.neigh = 15,eta.dmax = 3,rho.dmax = 3)
pred.rbf_ST<- rbf(z~x+y,df, eta=0.007119251, rho=0.5214308, newdata=IDSTA.OV,n.neigh=15, func="ST")
coordinates(pred.rbf_ST) = c("x","y")
x11()
gridded(pred.rbf_ST) <- TRUE
spplot(pred.rbf_ST,"var1.pred",cuts=30,col.regions=bpy.colors(100),
       main="Funcion de Base Radial Spline con Tension",scales = list(draw =T), xlab="Este (m)", 
       ylab = "Norte (m)", sp.layout=list(l2),key.space=list(space="rigth",cex=0.8))

#---8.2.5. Spline completamente regularizada (CRS)

datosdf<-as.data.frame(Temp_260208)
class(datosdf)
x<-Temp_260208$Coo_X
y<-Temp_260208$Coo_Y
z<-Temp_260208$Temperatura
data(IDSTA.OV)
df<-data.frame(x,y,z)
coordinates(df)=~x+y
proj4string(IDSTA.OV) <- CRS("+proj=utm +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824")
proj4string(df) <- CRS("+proj=utm +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824")
graph.rbf(z~x+y,df,eta.opt = TRUE,rho.opt = T,func = "CRS",iter = 20,n.neigh = 15,eta.dmax = 3,rho.dmax = 3)
pred.rbf_CRS<- rbf(z~x+y,df, eta=0.05672633, rho=0.4143771, newdata=IDSTA.OV,n.neigh=15, func="CRS")
coordinates(pred.rbf_CRS) = c("x","y")
x11()
gridded(pred.rbf_CRS) <- TRUE
spplot(pred.rbf_CRS,"var1.pred",cuts=30,col.regions=bpy.colors(100),
       main="Funcion de Base Radial Spline Completamente Regularizada",scales = list(draw =T), 
       xlab="Este (m)", ylab = "Norte (m)", sp.layout=list(l2),key.space=list(space="rigth",cex=0.8))

#-9. Validación Cruzada métodos deterministicos 

V_MQ<-rbf.tcv(z~x+y,df, eta=0.00001, rho=0,n.neigh=15, func="M")
V_TPS<-rbf.tcv(z~x+y,df, eta=1.649208, rho=0.6781803,n.neigh=15, func="TPS")
V_IM<-rbf.tcv(z~x+y,df, eta=3, rho=0,n.neigh=15, func="IM")
V_ST<-rbf.tcv(z~x+y,df, eta=0.007119251, rho=0.5214308,n.neigh=15, func="ST")
V_CRS<-rbf.tcv(z~x+y,df, eta=0.05672633, rho=0.4143771,n.neigh=15, func="CRS")
criterio.cv(V_MQ)
criterio.cv(V_TPS)
criterio.cv(V_IM)
criterio.cv(V_ST)
criterio.cv(V_CRS)
#-10. Diseño de red

source("G:/Mi unidad/11. Geoestadistica/Geoestadistica/Funciones/network.design.R")
library(maptools)
poly<-readShapePoly("C:/Users/Geraldine/Desktop/HRV_adm/HRV_adm0")

x11()
red_ko_50 <- network.design(z~1, vgm.model=vgm.matern, npoints=50, boundary=poly, 
                            type="stratified")
x11()
red_ko_100<- network.design(z~1, vgm.model=vgm.matern, npoints=100, boundary=poly,  
                            type="stratified")
x11()
red_ko_150<- network.design(z~1, vgm.model=vgm.matern, npoints=150, boundary=poly,  
                            type="stratified")

Networks.P <- rbind(red_ko_50,red_ko_100,red_ko_150)
colnames(Networks.P) <- c("ASEPE")
Networks.P
