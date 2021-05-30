                            #################################################
                            ########## Ejercicio Precipitaci?n  #############
                            #################################################

#Actividad: 
#Calcular el Semivariograma Experimental mediante los cuatro metodos: cl?sico, robusto, mediana, media recortada
#Ajustar el semivariograma Experimental en cada uno de los metodos mediante MCO, MCP, MV y MVR
#Ajustar los semivariogramas Teoricos Esferico, Exponencial, Gaussiano y Matern (kappa diferente de 0.5) al semivariograma experimental

 
#Iniciamos las librerias necesarias                                                   
library(gstat)
library(geoR)
library(sgeostat)
library(geospt)
#Cargamos lo datos "PRECIPITACION"
data(preci)

#Accedemos directamente a la base de datos
attach(preci)


#Continuando con el ajuste del Semivariograma

#A continuaci?n se calculan los semivariogramas experimentales: Cl?sico, Robusto, Mediana y Media recortada
#Teniendo en cuenta que para ello se utilizan 7 rezagos 
midataframe<-data.frame(preci[,2:4])
names(midataframe) <- c("x", "y", "p")
P.L.point <- point(midataframe) #Convertir en clase point
P.L.pair <- pair(P.L.point,num.lags=7,maxdist=7) #Convertirlo en clase parejas, el maxdist la idea es que ningun intervalo quede vacio para evitar errores

P.L.v <- est.variograms(P.L.point,P.L.pair,"p",trim=0.1)                               #El parametro trim se utiliza para calcular la media  recortada  
#Lo que tenemos es, la primer columna con el rezago que para el caso tenemos 7 rezagos y el centroide es bins
#Luego tenemo los semivariagramas clasico, robusto, median y el de la media recortada con el numero de parejas


#A continuacion se grafica los semivariogramas experimentales
plot(P.L.v$bins,P.L.v$robust,lty=1, col =1,main = "Ajuste de Modelos de Semivarianza",xlab="Distancia", ylab="Semivarianza", type="l")
lines(P.L.v$bins,P.L.v$med, col=2)
lines(P.L.v$bins,P.L.v$classic, col=3)
lines(P.L.v$bins,P.L.v$med.trim, col=4)
legend(locator(1), c("Clásico", "Robusto", "Mediana", "Media Recortada"), col=c(1,2,3,4), lty=c(1,1,1,1))
detach("package:sgeostat")      # Se inactiva sgeostat, dado que genera conflicto con gstat en algunas funciones

#A continuacion se ajustan algunos semivariogramas teoricos al semivariograma experimental  

dir.hor <- seq(0, 0, length.out=7) 
dir.ver <- seq(0, 0, length.out=7)
id <- seq (length.out=7) 
id <- rep("var1",7)

##ESFERICO CLASICO
y <- data.frame (P.L.v$n, P.L.v$bins,P.L.v$classic,dir.hor,dir.ver,id)
names(y) <- c("np", "dist", "gamma", "dir.hor","dir.ver","id")
class(y) <- c("variogram","gstatVariogram","data.frame")
preci.geoR <- as.geodata(preci, coords.col = 2:3, data.col = 4)    # Objeto del tipo geodata (coordenadas y datos)
plot(preci.geoR)

Sph.ml <-likfit(preci.geoR, trend = "1st", ini = c(250, 3), fix.nugget = F, cov.model="sph", lik.method = "ML")        # metodo MV
Sph.reml <-likfit(preci.geoR, trend = "1st", ini = c(250, 3), fix.nugget = F, cov.model="sph", lik.method = "REML")    # metodo MVR
# Sph.reml <- fit.variogram.reml (prec~1, ~x+y, preci, model = vgm(200, "Sph",4.2, 0))   # metodo MVR
Sph.ols <- fit.variogram(y, vgm(250, "Sph", 3, 5),fit.method = 6) #El metodo de ajusto se refiere al metodo para definir los pesos Wi                                               # metodo 6 MCO
Sph.wls <- fit.variogram(y, vgm(250, "Sph", 3, 5),fit.method = 7)                                                      # metodo 7 MCP
print(list(Sph.ml,Sph.reml,Sph.ols, Sph.wls))

dist.s <- P.L.v$bins       
#El variogramLine evalua el variograma y esto permite generar un semivarianza teorica
Sph.ML <- variogramLine(vgm(Sph.ml$cov.pars[1], "Sph", Sph.ml$phi,Sph.ml$nugget), min=0, dist_vector=dist.s)
Sph.RML <- variogramLine(vgm(Sph.reml$cov.pars[1], "Sph", Sph.ml$phi,Sph.ml$nugget), min=0, dist_vector=dist.s)
Sph.WLS <- variogramLine(vgm(Sph.wls$psill[2], "Sph", Sph.wls$range[2],Sph.wls$psill[1]), min=0, dist_vector=dist.s)
Sph.OLS <- variogramLine(vgm(Sph.ols$psill[2], "Sph", Sph.wls$range[2],Sph.ols$psill[1]), min=0, dist_vector=dist.s)
#Resta del semivariograma exprimental con el teorico
CME.Sph.ML <- sum((P.L.v$classic-Sph.ML$gamma)^2)/7
CME.Sph.RML <- sum((P.L.v$classic-Sph.RML$gamma)^2)/7
CME.Sph.OLS <- sum((P.L.v$classic-Sph.OLS$gamma)^2)/7
CME.Sph.WLS <- sum((P.L.v$classic-Sph.WLS$gamma)^2)/7
print(data.frame(CME.Sph.ML,CME.Sph.RML,CME.Sph.OLS,CME.Sph.WLS)) 
#Tenemos un numero de comparación por MCO, esta comparación es comparable porque lleva todo a un mismo marco

plot(P.L.v$bins,P.L.v$classic,lty=2,pch=1,lwd=2, bg="yellow",type = "p", ylim=c(0,350), col =1,font.main=3,main = ("AJUSTE DE MODELO ESFERICO CLASICO"),xlab="Distancia", ylab="Semivarianza")
lines ( Sph.ML, col =2,lty=6,lwd=2)
lines ( Sph.RML, col =3,lty=6,lwd=2)
lines ( Sph.WLS, col =4,lty=6,lwd=2)
lines ( Sph.OLS, col =5,lty=6,lwd=2)
legend (locator(1), legend = c("CLASICO", "ML","RML","WLS","OLS"), lwd=2, lty = c(0,6,6,6,6),  col=1:7)
#Esta imagen lo que nos indica que al parecer lo mejor son el de minimos cuadrados y ordinarios


#Podemos realizar todos los modelos programas 
vgm() #Ver los modelos que tiene la libreria