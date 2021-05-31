
                            #################################################
                            ###########  EJEMPLO:Precipitación ##############
                            #################################################
            
library(geoR)
library(gstat)                                                     # Iniciamos con el paquete gstat
library(geospt)
data(preci)
preci.geoR <- as.geodata(preci, coords.col = 2:3, data.col = 4)    # Objeto del tipo geodata (coordenadas y datos)
points.geodata(preci.geoR, x.leg=3, y.leg=5, main=c("Gráfico de Intensidades", "de Precipitación"), col.main=3, pt.div="quintile")
# Graficar los datos (puntos) y ver simbolos graduados, además la opción "add.to.plot" es útil para adicionar opciones de la instrucción "plot", como en este caso la grilla.
points.geodata(preci.geoR, x.leg=3, y.leg=5, main=c("Gráfico de Intensidades", "de Precipitación"), col.main=3, pt.div="quintile", add.to.plot = TRUE, panel.first = grid())          
                                   
# Valores de precipitación, desplazados para no superponerlos con los puntos
attach(preci)
text(x+0.12,y+0.12,prec,cex=0.6)
# 4 gráficas: configuración puntos en R2, perfiles E-W y S-N, e histograma de la función de densidad, utiles para ver la distribución espacial de los puntos, análizar tendencia en las direcciones E-W y S-N y analizar normalidad
plot(preci.geoR) 

# Para hacer gráficos en 3 dimensiones                                    
library(scatterplot3d)               
s3d<-scatterplot3d(x, y, prec, angle=45, main="Gráfico de Precipitación en 3D", col.main=4, xlab="Coordenada X", ylab="Coordenada Y", zlab="Precipitación")
# Adicionando un plano de regresión al gráfico anterior
my.lm <- lm(prec ~ x )
s3d$plane3d(my.lm, lty.box = "solid")        

# Resumen de Estadísticas:
summary(preci.geoR)
attach(preci)

#### Cálculo de la desviación estándar de la variable y ploteo de gráficos exploratorios
sd(prec)
par(mfrow=c(nrow=1, ncol=3))
hist(prec, freq=F, breaks=10, xlab="prec (mm)", ylab="norte")
curve(dnorm(x, mean=mean(prec), sd=sd(prec)), add=T)                  # X no funciona, solo minuscula
qqnorm(prec, ylab="prec (mm)", xlab="Cuantiles teóricos")
qqline(prec)                                                        #Agregarle al Q-Q Plot la línea media
boxplot(prec, main="BOX-PLOT", notch=F, horizontal=T, xlab="Precipitación (mm)")
points(mean(prec), y=1, pch=1, cex=2)
par(mfrow=c(nrow=1, ncol=1))

# Pruebas de normalidad de Shapiro-Wilk y Kolmogorov-Smirnov
# Verificación del supuesto de normalidad para los métodos lineales de kriging.
      # Opciones para tranformar:
      # boxcox.fit(geoR), para un 'vector' de datos
      # boxcox.geodata(), para un objeto de tipo 'geodata'

shapiro.test(prec)
n.q <- (prec-mean(prec))/sd(prec)            #Función prec
n.o <- order(n.q)                            # Lista con las posiciones de los datos ordenados
n.qo <- n.q[n.o]                             # Vector de cuantiles estandarizados y ordenados. Necesarios para prueba K-S
ks.test(n.qo, pnorm)                         # Le hago prueba K-S para saber si los datos provienen de una normal
ks.test(as.numeric(scale(sort(prec))),pnorm)

# Estimación del semivariograma:
detach(preci)
# La instrucción "vgm()", muestra las funciones de variograma disponibles en el programa geoR, se debe haber cargado el paquete gstat previamente.
vgm()
# Estimador clásico
prec.cl1<-variog(preci.geoR, option="cloud")              # opciones "bin" con uvec=1:7, utiles para semivariogramas promedio, que nos muestran con mayor claridad el modelo a usar
plot(prec.cl1, pch=20, main="Semivariograma Experimental Clásico")
grid()
prec.v1<-variog(preci.geoR,uvec=seq(0,7,0.5), bin.cloud = T )
plot(prec.v1,pch=20)
prec.vc1<-variog(preci.geoR,uvec=seq(0,7,1),bin.cloud=T)   # estimator.type ="modulus"
plot(prec.vc1,bin.cloud=T)


#### Nube de puntos del semivariograma
#Calcula el variograma de la muestra o residual (en el caso que se dé un modelo de tendencia)
                     
v <- variogram(prec~1, loc=~x+y, cutoff=1, cressie = F, data=preci, cloud=TRUE)      

           #########################################################################################
           ##    Kriging Simple: Funcionamiento Predicción un punto y estimación de su varianza   ##
           #########################################################################################

# En este caso la función de covarianza se debe construir a partir del vector de errores obtenidos de la diferencia entre las 
# observaciones y la media de la variable regionalizada poblacional, asumiendo previamente estacionariedad
# Para el ejemplo los parámetros de la función de covarianza son tomados del programa ArcGIS
# Funcionamiento con un vecindario de tamaño 5

xy <- preci[,2:3]
z <- preci$prec
mu = mean(z)                                                          # Media de la variable regionalizada
So <- c(3,4)                                                          # Coordenadas del punto a estimar
m.dist <- as.matrix(dist(rbind(xy,So)))                               # matriz de distancias
#D<-spDists(as.matrix(preci[,2:3]),longlat=F)
dist.So <- m.dist[nrow(m.dist),1:(ncol(m.dist)-1)]                    # vector de distancias al punto So.
vec.orden <- order(dist.So)                                           # vecinos ordenados
dist.vec.cerca <- dist.So[vec.orden[1:5]]                             # vecinos mas cercanos "5"
m.dist.vec <- as.matrix(dist(preci[,2:3]))[vec.orden[1:5], vec.orden[1:5]]
m.cov<-cov.spatial(m.dist.vec, cov.model= c("pure.nugget", "spherical"), cov.pars = rbind(c(9.2972,0),c(74.703, 3.573)))
v <- cov.spatial(dist.vec.cerca, cov.model= c("pure.nugget", "spherical"), cov.pars = rbind(c(9.2972,0),c(74.703, 3.573)))
cov.0 <- cov.spatial(0, cov.model= c("pure.nugget", "spherical"), cov.pars = rbind(c(9.2972,0),c(74.703, 3.573)))
Pesos.ks <- solve(m.cov)%*%v
KSpr <- mu+t(Pesos.ks)%*%(z[vec.orden[1:5]]-mu)
KSvr <- cov.0-t(Pesos.ks)%*%v

# Verificación con gstat
library(gstat)
library(sp)
pts = data.frame(xy, z=z)
coordinates(pts) <- c("x", "y")
So <- as.data.frame(cbind(x=So[1],y=So[2]))
coordinates(So) <- c("x", "y")
krige(z~1, pts, SpatialPoints(So), vgm(74.703, "Sph", 3.573, nugget=9.2972), nmax=5, beta = mu)
#krige(z~1, pts, SpatialPoints(So), vgm(74.703, "Sph", 3.573, nugget=9.2972))    # Kriging Ordinario
proc.time()


          ##########################################################################################
          ##  Kriging Ordinario: Funcionamiento Predicción un punto y estimación de su varianza   ##
          ##########################################################################################

one = rep(1,length(dist.vec.cerca))
cov.0 <- cov.spatial(0, cov.model= "spherical", cov.pars = c(112.33, 4.3441))
m.cov<-cov.spatial(m.dist.vec, cov.model= "spherical", cov.pars = c(112.33, 4.3441))
v <- cov.spatial(dist.vec.cerca, cov.model= "spherical", cov.pars = c(112.33, 4.3441))
m.cov.k <- rbind(cbind(m.cov, one),c(one,0))
m.cov.k.i <- solve(m.cov.k)
v.ko <- c(v,1)
#L<-solve(m.cov.k, v.ko)
Pesos.ko <- m.cov.k.i%*%v.ko
KOpr <- z[vec.orden[1:5]]%*%Pesos.ko[1:length(dist.vec.cerca)]
KOvr <- cov.0-v.ko%*%Pesos.ko
print(data.frame(So,KOpr,KOvr))                                            # Resultados punto So=(3,4)
krige(z~1, pts, SpatialPoints(So), vgm(112.33, "Sph", 4.3441), nmax=5)     
proc.time()


                             ###################################################
                             ######### Validación Cruzada KO y KS ##############
                             ###################################################

semiv.z <- variogram(z~1, ~x+y, data=preci, uvec=seq(0,7,1))            # Igual a "cressie=F" método clásico
plot(semiv.z, pl=T, ylab="Semivarianza", xlab="Distancia (m)")
sem.esf.sim <- vgm(112.33, "Sph", 4.3441)                               # Esférico para kriging Simple
sem.fit.esf.sim <- fit.variogram(semiv.z, sem.esf.sim, fit.method = 0)  # para dejar el mismo modelo de ArcGIS
sem.esf.ord <- vgm(74.703, "Sph", 3.573, nugget=9.2972)
sem.fit.esf.ord <- fit.variogram(semiv.z, sem.esf.ord, fit.method = 0)  # para dejar el mismo modelo de ArcGIS
precipi <- data.frame(preci[,2:4])
KO.esf.cv.z <- krige.cv(prec~1, ~x+y, precipi, sem.fit.esf.ord, nmin=0, nmax=5)    
KS.esf.cv.z <- krige.cv(prec~1, ~ x+y, precipi, sem.fit.esf.sim, beta=mu, nmin=0, nmax=5)

# Correr primero la función criterio.cv

resultados.cv.z <- rbind(criterio.cv(KO.esf.cv.z), criterio.cv(KS.esf.cv.z))
rownames(resultados.cv.z) <- c("KO.esf.cv.z", "KS.esf.cv.z")
resultados.cv.z


            ########################################################################################
            ## Generación de los mapas de pronóstico y de errores para Kriging Simple y Ordinario ##
            ########################################################################################

puntos<-expand.grid(x=seq(min(x),max(x),0.02), y=seq(min(y),max(y),0.02))
plot(puntos)
coordinates(puntos) = c("x", "y")
gridded(puntos) <- TRUE
pron.pts.ko <- krige(z~1, pts, puntos, vgm(112.33, "Sph", 4.3441), nmax=5)     # Kriging Ordinario
pron.pts.ks <- krige(z~1, pts, puntos, vgm(74.703, "Sph", 3.573, nugget=9.2972), nmax=5, beta = mu)                                               # Kriging Simple

l2 = list("sp.points", So, pch = 3, col = "grey")
p1 <- spplot(pron.pts.ko, "var1.pred", main="Precipitación Media Anual \nPredicciones Kriging Ordinario", col.regions=bpy.colors(100), cuts=10, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", sp.layout=list(l2), key.space=list(space="right", cex=0.6))
p2 <- spplot(pron.pts.ko, "var1.var", main="Precipitación Media Anual \nErrores Kriging Ordinario", col.regions=bpy.colors(100), cuts=10, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
p3 <- spplot(pron.pts.ks, "var1.pred", main="Precipitación Media Anual \nPredicciones Kriging Simple", col.regions=bpy.colors(100), cuts=10, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", sp.layout=list(l2), key.space=list(space="right", cex=0.6))
p4 <- spplot(pron.pts.ks, "var1.var", main="Precipitación Media Anual \nErrores Kriging Simple", col.regions=bpy.colors(100), cuts=10, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
print(p1, split = c(1, 1, 2, 2), more = T)
print(p2, split = c(2, 1, 2, 2), more = T)
print(p3, split = c(1, 2, 2, 2), more = T)
print(p4, split = c(2, 2, 2, 2), more = F)

library(maptools)
writeAsciiGrid(pron.pts.ko, "ko.asc")
write.asciigrid(pron.pts.ko, "ko.asc")

          ##########################################################################################
          ##  Kriging Universal: Funcionamiento Predicción un punto y estimación de su varianza   ##
          ##########################################################################################

# En este caso la función de covarianza se debe construir a partir del vector de errores obtenidos de la diferencia entre las observaciones y los pronósticos obtenidos de un modelo de deriva (o tendencia, es decir se debe correr una regresión de un polinomio de orden 1, 2 o 3 "depende de la cantidad de puntos que se disponga", y en esta regresión la variable regionalizada Z es la variable respuesta y las variables asociadas a las coordenadas X y Y son las variables independientes), asumiendo previamente no estacionariedad fuerte, pero si estacionariedad intrinseca o débil.
# Para el ejemplo los parámetros de la función de covarianza son tomados del programa ArcGIS, esto solo con el objetivo de realizar la comparación entre R y ArcGIS.


So <- c(3,4)
tendencia <- xy[vec.orden[1:5],]
one = rep(1,length(dist.vec.cerca))
cov.0 <- cov.spatial(0, cov.model= "spherical", cov.pars = c(53.064, 2.8858))
m.cov<-cov.spatial(m.dist.vec, cov.model= "spherical", cov.pars = c(53.064, 2.8858))
v <- cov.spatial(dist.vec.cerca, cov.model= "spherical", cov.pars = c(53.064, 2.8858))
m.cov.k <- rbind(as.matrix(cbind(m.cov, one, tendencia)),as.matrix(cbind(rbind(one,t(tendencia)),matrix(0,ncol=(ncol(tendencia)+1),nrow=(ncol(tendencia)+1)))))
m.cov.k.i <- solve(m.cov.k)
v.ko <- c(v,1,So)
#L<-solve(m.cov.k, v.ko)
Pesos.ko <- m.cov.k.i%*%v.ko
KUpr <- z[vec.orden[1:5]]%*%Pesos.ko[1:length(dist.vec.cerca)]
KUvr <- cov.0-v.ko%*%Pesos.ko
So = as.data.frame(cbind(x=So[1],y=So[2]))
coordinates(So)<- c("x", "y")
print(data.frame(So,KUpr,KUvr))                                                 # Resultados punto So=(3,4)

# Kriging Universal Orden 1
krige(z~x+y, pts, SpatialPoints(So), vgm(53.064, "Sph", 2.8858), nmax=5)
# Kriging Universal Orden 2        
krige(z~x + y + x*y + I(x^2)+I(y^2), pts, SpatialPoints(So), vgm(19.201, "Sph", 1.5823), nmax=6)            
# A partir del orden 3 tendríamos más parámetros que ecuaciones, o si consideramos menos observaciones en el vecindario con un orden 2 tendriamos varianzas negativas.
krige(z~x + y + x*y + I(x^2)+I(y^2), pts, SpatialPoints(So), vgm(19.201, "Sph", 1.5823), nmax=5)



                         ###########################################################
                         ########### Validación Cruzada KU Ordenes 1 y 2 ###########
                         ###########################################################

# Esférico para kriging Universal orden 1
sem.esf.univ1 <- vgm(53.064, "Sph", 2.8858)                           
sem.fit.esf.univ1 <- fit.variogram(semiv.z, sem.esf.univ1, fit.method = 0)  
# Esférico para kriging Universal orden 2
sem.esf.univ2 <- vgm(19.201, "Sph", 1.5823)
sem.fit.esf.univ2 <- fit.variogram(semiv.z, sem.esf.univ2, fit.method = 0)  
precipi<-data.frame(preci[,2:4])
KU1.esf.cv.z <- krige.cv(prec~x+y, ~x+y, precipi, sem.fit.esf.univ1, nmin=0, nmax=5)    
KU2.esf.cv.z <- krige.cv(prec~x + y + x*y + I(x^2)+I(y^2), ~ x+y, precipi, sem.fit.esf.univ2, nmin=0, nmax=6)


resultados.cv1.z <- rbind(criterio.cv(KU1.esf.cv.z), criterio.cv(KU2.esf.cv.z))
rownames(resultados.cv1.z) <- c("KU1.esf.cv.z", "KU2.esf.cv.z")
resultados.cv1.z
# Se recuerda que 0 <= R2 <= 1, debido a la validación cruzada en este caso fue menor de 0, 
# dicho criterio no se aplica para la decisión. Además se evidencia que entre más se aleja del intervalo [0,1], 
# mayor es el error.



          ###########################################################################################
          ## Generación de los mapas de pronóstico y de errores para Kriging Universal Orden 1 y 2 ##
          ###########################################################################################

puntos<-expand.grid(x=seq(min(x),max(x),0.02), y=seq(min(y),max(y),0.02))
plot(puntos)
coordinates(puntos) = c("x", "y")
gridded(puntos) <- TRUE
# Kriging Universal Orden 1
pron.pts.ku1 <- krige(z~x+y, pts, SpatialPoints(puntos), vgm(53.064, "Sph", 2.8858), nmax=5)
# Kriging Universal Orden 2     
pron.pts.ku2 <- krige(z~x + y + x*y + I(x^2)+I(y^2), pts, SpatialPoints(puntos), vgm(19.201, "Sph", 1.5823), nmax=6)                                               



l2 = list("sp.points", So, pch = 3, col = "grey")
rw.colors <- colorRampPalette(c("red", "yellow"))
p1 <- spplot(pron.pts.ku1, "var1.pred", main="Prec. Media Anual Predicciones \nKriging Universal Orden 1", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", sp.layout=list(l2), key.space=list(space="right", cex=0.6))
p2 <- spplot(pron.pts.ku1, "var1.var", main="Prec. Media Anual Errores \nKriging Universal Orden 1", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
p3 <- spplot(pron.pts.ku2, "var1.pred", main="Prec. Media Anual Predicciones \nKriging Universal Orden 2", col.regions=rw.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", sp.layout=list(l2), key.space=list(space="right", cex=0.6))
p4 <- spplot(pron.pts.ku2, "var1.var", main="Prec. Media Anual Errores \nKriging Universal Orden 2", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
print(p1, split = c(1, 1, 2, 2), more = T)
print(p2, split = c(2, 1, 2, 2), more = T)
print(p3, split = c(1, 2, 2, 2), more = T)
print(p4, split = c(2, 2, 2, 2), more = F)

##########################################################################################
##  Kriging Indicador: Funcionamiento Predicción un punto y estimación de su varianza   ##
##########################################################################################
detach(preci)
puntos<-expand.grid(x=seq(min(preci$x),max(preci$x),0.02), y=seq(min(preci$y),max(preci$y),0.02))
plot(puntos)
coordinates(puntos) = c("x", "y")
gridded(puntos) <- TRUE
pron.pts.ki <- krige(I(preci$prec<=420)~1, pts, puntos, vgm(0.2831, "Sph", 6.9642, nugget=0.0929), nmax=10)     # Kriging Ordinario
pki <- spplot(pron.pts.ki, "var1.pred", main="Prec. Media Anual Predicciones \nKriging Indicador", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
pkiv <- spplot(pron.pts.ki, "var1.var", main="Prec. Media Anual Errores \nKriging Indicador", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))



                             ######################################################
                             #######        Métodos Deterministicos:        #######
                             #######    IDW (Distancia Inversa Ponderada)   #######
                             ######################################################

library(sp)                                                   
So <- c(3,4)                                                 # Coordenadas del punto a estima
So <- as.data.frame(cbind(x=So[1],y=So[2]))
coordinates(So) <- c("x", "y")
m.dist <- as.matrix(dist(rbind(xy,So)))                      # matriz de distancias
# dist.So <- spDists(So,xy)                                  # vector de distancias al punto So.
dist.So <- m.dist[nrow(m.dist),1:(ncol(m.dist)-1)]           # vector de distancias al punto So.
vec.orden <- order(dist.So)                                  # vecinos ordenados
dist.vec.cerca <- dist.So[vec.orden[1:5]]                    # vecinos mas cercanos "5"
factor.p <- 2
Peso.i <- dist.vec.cerca^(-factor.p)/sum(dist.vec.cerca^(-factor.p))
Pronostico.So <- z[vec.orden[1:5]]%*%Peso.i

library(gstat)
pts = data.frame(xy, z=z)
coordinates(pts) = c("x", "y")
So = as.data.frame(cbind(x=So[1],y=So[2]))
coordinates(So)<- c("x", "y")
# P=2, por defecto para idw 
prec.idw <- krige(z~1, pts, SpatialPoints(puntos), nmax=5)                   
prec.idwSo <- krige(z~1, pts, So, nmax=5)
spplot(prec.idw["var1.pred"], col.regions=bpy.colors(100), main = "Interpolaciones de Distancia Inversa\n Ponderada de la Precipitación", cuts=10, key.space=list(space="right", cex=0.8))

# por defecto considera todas
idw.int1 <- idw(prec~ 1, ~ x+y, preci, puntos, nmax=10, nmin=10, idp=1)      
idw.intSo1 <- idw(prec~ 1, ~ x+y, preci, So, nmax=10, nmin=10, idp=1)            
spplot(idw.int1[1], col.regions=bpy.colors(100), scales = list(draw =T), xlab ="Este (m)", ylab = "Norte (m)", main = "Interpolaciones de Distancia Inversa\n Ponderada de la Precipitación (P=1)")  
# "rev", util para invertir el orden de los colores


idw.int2 <- idw(prec~ 1, ~ x+y, preci, puntos, nmax=5, nmin=5, idp=2)
idw.intSo2 <- idw(prec~ 1, ~ x+y, preci, So, nmax=5, nmin=5, idp=2)           
spplot(idw.int2[1], col.regions=bpy.colors(100), scales = list(draw =T), xlab ="Este (m)", ylab = "Norte (m)", main = "Interpolaciones de Distancia Inversa\n Ponderada de la Precipitación (P=2)")

                           ####################################################################
                           ########       OPTIMIZACIÓN PARA SELECCIONAR P EN IDW     ##########
                           ####################################################################

p.optimo <- function(p, formula, locations, data, newdata, nmax, nmin, maxdist, var.reg){
idw.pred <- as.data.frame(matrix(NA,nrow= nrow(data), ncol=4))
colnames(idw.pred) <- c("x","y","var1.pred","var1.var")
for(i in 1:(nrow(data))){
idw.pred[i,] <- idw(formula, locations, data[-i,], newdata[i,], nmax, nmin, maxdist, idp=p)
} 
RMSPE <-  sqrt(sum((idw.pred$var1.pred-var.reg)^2)/nrow(data))
RMSPE
}

P <- optimize(p.optimo, c(0,10), formula=prec~1, locations=~x+y, data=preci, newdata=preci, nmax=10, nmin=10, maxdist=Inf, var.reg=prec)
cat("Parámetro Óptimo IDW: ", "\n", "P       = ", P$minimum, "\n", "RMSPE   = ", P$objective, "\n")

P.opt <- optim( par=1 , fn=p.optimo, gr= "Nelder-Mead", method = "Nelder-Mead", formula=prec~1, locations=~x+y, data=preci, newdata=preci, nmax=10, nmin=10, maxdist=Inf, var.reg=prec)
P.opt <- optim( par=1 , fn=p.optimo, method = "L-BFGS-B", formula=prec~1, locations=~x+y, data=preci, newdata=preci, nmax=10, nmin=10, maxdist=Inf, var.reg=prec)


p1 <- p.optimo(p=1, formula=prec~1, locations=~x+y, data=preci, newdata=preci, nmax=10, nmin=10, maxdist= Inf, var.reg=prec)
p2 <- p.optimo(p=2, formula=prec~1, locations=~x+y, data=preci, newdata=preci, nmax=10, nmin=10, maxdist=Inf, var.reg=prec)
p3 <- P$objective
p4 <- p.optimo(p=3, formula=prec~1, locations=~x+y, data=preci, newdata=preci, nmax=10, nmin=10, maxdist=Inf, var.reg=prec)
p5 <- p.optimo(p=4, formula=prec~1, locations=~x+y, data=preci, newdata=preci, nmax=10, nmin=10, maxdist=Inf, var.reg=prec)
p6 <- p.optimo(p=5, formula=prec~1, locations=~x+y, data=preci, newdata=preci, nmax=10, nmin=10, maxdist=Inf, var.reg=prec)


RMSPE <- c(p1, p2, p3, p4, p5, p6)
plot(c(1,2,P$minimum,3,4,5),RMSPE, main="Gráfico de optimización del Parámetro (P)\n Distancia Inversa Ponderada", ylab="RMSPE", xlab="p óptimo = 2.96997", type="l")

rmspe <- function(p){
idw.pred <- as.data.frame(matrix(NA,nrow= nrow(preci), ncol=4))
colnames(idw.pred) <- c("x","y","var1.pred","var1.var")
for(i in 1:(nrow(preci))){
idw.pred[i,] <- idw(prec~ 1, ~ x+y, preci[-i,], preci[i,], nmax=10, nmin=10, idp=p)    # optimo 2.9701
} 
RMSPE <-  sqrt(sum((idw.pred$var1.pred-prec)^2)/nrow(preci))
RMSPE
}
P <- optimize(rmspe,c(1,10))
p1 <- rmspe(p=0.5); p2 <- rmspe(p=1); p3 <- rmspe(p=1.5); p4 <- rmspe(p=2); p5 <- rmspe(p=2.5); p6 <- rmspe(p=2.96997); p7 <- rmspe(p=3); p8 <- rmspe(p=3.5); p9 <- rmspe(p=4); p10 <- rmspe(p=4.5); p11 <- rmspe(p=5); p12 <- rmspe(p=5.5); p13 <- rmspe(p=6)
RMSPE <- c(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13)
plot(c(0.5, 1, 1.5, 2, 2.5, P$minimum, 3, 3.5, 4, 4.5, 5, 5.5, 6),RMSPE, main="Gráfico de optimización del Parámetro (P)\n Distancia Inversa Ponderada", ylab="RMSPE", xlab="p óptimo = 2.96997", type="l")




          #########################################################################################
          ##############                  MÉTODOS DETERMINISTICOS:                  ###############
          ##############                  FUNCIONES DE BASE RADIAL                  ###############
          #########################################################################################


puntos<-expand.grid(x=seq(min(x),max(x),0.02), y=seq(min(y),max(y),0.02))
plot(puntos)

pred.rbf.m <- rbf(eta=0.2589, z=preci$prec, coordinates=preci[,2:3], newdata=puntos, n.neigh=9, func="M")
coordinates(pred.rbf.m) = c("x", "y")
gridded(pred.rbf.m) <- TRUE
# muestra el mapa de prediccion con la rbf multicuadratica
map.m <- spplot(pred.rbf.m["var1.pred"], cuts=40, scales = list(draw =T), col.regions=bpy.colors(100), key.space=list(space="right", cex=0.8))

pred.rbf.tps <- rbf(eta=0.1461, z=preci$prec, coordinates=preci[,2:3], newdata=puntos, n.neigh=9, func="TPS")
coordinates(pred.rbf.tps) = c("x", "y")
gridded(pred.rbf.tps) <- TRUE
# muestra el mapa de prediccion con la rbf TPS
map.tps <- spplot(pred.rbf.m["var1.pred"], cuts=40, scales = list(draw =T), col.regions=bpy.colors(100), key.space=list(space="right", cex=0.8))

load("D:/CARLOS/Doctorado/Investigación 2008-2009/Diseño e Implementación de Software Estadístico/Resultados/Metodos Kriging/Final/Envio/Kriging y FBR")
# save.image("D:/CARLOS/Doctorado/Investigación 2008-2009/Diseño e Implementación de Software Estadístico/Resultados/Metodos Kriging/Final/Envio/Kriging y FBR")