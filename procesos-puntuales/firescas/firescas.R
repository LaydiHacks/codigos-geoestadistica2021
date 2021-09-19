#================================================
#Cambiar el directorio de trabajo
#================================================
setwd(dir = 'firescas/')

#================================================
#paquetes necesarios
#================================================
library(splancs)
library(spatstat)
library(readr)
library(geospt)
library(RANN)
library("plotrix")
library(spdep)
library(SpatialTools)

#================================================
#Lectura Datos
#================================================
### Contorno de Castellón
incendiosCS <- read.csv("inc_01_06.csv",sep = ";")
contornoCS <- read.csv("cast_cont.csv",sep = ";")

xC <- contornoCS$X[633:1]
yC <- contornoCS$Y[633:1]
conpolyCast <- owin(poly = list(x = xC, y = yC))
### Filtro de incendios 2003 a 2006
ano2003a6CS <- incendiosCS[incendiosCS$fecha %in% c("2003","2004","2005","2006"), ]
### Patrón puntual en spatstat
Wildfires <- ppp(x = ano2003a6CS$x, y = ano2003a6CS$y, window = conpolyCast)



plot(Wildfires, main="Mapa de Eventos", pch=20)
points(Wildfires, pch=20, col='red')
summary(Wildfires)
contour(density(Wildfires), add=TRUE, main="Estimación kernel de la intensidad")



##################################################################
# Métodos basados en cuadrantes
##################################################################
linesX = 6
linesY = 6
Q<-quadratcount(Wildfires, nx=linesX,ny=linesY)
plot(Wildfires, main="Mapa de Eventos", col='white')
points(Wildfires, pch=20, col='gray')
plot(Q, add=TRUE, col='red')

plot(density(Wildfires), main="Estimación kernel de la intensidad")
plot(Q, add=TRUE, col=4)

plot(as.im(Wildfires,dimyx=linesX), main="Conteo por cuandrantes - imagen")

M<-quadrat.test(Wildfires, nx=linesX, ny=linesY)
plot(Wildfires, main="Mapa de Eventos", col='white')
points(Wildfires,  pch=20, col='gray')
plot(M, add=TRUE,  col='red')


par(mfrow=c(1,1))
persp(density(Wildfires),theta=-50,phi=30,shade=0.75,xlab="x",ylab="y",zlab="",cex.axis=0.7,cex.lab=0.7)



##################################################################
# Estimación Kernel Cuártico
##################################################################

#Generar puntos aleatorios para realizar predicción
xK <- runif(200, min(xC), max(xC))
yK <- runif(200, min(yC), max(yC))
aleatorios <- ppp(x = xK, y = yK, window = conpolyCast)
aleatorios <- as.ppp(aleatorios) # Eliminar los puntos fuera del contorno
plot(aleatorios)
points(aleatorios, pch=20, col='green')

data_xy <- ano2003a6CS[,8:9]

x11()
plot(Wildfires, pch=20)
points(Wildfires, pch=20, col='gray')
for (i in 1:length(aleatorios$x)) {
  So <- data.frame(t(c(aleatorios$x[i],aleatorios$y[i])))
  names(So) <- c("x","y")
  knn <- nn2(data_xy,So,k=6)
  dist_max = max(knn$nn.dists)
  draw.circle(aleatorios$x[i],aleatorios$y[i], dist_max, nv = 10000, border=2, lty=1, lwd=1 )
  points(So[1],So[2],col=3,pch=19)
  
  PP1 <- data.frame(data_xy[sort(knn$nn.idx),1:2],Dist=knn$nn.dists[order(knn$nn.idx)])
  LHS <-  3/(pi*dist_max^2)
  RHS <- (1-(PP1$Dist/(dist_max^2)))^2
  PP1$LHS.RHS <- LHS*RHS
  KE <- sum(PP1$LHS.RHS)# Intensidad kernel 
  aleatorios$IKE[i]  <- KE
}

summary(aleatorios$IKE)
hist(aleatorios$IKE)

##################################################################
# Knn
##################################################################

coordinates(ano2003a6CS) <- ~x+y
coords <- coordinates(ano2003a6CS)    
IDs <- row.names(as(ano2003a6CS, "data.frame"))

#======================================
# Construcción Función G: evento-evento
#======================================
M_kn1 <- knn2nb(knearneigh(coords, k=1), row.names=IDs)
plot(M_kn1, coords, col=4, lwd=2)
points(coords, pch = 19, col = "Red")
n = length(ano2003a6CS$x)
#Generar Cortes
M <- as.matrix(dist(coords))+diag(NA,n,n)
D <- apply(M,1,min,na.rm=T)

c1 <- 5000
c2 <- 10000
c3 <- 15000
c4 <- 20000
c5 <- 25000
c6 <- 30000
c7 <- 35000
c8 <- 40000
c9 <- 45000
c10 <- 50000
c11 <- 55000

I1  <- mean(ifelse(D<c1,1,0))
I2 <- mean(ifelse(D<c2,1,0))
I3 <- mean(ifelse(D<c3,1,0))
I4 <- mean(ifelse(D<c4,1,0))
I5 <- mean(ifelse(D<c5,1,0))
I6 <- mean(ifelse(D<c6,1,0))
I7 <- mean(ifelse(D<c7,1,0))
I8 <- mean(ifelse(D<c8,1,0))
I9 <- mean(ifelse(D<c9,1,0))
I10 <- mean(ifelse(D<c10,1,0))
I11 <- mean(ifelse(D<c11,1,0))

G <- data.frame(Distance=c(0,c1,c2,c3,c4,c5,c6, c7, c8, c9, c10,c11),Gest=c(0,I6,I2,I2,I3, I4, I6, I7, I8, I9,I10,I11))
plot(G,type="l", lwd=2.5)


#======================================
# Construcción Función F: punto-evento
#======================================

set.seed(125)
ptox <- runif(150, min(xC), max(xC))
ptoy <- runif(150, min(yC), max(yC))

plot(M_kn1, coords, col=2, lwd=2)
points(coords, pch = 19, col = "Red")

aleatorios_f <- ppp(x = ptox , y = ptoy, window = conpolyCast)
aleatorios_f <- as.ppp(aleatorios)
Puntos <- data.frame(x=aleatorios_f$x,y=aleatorios_f$y)
points(Puntos, pch = 19, col = "Black")

D1 <- dist2(as.matrix(Puntos),coords)

X <- apply(D1,1,function(x) {which(x == min(x))})
D2 <- apply(D1,1,min)

segments(Puntos[,1],Puntos[,2],coords[X,1],coords[X,2],col=4, lwd=2.5)

c1 <- 5000
c2 <- 10000
c3 <- 15000
c4 <- 20000
c5 <- 25000
c6 <- 30000
c7 <- 35000
c8 <- 40000
c9 <- 45000
c10 <- 50000
c11 <- 55000

In1  <- mean(ifelse(D2<c1,1,0))
In2 <- mean(ifelse(D2<c2,1,0))
In3 <- mean(ifelse(D2<c3,1,0))
In4 <- mean(ifelse(D2<c4,1,0))
In5 <- mean(ifelse(D2<c5,1,0))
In6 <- mean(ifelse(D2<c6,1,0))
In7 <- mean(ifelse(D2<c7,1,0))
In8 <- mean(ifelse(D2<c8,1,0))
In9 <- mean(ifelse(D2<c9,1,0))
In10 <- mean(ifelse(D2<c10,1,0))
In11 <- mean(ifelse(D2<c11,1,0))

F <- data.frame(Distance=c(0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11),Fest=c(0,In1,In2,In3, In4,In5,In6,In7,In8,In9,In10,In11))
plot(F,type="l",lwd=2.5)


##################################################################
# K-Funcion
##################################################################

kest(Wildfires)
plot(Kest(Wildfires), main = "K-estimada (Incendios)")
plot(pcf(Wildfires), main = "pair-correlation estimada (Incendios)")

     
     