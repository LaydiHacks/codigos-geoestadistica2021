
#=======================================
# Construcción Función G: evento-evento
#=======================================

Muestra <- read.table("D:/CARLOS/Estadistica Espacial/Geoestadistica/Patrones Espaciales/Ejercicios R/DatosPP.txt", head=T)

library(spdep)
coordinates(Muestra) <- ~x+y
coords1 <- coordinates(Muestra)    
IDs <- row.names(as(Muestra, "data.frame"))

M_kn1 <- knn2nb(knearneigh(coords1, k=1), row.names=IDs)

#par(mar = c(0, 0, 0, 0), pty = "s")
plot(M_kn1, coords1, col=4, lwd=2)
points(coords1, pch = 19, col = "Red")
text(coords1[,1]+2,coords1[,2]+2,IDs,cex=1.2)

M <- as.matrix(dist(coords1))+diag(NA,12,12)
D <- apply(M,1,min,na.rm=T)

I6  <- mean(ifelse(D<6,1,0))
I12 <- mean(ifelse(D<12,1,0))
I18 <- mean(ifelse(D<18,1,0))
I24 <- mean(ifelse(D<24,1,0))
I30 <- mean(ifelse(D<30,1,0))
I35 <- mean(ifelse(D<36,1,0))

G <- data.frame(Distance=c(0,6,12,18,24,30,36),Gest=c(0,I6,I12,I18,I24,I30,I35))
plot(G,type="l", lwd=2.5)


#======================================
# Construcción Función F: punto-evento
#======================================


set.seed(125)
ptox <- runif(10,Muestra@bbox[1,1],Muestra@bbox[1,2])
ptoy <- runif(10,Muestra@bbox[2,1],Muestra@bbox[2,2])

plot(M_kn1, coords1, col=2, lwd=2)
points(coords1, pch = 19, col = "Red")
text(coords1[,1]+2,coords1[,2]+2,IDs,cex=1.2,col="Red")
Puntos <- data.frame(x=ptox,y=ptoy)
points(Puntos, pch = 19, col = "Black")
text(Puntos[,1]+2,Puntos[,2]+2,rownames(Puntos),cex=1.2)

library(SpatialTools)
D1 <- dist2(as.matrix(Puntos),coords1)

X <- apply(D1,1,function(x) {which(x == min(x))})
D2 <- apply(D1,1,min)

segments(Puntos[,1],Puntos[,2],coords1[X,1],coords1[X,2],col=4, 
lwd=2.5)

In3  <- mean(ifelse(D2<3,1,0))
In6 <- mean(ifelse(D2<6,1,0))
In9 <- mean(ifelse(D2<9,1,0))
In12 <- mean(ifelse(D2<12,1,0))
In15 <- mean(ifelse(D2<15,1,0))
In18 <- mean(ifelse(D2<18,1,0))
In21 <- mean(ifelse(D2<21,1,0))

F <- data.frame(Distance=c(0,3,6,9,12,15,18,21),Fest=c(0,In3,In6,In9,In12,In15,In18,In21))
plot(F,type="l",lwd=2.5)
