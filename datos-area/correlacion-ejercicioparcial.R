########################################################################
############  SOLUCI?N SEGUNDO PARCIAL GEOESTADISTICA 2014  ############
########################################################################

#save.image("D:/CARLOS/Estadistica Espacial/Geoestadistica/Parciales/p2014.I")
#load("D:/CARLOS/Estadistica Espacial/Geoestadistica/Parciales/p2014.I")


library(shapefiles)
P1 <- data.frame(c(2.5,8.5,8.5,3,2.5,2.5),c(12,12,8,8,8,12))
colnames(P1)<-c("x","y")
P2 <- data.frame(c(8.5,14,14,8.5,8.5,8.5),c(10,10,4,4,8,10))
colnames(P2)<-c("x","y")
P3 <- data.frame(c(0,3,3,0,0),c(8,8,3,3,8))
colnames(P3)<-c("x","y")
P4 <- data.frame(c(3,8.5,8.5,8.5,3,3),c(8,8,4,3,3,8))
colnames(P4)<-c("x","y")
P5 <- data.frame(c(8.5,14,14,8.5,8.5,8.5,8.5),c(4,4,0,0,3,1,4))
colnames(P5)<-c("x","y")
P6 <- data.frame(c(5,8.5,8.5,5,5),c(3,3,1,1,3))
colnames(P6)<-c("x","y")

poligonos <- data.frame(c(rep(1,nrow(P1)),rep(2,nrow(P2)),rep(3,nrow(P3)),rep(4,nrow(P4)),rep(5,nrow(P5)),rep(6,nrow(P6))),
                        rbind(P1,P2,P3,P4,P5,P6))
names(poligonos) <- c("Id","x","y")
ddTable <- data.frame(Id=c(1,2,3,4,5,6),Name=c("P1","P2","P3","P4","P5","P6"))
ddShapefile <- convert.to.shapefile(poligonos, ddTable, "Id", 5)
write.shapefile(ddShapefile, "db/pol", arcgis=T)

library(spdep)
library(maptools)
library(rgdal)
par.poly1 <- readOGR("db", "pol")
plot(par.poly1)

# Contiguidad Efecto Reina: Orden 1 

par_nbq1 <- poly2nb(par.poly1)
par_nbr1 <- poly2nb(par.poly1,queen=F)           # torre
e.lwq <- nb2listw(par_nbq1, style="W")

library(spatialreg)
W.q <- as.matrix(as_dgRMatrix_listw(e.lwq))
e.lwr <- nb2listw(par_nbr1, style="W",zero.policy=T)  
W.r <- as.matrix(as_dgRMatrix_listw(e.lwr))
par.lags1 <- nblag(par_nbq1, 2)                  # Orden 2
e.lw2 <- nb2listw(par.lags1[[2]], style="W",zero.policy=T)
W.q2 <- as.matrix(as_dgRMatrix_listw(e.lw2))
# set.ZeroPolicyOption(T), para no escribir en todas las funciones: zero.policy=T 
# nb2listw(par.lags1[[2]], style="W")
# print(e.lw2)
set.ZeroPolicyOption(F)
e.lw2 <-nb2listw(par.lags1[[2]], style="W", zero.policy=T)
print(e.lw2, zero.policy=T)

x11()
plot(par.poly1)
coords1 <- coordinates(par.poly1)
plot(par_nbq1, coords1, add=T, col=4, lwd=2)

plot(par.poly1)
plot(par.lags1[[2]], coords1, add=T, col=4, lwd=2)

par.poly1@data$VarX <-  c(14,12,11,25,28,30)                         # par.poly@data$VarX
par.poly1@data$VarY <-  c(70,65,62,95,40,32)                         # par.poly@data$VarY

ddShapefile1 <- convert.to.shapefile(poligonos, par.poly1@data, "Id", 5)
write.shapefile(ddShapefile1, "D:/CARLOS/Estadistica Espacial/Geoestadistica/Parciales/pol1", arcgis=T)

aw <- print(nb2listw(par.lags1[[2]], zero.policy=TRUE),zero.policy=TRUE)
                                # aw
moran.test(par.poly1@data$VarX, e.lw2,zero.policy=TRUE)
W <- as.matrix(as_dgRMatrix_listw(e.lw1))
W1 <- as.matrix(as_dgRMatrix_listw(e.lw2))

source("D:/CARLOS/Estadistica Espacial/Funciones/Lattices/moran.bi.R")
moran.bi(par.poly1@data$VarX,par.poly1@data$VarX,e.lw2,zero.policy =T)
moran.bi(par.poly1@data$VarX,par.poly1@data$VarY,e.lw2,zero.policy =T)

source("D:/CARLOS/Estadistica Espacial/Funciones/Lattices/moranbi.test.R")
source("D:/CARLOS/Estadistica Espacial/Funciones/Lattices/randomize_vector.R")

set.seed(123)
MBP <- moranbi.test(par.poly1@data$VarX,par.poly1@data$VarY,e.lw1,999,graph=T,zero.policy =T,N=1000)
moranbi.test(par.poly1@data$VarX,par.poly1@data$VarX,e.lw2,999,graph=T,zero.policy =T,N=1000)

LM <- localmoran(par.poly1@data$VarX, e.lw2, zero.policy =T)

source("D:/CARLOS/Estadistica Espacial/Funciones/Lattices/moran.cluster.R")
# LISA Cluster Map
moran.cluster(par.poly1@data$VarX, e.lw2, zero.policy = T, par.poly1, significant=T)




source("E:/CARLOS/Estadistica Espacial/Funciones/Lattices/localmoran.bi.R")
LMB <- localmoran.bi(par.poly1@data$VarX, par.poly1@data$VarY, e.lw2, zero.policy =T)


