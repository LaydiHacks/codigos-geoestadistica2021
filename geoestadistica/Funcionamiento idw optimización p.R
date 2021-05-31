###########################################################################
################    Funcionamiento IDW con optimizaci?n    ################
###########################################################################

source("funciones/geoestadistica/idw.cv.R")  
source("funciones/geoestadistica/graph.idw.R")        
library(pryr)
library(geospt)
library(sp)
library(gstat)
library(fields)
library(maptools)
library(rgdal)
library(sf)
                       
ejercIDW <- read.table("db/ejercIDW1.txt",header =T, sep="")
mapa <- readOGR("db/ARIARI.shp")

plot(ariari)
puntos <- spsample(ariari,20000,"regular")
plot(puntos)


idw.p <- idw(PRECI_TOT~ 1, ~ x+y, ejercIDW, puntos, nmax=15, nmin=15, idp=2)

pal2 <- colorRampPalette(c("snow3","royalblue1", "blue4"))

#Interpolaciones de Distancia Inversa\n Ponderada de la Precipitaci?n (P=2)
p1 <- spplot(idw.p[1], col.regions=pal2(100), cuts =60, scales = list(draw =T), xlab ="Este (m)", ylab = "Norte (m)", 
       main = "", auto.key = F)

split.screen( rbind(c(0, 1,0,1), c(1,1,0,1)))
split.screen(c(1,2), screen=1)-> ind
screen( ind[1])
p1
screen( ind[2])
image.plot(legend.only=TRUE, legend.width=0.5, col=pal2(100), smallplot=c(0.6,0.68, 0.5,0.75), 
           zlim=c(min(idw.p$var1.pred),max(idw.p$var1.pred)), axis.args = list(cex.axis = 0.7))
close.screen( all=TRUE)

idw.cv(PRECI_TOT~ 1, ~ x+y, ejercIDW, nmax=15, nmin=15, 2)
p.opt <- optimize(idw.cv, c(0,10), formula=PRECI_TOT~1, locations=~x+y, data=ejercIDW, nmax=15, nmin=15, progress=F)
# $minimum
#[1] 2.073498
#$objective
#[1] 26.16721

idw.graph <- graph.idw(PRECI_TOT~ 1, ~x+y, data=ejercIDW, np=50, p.dmax=4, nmax=15, nmin=15, P.T=T)


# Ejemplo preci
data(preci)
idw.graph <- graph.idw(prec~ 1, ~x+y, data=preci, np=50, p.dmax=4, nmax=10, nmin=10, P.T=T)
