#######################################################################
############       EJERCICIO DATOS AREAS "COLUMBUS"        ############
#######################################################################
#https://arcruz0.github.io/libroadp/mapas.html
#http://rspatial.org/analysis/rst/1-introduction.html
#http://www.econ.uiuc.edu/~lab/workshop/Spatial_in_R.html
#http://rstudio-pubs-static.s3.amazonaws.com/5009_aa49bd54e40f41f19352af3d9e00b036.html
# clear the workspace, plots and console
rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014")

#save.image("D:/CARLOS/Econometria Espacial/Doctorado Econom�a UN/Clases/SAR")
#load("D:/CARLOS/Econometria Espacial/Doctorado Econom�a UN/Clases/SAR")
   
#####################################################
# 0. ACCESO A LOS DATOS Y CALCULO PESOS ESPACIALES  #
#####################################################
library(spgwr)
library(adespatial)
library(raster)
library(tmap)
library(lmtest)
library(RColorBrewer)
library(classInt)
library(spdep)
library(sphet)
library(pgirmess)
library(rgdal)
library(spatialreg)
#library(RcmdrPlugin.epack)   #no cargar ya que:     
# R tiene el molesto problema de que no permite cargar m�s de 100 DLL.

columbus <- sf::read_sf(system.file("etc/shapes/columbus.shp",package="spdep"))
plot(columbus$geometry)
plot(st_geometry(columbus), border="blue")
compassRose(6.5,14,rot=0,cex=0.5)

columbus <- readOGR(system.file("shapes/columbus.shp", package="spData")[1])
proj4string(columbus) <- CRS("+proj=longlat +datum=WGS84")
columbus <- spTransform(columbus, "+proj=longlat")
plot(columbus)

coords <- coordinates(columbus)
rn <- sapply(slot(columbus, "polygons"), function(x) slot(x, "ID"))
k1 <- knn2nb(knearneigh(coords))
all.linked <- max(unlist(nbdists(k1, coords)))
col.nb.0.all <- dnearneigh(coords, 0, all.linked, row.names=rn)
summary(col.nb.0.all, coords)
plot(columbus, border="blue")
plot(col.nb.0.all, coords, add=TRUE)
title(main=paste("Distance based neighbours 0-", format(all.linked), " distance units", sep=""))

# Identificaci�n vecinos
par(mai=c(0,0,0,0))
set.seed(1)
plot(columbus, col=sample(rainbow(500))[1:49])
xy <- coordinates(columbus)
points(xy, cex=3, pch=20, col='white')
text(columbus, 'POLYID', cex=0.6)
detach(package:raster)

#######################################################
# 1. ACCESO A LOS DATOS Y CALCULO PESOS ESPACIALES
#######################################################
col_nbq1 <- poly2nb(columbus) 
col.lags <- nblag(col_nbq1, 9)
 
col.poly <- columbus

# col.poly <- readShapePoly("D:/GeoDA/columbus.shp")       # Seleccionar el directorio previamente de GeoDA
brks <- quantile(col.poly$CRIME, seq(0, 1, 1/max(col.poly$CRIME)))
cols <- grey((length(brks):2)/length(brks))
# my.colors <- rev(c("red", "green", "blue","limegreen","brown","darkgrey","yellow","purple","orange","cyan","pink", "white"))    
# asigna colores diferentes, pero sin quantil
l2 = list("SpatialPolygonsRescale", layout.north.arrow(), offset =  c(6.1,13.5), scale = 0.5)   # Layout
l3 = list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(6.6,13.5),  scale = 0.5, fill=c("transparent","black"))
l4 = list("sp.text", c(6.6,13.7), "0")
l5 = list("sp.text", c(7.1,13.7), "500 m")
spplot(col.poly["CRIME"], scales=list(draw=TRUE), col.regions=cm.colors(20),  sp.layout=list(l2,l3,l4,l5))
# rainbow(n, start=.7, end=.1), heat.colors(n), terrain.colors(n), topo.colors(n), cm.colors(n)

source("Geoestadistica/codigos-geoestadistica2021/funciones/northarrow.r")
source("Geoestadistica/codigos-geoestadistica2021/funciones/scalebar.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/moran.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/geary.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/moranbi.test.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/randomize_vector.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/moran.cluster.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/moran.bi.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/moran.cluster.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/getis.cluster.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/localmoran.bi.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/moranbi.plot.R") 
source("Geoestadistica/codigos-geoestadistica2021/funciones/quantile.e.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/sp.na.omit.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/correlogram.d.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/sp.correlogram.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/spcorrelogram.bi.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/moran.bi1.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/moranbi1.test.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/geary.bi.R")
source("Geoestadistica/codigos-geoestadistica2021/funciones/test.w.R")
       
#https://rpubs.com/chrisbrunsdon/114718
#http://stackoverflow.com/questions/36251653/r-program-map-is-small
#https://www.rdocumentation.org/packages/spdep/versions/0.6-9/topics/trW
#http://www.sthda.com/english/wiki/colors-in-r
#display.brewer.all()
tm_shape(columbus) +
    tm_polygons(c("CRIME", "HOVAL","INC"), 
        style=c("pretty", "kmeans","pretty"),
        palette=list("RdYlGn", "RdBu","YlOrRd"),
        auto.palette.mapping=FALSE,
        title=c("CRIME RATE", "HOVAL", "INC"))
detach(package:tmap,unload=TRUE)

# bubble plot equal-frequency class intervals
plotvar <- round(col.poly@data$CRIME,3)
nclr <- 5
plotclr <- brewer.pal(nclr,"YlOrBr")
#plotclr <- plotclr[nclr:1] # reorder colors if appropriate
max.symbol.size=6
min.symbol.size=1
class <- classIntervals(plotvar, nclr, style="quantile")
colcode <- findColours(class, plotclr)
symbol.size <- ((plotvar-min(plotvar))/
   (max(plotvar)-min(plotvar))*(max.symbol.size-min.symbol.size)
   +min.symbol.size)

Data <- col.poly@data
plot(col.poly)
col.poly.cntr <- coordinates(col.poly)
points(col.poly.cntr, pch=16, col=colcode, cex=symbol.size)
points(col.poly.cntr, cex=symbol.size)
text(8.1,15.2, "Area: Equal-Frequency Class Intervals")
legend(locator(1), legend=names(attr(colcode, "table")), fill=attr(colcode, "palette"), cex=0.8, bty="n")


coordinates(Data) <- c("X","Y")
bubble(Data,"CRIME",maxsize=2,xlim = col.poly@bbox[,1],	ylim = col.poly@bbox[,2], axes=F) 
#plot(columbus,add=T)          
           
plotvar <- round(col.poly@data$CRIME,1)    
nclr <- 5
plotclr <- brewer.pal(nclr,"YlOrBr")
class <- classIntervals(plotvar, nclr, style="quantile")
colcode <- findColours(class, plotclr)
plot(col.poly)   
plot(col.poly, col=colcode, add=T)
# title(main="Tasa de Criminalidad", sub="Cuantiles (Igual-Frecuencia) Intervalos de Clase")
# mtext("Cuantiles (Igual-Frecuencia) Intervalos de Clase", side=1)
legend(locator(1), legend=names(attr(colcode, "table")), fill=attr(colcode, "palette"), cex=0.8, bty="n")
northarrow(c(10.1,14.5),0.2,cex=0.6)
scalebar(c(6.1,10.7),length=1,unit="m",division.cex=0.6)

quantile.e(x=col.poly@data$CRIME,ic=5,digits=1,style="quantile",border=col.poly,size=0.3, xn=10.1,yn=14.5,xsb=6.1,ysb=10.7,length=1,unit="m")

#scale <- sqrt(columbus$CRIME/max(columbus$CRIME))*70
#scale <- scale[order(columbus$CRIME, decreasing = TRUE)]
#symbols(coords[,1:2][order(columbus$CRIME, decreasing = TRUE),],circle=scale/200,bg=grey(.8),fg=grey(.5),
#inches=FALSE,add=TRUE)

                      
#o.nb <- data(columbus.gal)                      # Archivo de pesos generado en GeoDA             
o.nb <- read.gal("C:/Program Files/R/R-4.0.0/library/spdep/etc/weights/columbus.gal")    
o.nb[[1]]                                                 # Vecinos al poligono 1

o.nb[[49]]                                                # Vecinos al poligono 49

# Correlograma Moran a partir de matriz contiguidad espacial
sp.cr <- sp.correlogram(o.nb, col.poly@data$CRIME, order=7, method="C", style="W", zero.policy=T)
x11()
plot(sp.cr)
cor.s <- sp.correlogram(o.nb, col.poly@data$CRIME, order=9, method="I", style="W", zero.policy=T)
plot(cor.s)

a.lw <- nb2listw(o.nb, style="W")
mbi <- moranbi1.test(x=columbus@data$CRIME,y=columbus@data$INC,a.lw,zero.policy =T,randomisation =T,
              alternative="two.sided",adjust.n=TRUE)
mbi

corbi.ci <- spcorrelogram.bi(col_nbq1, col.poly@data$CRIME, col.poly@data$INC, order=9, 
                             method="I", style="W", randomisation =T, zero.policy=T,alternative="two.sided")
corbi.ci
plot(corbi.ci)
plot.spcorbi(corbi.ci,main="Bivariate CRIME-INCOME")

corbi.ch <- spcorrelogram.bi(o.nb, col.poly@data$CRIME, col.poly@data$HOVAL, order=9, 
                             method="I", style="W", zero.policy=T)
corbi.ch
plot(corbi.ch)
plot(corbi.ch,main="Bivariate CRIME-HOVAL")

corbig.ci <- spcorrelogram.bi(o.nb, col.poly@data$CRIME, col.poly@data$INC, order=7, 
                             method="C", style="W", zero.policy=T)
plot(corbig.ci,main="Bivariate CRIME-INC")

# Correlograma Moran a partir de matriz distancias 
corD <- correlog(coords,col.poly@data$CRIME,method="Moran",nbclass=20)
corD
plot(corD)

cord <- correlogram.d(coords,col.poly@data$CRIME,method="Moran",zero.policy=T,nbclass=20)
cord$res
plot(cord$res)

#####pesos K-Vecino####
#k.nb <- read.gwt2nb("atl_hom.gwt") 
## Generar opciones para poder escojer despues la mejor matriz W
IDs <- row.names(as(col.poly, "data.frame"))
k1 <- knn2nb(knearneigh(coords))
col.k1 <-nb2listw(k1, style="W")
k2 <- knn2nb(knearneigh(coords,2))
col.k2 <-nb2listw(k2, style="W")
k3 <- knn2nb(knearneigh(coords,3))
col.k3 <-nb2listw(k3, style="W")
k4 <- knn2nb(knearneigh(coords,4))
col.k4 <-nb2listw(k4, style="W")
k5 <- knn2nb(knearneigh(coords,5))
col.k5 <-nb2listw(k5, style="W")
k6 <- knn2nb(knearneigh(coords,6))
col.k6 <-nb2listw(k6, style="W")

## Generar el correlograma, recordar que cuando las barras tocan el cero, no son significativas
## Esto es importante porque podemos identificar que terminos podemos ingresar en el modelo, 
## en el analisis confirmatario podemos probar esos terminos. 
col.kn7<-knn2nb(knearneigh(coords, k=7), row.names=IDs) 
col.knlags <- nblag(col.kn7, 7)
sp.cr <- sp.correlogram(col.kn7, col.poly@data$CRIME, order=7, method="corr", style="W", zero.policy=T)
cor <- sp.correlogram(k6, col.poly@data$CRIME, order=6, method="I", style="W", zero.policy=T)
plot(cor)

# Adem�s del estilo "W", es posible crear los pesos de acuerdo con los estilos "B", "C", "U" 
# y "S". El estilo "B" es la codificaci�n binaria b�sica, "C" es la estandarizaci�n global, 
# "U" es igual a "C" dividido el n�mero de vecinos, y "S" es el esquema de codificaci�n en 
# que se estabiliza la varianza (Tiefelsdorf et al. 1999).
                  
o.nb <- read.gal("C:/Program Files/R/R-4.0.0/library/spdep/etc/weights/columbus.gal")    
# Si hacemos esto es necesario revizar que las conexiones quedaron Bien
# porque hay casos, en donde las importantes no quedan tan bien 
a.lw <- nb2listw(o.nb, style="W")
a.lwb <- nb2listw(o.nb, style="B")
a.lw$weights[[1]]                               # Pesos de los vecinos del poligono 1
a.lw$weights[[49]]                              # Pesos de los vecinos del poligono 49

# matriz de vecindad en t�rminos de unos y ceros

bin.matrix <- matrix(0, ncol=length(o.nb), nrow=length(o.nb))
              for (i in 1:length(o.nb)){
              for (k in 1:length(o.nb[[i]]))
              {
              bin.matrix[i, o.nb[[i]][k]] <- 1
              }
              }
# �C�mo ser�a a partir de lo anteriormente visto?
bin.matrix.1<-as_dgRMatrix_listw(a.lwb)
              
# A continuaci�n se presentan solo las primeras seis filas y las primeras 14 columnas de esa matriz: 

bin.matrix[1:6, 1:14]

# Visualmente se obtiene una imagen como la presentada mediante este c�digo:
image(as_dgRMatrix_listw(a.lw))                   # Ponderaci�n estandarizada por fila
image(as_dgRMatrix_listw(a.lwb))                  # Ponderaci�n binaria
image(bin.matrix)

#########################################################
# 2. C�LCULO DE LA AUTOCORRELACI�N ESPACIAL Y GR�FICO
#########################################################

# El �ndice I se obtiene a partir de la Ecuaci�n (1):
CRIME <- as.data.frame(col.poly)$CRIME
x <- c(CRIME)
# Se calcula el valor promedio de la variable, y se resta a cada observaci�n la media,
# se eleva al cuadrado y se calcula la suma:

x.media <- mean(x)
x.centrado <- x - x.media
scx <- sum(x.centrado^2)                          # Suma de Cuadrados de x (variable crimen)
scx1 <- var(x)*(length(x)-1)                      # Suma de Cuadrados de x (variable crimen), m�s r�pido.
crossprod(x.centrado)
            
# c�lculo de la variable de rezagos, WY, que se obtiene a partir de la matriz de pesos. Para ello, 
# en cada pol�gono se multiplica el valor de la variable CRIME de cada uno de sus vecinos por su 
# respectivo peso y luego se suman dichos productos.

wy <- matrix(0, ncol=1, nrow=length(CRIME))
      for (i in 1:length(CRIME)){
      vec <- matrix(0, nrow=length(o.nb[[i]]))
      for (j in 1:length(o.nb[[i]])){
      vec[j,]<-  CRIME[o.nb[[i]][j]]*a.lw$weights[[i]][j]
      wy[i,] <- sum(vec)
      }
      }

# Otra forma m�s simple ser�a:
W <- as.matrix(as_dgRMatrix_listw(a.lw))
Y <- CRIME
WY <- W%*%Y

lag.listw(a.lw, CRIME)                      # Directamente en el paquete "spdep"
 
# Luego se define el N y el S, que para este ejercicio son ambos iguales a 49, y obtenemos el �ndice I de Moran aplicando la Ecuaci�n (1):
N <- length(CRIME)
S <- sum(unlist(a.lw$weights))
I <- (N/S)*((sum(x.centrado * wy))/scx)
I <- as.vector(scale(Y))%*%W%*%as.vector(scale(Y))/as.vector(scale(Y))%*%as.vector(scale(Y))

# El estad�stico de Moran obtenido para la variable de inter�s indica que existe una correlaci�n espacial de 0.5 
# Luego se calcula la kurtosis de la muestra:
K <- (length(x) * sum(x.centrado^4))/(scx^2)


# Directamente el �ndice de Moran y la Kurtosis mediante la funci�n moran() son:
#x1<-(x-mean(x))/sd(x)                           # Estandarizando, para obtener el mismo gr�fico de GeoDA
moran(x, listw=a.lw, n=length(o.nb), S0=S, zero.policy=F)       # Cambiar por "x1", para la estandarizaci�n
set.seed(123)
moran.test(x, nb2listw(o.nb, style="W"), zero.policy=F)
Z.I <- (0.485770914-(-0.020833333))/sqrt(0.008991121)

wc <- spweights.constants(a.lw, zero.policy = F, adjust.n = T)
W <- as.matrix(as_dgRMatrix_listw(a.lw))
S0 <- sum(W)
S02 <- S0^2
nc <- apply(W,2,sum)
nr <- apply(W,1,sum)
S2 <- sum((nc+nr)^2)
n <- sum(nc)
nn <- n*n
sum(W^2)
W1 <- W
W2 <- W
W1[lower.tri(W1, diag=T)] <- 0
W2[upper.tri(W2, diag=T)] <- 0
W3 <- (t(W1)+W2)^2
#S1 <- sum(((t(W)+W)^2)/2)
S1 <- sum(W3)
a.lw
EI <- -1/(n-1) 
# sin randomisation
VI <- (nn*S1-n*S2+3*S02)/(S02*(nn-1))- (EI)^2
# con randomisation
VI <- (n*((n^2-3*n+3)*S1-n*S2+3*S02)-K*(S1*(nn-n)-2*n*S2+6*S02))/((n-1)*(n-2)*(n-3)*S02)-1/((n-1)^2)
ZI <- (I-EI)/sqrt(VI)

set.seed(123)
geary.test(x, nb2listw(o.nb, style="W"), zero.policy=F)
(0.547803377-1)/sqrt(0.009804108)                               # Autocorrelaci�n espacial, concentraci�n de valores similares
globalG.test(x, nb2listw(o.nb, style="B"))

# Gr�fico de dispersi�n del �ndice de Moran
mp<-moran.plot(x, a.lw, main="Gr�fico de Dispersi�n de Moran", ylim=c(-2,2))      # Cambiar por "x1", para la estandarizaci�n
mp$infmat
set.seed(127)
mc<-moran.mc(x, a.lw, nsim=1000, zero.policy=F)     # Se rechaza Ho:rho=0, es decir hay autocorrelaci�n espacial.
mc
plot(mc)                                            # Funci�n de densidad de las simulaciones de Monte Carlo.
gt <- geary.mc(x, a.lw, nsim=1000, zero.policy=F)
gt                                                              
plot(gt)

moran.bi(col.poly@data$CRIME,col.poly@data$INC,a.lw,zero.policy =T)
set.seed(123)
moranbi.test(X=col.poly@data$CRIME,Y=col.poly@data$INC,a.lw,zero.policy =T,adjust.n=TRUE,N=999,graph=T,print.results=T)

a.lwq1 <- nb2listw(col.lags[[1]], style="W")
a.lwq2 <- nb2listw(col.lags[[2]], style="W")
a.lwq3 <- nb2listw(col.lags[[3]], style="W")
a.lwq4 <- nb2listw(col.lags[[4]], style="W")
a.lwq5 <- nb2listw(col.lags[[5]], style="W")
a.lwq6 <- nb2listw(col.lags[[6]], style="W", zero.policy =T)

#============================================================
# Criterios basados en gr�ficas
#============================================================
op=par(mfrow=c(2,2))
trinb=tri2nb(coords)
delaunay <-nb2listw(trinb, style="W")
plot(col.poly,border="gray")
plot(trinb,coords,add=T,col="blue")
title(main="Triangulaci�n Delaunay")
soinb=graph2nb(soi.graph(trinb,coords))
esf.influencia <-nb2listw(soinb, style="W")
plot(col.poly,border="gray")
plot(soinb,coords,add=T,col="green")
title(main="Esfera de influencia")
gabrielnb=graph2nb(gabrielneigh(coords),sym=TRUE)
gabriel <-nb2listw(gabrielnb, style="W")
plot(col.poly,border="gray")
plot(gabrielnb,coords,add=T,col="red")
title(main="Gr�fica de Gabriel")
relativenb=graph2nb(relativeneigh(coords),sym=TRUE)
vec.relative <-nb2listw(relativenb, style="W")
plot(col.poly,border="gray")
plot(relativenb,coords,add=T,col="orange")
title(main="Vecinos relativos")
par(op)


#==========================================================================
############### SELECCI�N DE MATRICES DE VECINDAD POR (PCNM) ##############
#==========================================================================
#principal coordinates of neighbour matrices (PCNM, Borcard and Legendre (2002))
col_nbq1 <- poly2nb(col.poly)               # Efecto Reina
colt.lags <- nblag(col_nbq1, 10)
q.lw1  <- nb2listw(colt.lags[[1]], style="W")
q.lw2  <- nb2listw(colt.lags[[2]], style="W")

col_nbr1 <- poly2nb(col.poly,queen=FALSE)
colr.lags <- nblag(col_nbr1, 10)
r.lw1  <- nb2listw(colr.lags[[1]], style="W")
r.lw2  <- nb2listw(colr.lags[[2]], style="W")
r.lw3  <- nb2listw(colr.lags[[3]], style="W")
r.lw4  <- nb2listw(colr.lags[[4]], style="W")
r.lw5  <- nb2listw(colr.lags[[5]], style="W")
r.lw6  <- nb2listw(colr.lags[[6]], style="W", zero.policy =T)

summary(test.W(col.poly@data$CRIME,colr.lags[[2]],coords))
summary(test.W(col.poly@data$CRIME,col_nbq1))
summary(test.W(col.poly@data$CRIME,col.knlags[[5]]))
summary(test.W(col.poly@data$CRIME,trinb))
summary(test.W(col.poly@data$CRIME,soinb))
summary(test.W(col.poly@data$CRIME,gabrielnb))
summary(test.W(col.poly@data$CRIME,relativenb))
summary(test.W(col.poly@data$CRIME,col.kn7))
nb2listw(col.knlags[[6]], style="B", glist = NULL)

#=============================================================
############### SELECCI�N DE LA MATRIZ DE PESOS ##############
#=============================================================

Pesos.list<-list(reina1=a.lwq1,reina2=a.lwq2,reina3=a.lwq3,reina4=a.lwq4,reina5=a.lwq5,reina6=print.listw(a.lwq6,zero.policy=TRUE),
                             torre1=r.lw1,torre2=r.lw2,torre3=r.lw3,torre4=r.lw4,torre5=r.lw5,torre6=print.listw(r.lw6,zero.policy=TRUE), 
                 kvecinos1=col.k1,kvecinos2=col.k2,kvecinos3=col.k3,kvecinos4=col.k4,kvecinos5=col.k5,gabriel=gabriel,
                 delaunay=delaunay,esfera.inf=esf.influencia,vec.relativos=vec.relative)

class(Pesos.list)
# Bauman et al. (2018) observo que optimizar la elecci�n de la Matriz de Pesos Espaciales (SWM, por sus siglas 
# en ingl�s) conduc�a a tasas de error de tipo I infladas si no se aplicaba un control expl�cito del n�mero de 
# SWM probadas. Por lo tanto, se debe aplicar una correcci�n de Sid�k (1967) para m�ltiples pruebas al P_value 
# de la prueba global de cada SWM (es decir, el modelo que integra todo el conjunto de predictores espaciales). 
# La correcci�n de Sidak se calcula como Pcorregido = 1 - (1 - P)^(1/n) donde n es el n�mero de pruebas realizadas, 
# P es el P_value observado y Pcorrected es el nuevo valor p despu�s de la correcci�n. El valor p se calcula 
# primero usando n permutaciones y luego se corrige de acuerdo con el n�mero total de SWM probadas.
nbw <- length(Pesos.list)
# al 5% de significanc�a
1 - (1 - 0.05)^(nbw)
# la optimizaci�n para la selecci�n se realiza maximizando el R cuadrado ajustado (R2 Adjust) o minimizando 
# la autocorrelaci�n espacial residual.
W_sel <- listw.select(col.poly$CRIME, Pesos.list, MEM.autocor = "all", p.adjust = TRUE, nperm = 50)
W_sel$candidates
W_sel$best.id

# Para copiar y editar en Latex
xtable(W_sel$candidates, digits = 6)
W_sel$best$MEM.select
# El criterio que m�s valor tiene para la selecci�n definitiva de la matriz de contig�idad
# es el R cuadrado ajustado (RAdjust 2), que minimiza la autocorrelaci�n. Por esto la matriz
# elegida ser� K vecinos con 5 vecinos m�s cercanos.
Best.SWM <- Pesos.list[W_sel$best.id]             # Para el ejemplo SWM es col.k2
detach(package:adespatial)

# I Moran bivariate
set.seed(123)
MBCI1 <- moranbi.test(col.poly@data$CRIME,col.poly@data$INC,col.k5,N=999,zero.policy =T,graph=T) # Rechazo Ho
MBCI2 <- moranbi.test(col.poly@data$CRIME,col.poly@data$INC,a.lwq1,N=999,zero.policy =T,graph=T) # Rechazo Ho
MBCI3 <- moranbi.test(col.poly@data$CRIME,col.poly@data$INC,a.lwq2,N=999,zero.policy =T,graph=T) 
# No Rechazo Ho
moran.bi(col.poly@data$CRIME,col.poly@data$INC,a.lwq6,zero.policy =T)
moranbi.plot(col.poly@data$CRIME,col.poly@data$INC,quiet =F,zero.policy =F,listw=a.lwq1)

set.seed(123)
MBCH1 <- moranbi.test(col.poly@data$CRIME,col.poly@data$HOVAL,a.lwq1,N=999,zero.policy =T,graph=T) 
# Rechazo Ho
MBCH2 <- moranbi.test(col.poly@data$CRIME,col.poly@data$HOVAL,a.lwq2,N=999,zero.policy =T,graph=T) 
# No Rechazo Ho
MBCH3 <- moranbi.test(col.poly@data$CRIME,col.poly@data$HOVAL,a.lwq3,N=999,zero.policy =T,graph=T) 
# No Rechazo Ho


# LISA Cluster Map
moran.cluster(col.poly@data$CRIME, a.lwq1, zero.policy = T, col.poly, significant=T)

LMCI <- localmoran.bi(col.poly@data$CRIME, col.poly@data$INC, a.lwq1, zero.policy =T)
LMCH <- localmoran.bi(col.poly@data$CRIME, col.poly@data$HOVAL, a.lwq1, zero.policy =T)

# Getis Cluster Map
getis.cluster(col.poly@data$CRIME, a.lwq1, zero.policy = T, col.poly, significant=T)

###############################################
#  Mapeando outliers e influyentes locales   ##
###############################################

x1<-(x-mean(x))/sd(x)
mp<-moran.plot(x1, a.lwq1, main="Gr�fico de Dispersi�n de Moran")      
# Cambiar por "x1", para la estandarizaci�n

infl1 <- apply(mp$is.inf, 1, any)
lhx1 <- cut(col.poly@data$CRIME, breaks=c(min(col.poly@data$CRIME), mean(col.poly@data$CRIME), 
                          max(col.poly@data$CRIME)), labels=c("L", "H"), include.lowest=TRUE)
wx1 <- lag(a.lwq1, columbus$CRIME)
lhwx1 <- cut(wx1, breaks=c(min(wx1), mean(wx1), max(wx1)), labels=c("L", "H"), include.lowest=TRUE)
lhlh1 <- interaction(lhx1, lhwx1, infl1, drop=TRUE)
cols1 <- rep(1, length(lhlh1))
cols1[lhlh1 == "L.L.TRUE"] <- 2
cols1[lhlh1 == "L.H.TRUE"] <- 3
cols1[lhlh1 == "H.L.TRUE"] <- 4
cols1[lhlh1 == "H.H.TRUE"] <- 5
plot(columbus, col=colcode[cols1])          # gray.colors(4, 0.95, 0.55, 2.2)[cols])
legend(locator(1), legend=c("None", "LL", "LH", "HL", "HH"), fill=attr(colcode, "palette"), bty="n", cex=1, y.intersp=0.8)

##################################################################################################


#############################
# Ajuste Modelo Cl�sico
#############################

# columbus.data <- as.data.frame(col.poly)
lm_fit <- lm(CRIME ~ INC + HOVAL, data=columbus)
summary(lm_fit)
columbus$lm_res <- residuals(lm_fit)
# Mapa residuales: Opci�n 1
spplot(columbus["lm_res"], col.regions = rev(terrain.colors(20)))

# Mapa residuales: Opci�n 2
pal2 <- colorRampPalette(c("red3", "wheat1", "blue3"))
spplot(columbus,"lm_res", col.regions = pal2(20))

# Mapa residuales: Opci�n 3
my.palette <- brewer.pal(n = 9, name = "YlOrRd")
spplot(columbus, "lm_res", col.regions = my.palette, cuts = 8, col = "transparent")

# Mapa residuales: Opci�n 4
breaks.ci <- classIntervals(columbus$lm_res, n = 9, style = "quantile", intervalClosure = "right")$brks
# se amplian los limites inferior y superior en: .Machine$double.eps*5000000000000,  con el fin que un poligono no quede siempre en blanco.
breaks.ci[1] <- breaks.ci[1] - .Machine$double.eps*5000000000000
breaks.ci[length(breaks.ci)] <- breaks.ci[length(breaks.ci)] + .Machine$double.eps*5000000000000
spplot(columbus, "lm_res", col = "transparent", col.regions = my.palette,  at = breaks.ci)

lm.morantest(lm_fit, a.lwq1)
Z.Ir <- (0.222109407-(-0.033418335))/sqrt(0.008099305)                                                  
# Por lo tanto, los residuos estan autocorrelacionados

#set.seed(123)
#moran.test(residuals(lm_fit), nb2listw(o.nb, style="W"), zero.policy=F)

# dwtest(CRIME ~ INC + HOVAL, data=columbus, alternative = "two.sided")         
# Los residuos no estan autocorrelacionados prueba no espacial
# acf(columbus$lm_res)
# pacf(columbus$lm_res)

LM_res <- lm.LMtests(lm_fit, a.lw, test = "all")
t(sapply(LM_res, function(x) unlist(x[1:3])))

u <- resid(lm_fit)
sigma2 <- (t(u) %*% u)/nrow(W)
T1 <- sum(diag(t(W)%*%W+W%*%W))
LM.ERR <- (((t(u)%*%W%*%u)/sigma2)^2)/T1
p.value <- 1-pchisq(LM.ERR,1)


# Validaci�n supuestos
bptest(lm_fit)
resettest(lm_fit)
raintest(lm_fit)
shapiro.test(residuals(lm_fit))
library(car)
vif(lm_fit)

library(RcmdrPlugin.epack)
# Transformaci�n Box-Cox
BC <- boxcox(CRIME ~ INC+HOVAL, data = columbus, lambda = seq(-3, 3, len = 20))
bc2(columbus$CRIME)
bc2(columbus$HOVAL)
bc2(columbus$INC)
detach(package:RcmdrPlugin.epack)

##################################################################################
# (No es importante correr esto)
# LISA residuales  
lm_sad <- localmoran.sad(lm_fit, nb = col.nb.0.all)                    # local de Moran
lm_sad_df <- as.data.frame(lm_sad)
lm_rand <- localmoran(residuals(lm_fit), a.lw)
plot(lm_sad_df[,5], lm_rand[,5] , type="p")

# Getis-Ord global G statistic
globalG.test(x, listw=a.lwq1)                                            

# Local de Getis y Ord
lg <- localG(x, listw=a.lw)                                                     
lg
plot(lg)
##################################################################################


###########################################
# 3. AJUSTE DEL MODELO SAR Y VALIDACI�N
###########################################

# primero se seleccionan las variables que son utiles para el modelo:
datos.col <- as.data.frame(col.poly)
d.col <- as.data.frame(cbind(datos.col$POLYID, datos.col$HOVAL, datos.col$INC, datos.col$CRIME))
colnames(d.col) <- c("POLYID", "HOVAL", "INC", "CRIME")
attach(d.col)

# La variable dependiente de nuestro modelo ser� la variable CRIME, y las variables 
# explicativas son INC y HOVAL. Se procede a separarlas:
Y <- CRIME
INC <- datos.col$INC
HOVAL <- datos.col$HOVAL
X <- cbind(INC, HOVAL)

# Adicionamos una columna de unos (1) al conjunto de las variables independientes en el objeto x 
# previamente creado. Con base en los objetos x y y se har� posteriormente la estimaci�n del 
# par�metro (rho) del modelo. Observamos las primeras 5 filas del objeto x:

X <- cbind(1, X)
X[1:5,]

cardnb <- card(a.lw$neighbours)                   # cuenta los vecinos de cada poligono
str(sd1 <- 1/sqrt(cardnb))

glist <- vector(mode = "list", length = length(a.lw$neighbours))
 for (i in 1:length(a.lw$neighbours)) glist[[i]] <- cardnb[i]*a.lw$weights[[i]]
 for (i in 1:length(a.lw$neighbours)) {
     inb <- a.lw$neighbours[[i]]
     icd <- cardnb[i]
     if (icd > 0) {
 for (j in 1:icd) {
glist[[i]][j] <- sd1[i] * glist[[i]][j] * sd1[inb[j]]
}
}
}

W <- as_dgRMatrix_listw(similar.listw(a.lw))
str(W)
W@x
W <- as.matrix(W)

I <- as_dsCMatrix_I(nrow(X))
# I <- as_dsCMatrix_I(dim(W)[1])
I@x                                              # Slot "x":
I@i                                              # Slot "i":
I@p                                              # Slot "p":
I@Dim                                            # Slot "Dim":
I <- as.matrix(I)

lm.null <- lm(Y ~ X - 1)                         # Modelo sin intercepto
lm.w <- lm.fit(X, WY)
e.null <- lm.null$residuals
e.w <- lm.w$residuals
e.a <- t(e.null) %*% e.null
e.b <- t(e.w) %*% e.null
e.c <- t(e.w) %*% e.w


sar.lag.mix.f.sM <- function(rho, W, I, e.a, e.b, e.c, n, tmpmax){
 	SSE <- e.a - 2*rho*e.b + rho*rho*e.c
 	s2 <- SSE/n
 	Jacobian <- log(det(chol((I - rho * W), LINPACK=TRUE))^2)
		# tmpmax=(3*((sum(card(a.lw$neighbours))) + n))))^2)
 	ret <- (Jacobian - ((n/2)*log(2*pi)) - (n/2)*log(s2) - 	(1/(2*s2))*SSE)
 	ret
 }

opt <- optimize(sar.lag.mix.f.sM, interval = c(-1,1), maximum =  TRUE, 
tol = .Machine$double.eps^0.5, W = W, I = I, e.a = e.a, e.b =  e.b, e.c = e.c, n = nrow(X))
 
rho <- opt$maximum

lm.lag <- lm((Y - rho * WY) ~ X - 1)
summary(lm.lag)

#lm.lag1 <- nls(Y ~ X%*%c(a,b,c)+rho*WY, start=c(a=0,b=1.5,c=1,rho=0),lower=c(a=-Inf,b=-Inf,c=-Inf,
#              rho=-1),upper=c(a=Inf,b=Inf,c=Inf,rho=1))
#summary(lm.lag1)

####################################################################################
# Lo mismo pero simplificando en R, a partir de la instrucci�n del paquete "spdep":
####################################################################################

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

#Nagelkerke NJD (1991) A note on a general definition of the coefficient of determination. Biometrika 78: 691�C692. 
NK <- function(obj, y) { 
     n <- length(obj$residuals) 
     nullLL <- logLik(lm(y ~ 1)) 
     c(1 - exp(-(2/n)*(logLik(obj) - nullLL))) 
      } 
NK(col.lag.sm,CRIME)

library(car)
bc <- boxCox(CRIME ~ INC + HOVAL, data = columbus, lambda = seq(-1, 2, length = 20))
bc$x[which.max(bc$y)]

# Modelo Spatial Error:
col.error.sm <- errorsarlm(CRIME ~ INC + HOVAL + X + Y, data=as.data.frame(col.poly), listw=col.k5)
summary(col.error.sm, Nagelkerke=T)
NK(col.error.sm,CRIME)
predict.sarlm(col.error.sm)
AIC(col.error.sm)
deviance.sarlm(col.error.sm)
residuals.sarlm(col.error.sm)
coef.sarlm(col.error.sm)
fitted.sarlm(col.error.sm)
bptest.sarlm(col.error.sm)
hetero.plot(col.error.sm)

# Comparaci�n
anova(col.lag.sm,col.error.sm)            

# Verificaci�n de Supuestos: "Spatial Lag"

# Spatial Lag:
shapiro.test(col.lag.sm$residuals)
# ks.test(col.lag.sm$residuals, alternative = "two.sided", pnorm) 
plot(col.lag.sm$residuals)                                  # Gr�fico de los residuales
dens.sl <- density(col.lag.sm$residuals)                    # Gr�fico de densidad
plot(dens.sl, main="Gr�fico de la funci�n de densidad\n de los residuos (Spatial Lag)") 
bptest.sarlm(col.lag.sm)


# Spatial Error:
shapiro.test(col.error.sm$residuals)
# ks.test(col.error.sm$residuals, alternative = "two.sided", pnorm) 
plot(col.error.sm$residuals)                                # Gr�fico de los residuales
dens.se <- density(col.error.sm$residuals)                  # Gr�fico de densidad
plot(dens.se, main="Gr�fico de la funci�n de densidad\n de los residuos (Spatial Error)") 
 
# SARAR
col.sarar <- sacsarlm(CRIME ~ INC + HOVAL, data=as.data.frame(col.poly), listw=col.k5)    # CRIME ~1
summary(col.sarar, Nagelkerke=T)
predict.sarlm(col.lag.sm)
AIC(col.sarar)
deviance.sarlm(col.sarar)
residuals.sarlm(col.sarar)
coef.sarlm(col.sarar)
fitted.sarlm(col.sarar)
bptest.sarlm(col.sarar)
hetero.plot(col.lag.sm)
shapiro.test(col.sarar$residuals)

# spautolm:  "SAR" or "CAR" for simultaneous or conditional autoregressions, "SMA" for spatial moving average
col.autSAR <- spautolm(CRIME~INC + HOVAL + X + Y, data=as.data.frame(col.poly), listw=col.k5, family="SAR") 
summary(col.autSAR, Nagelkerke=T)
hetero.plot(col.autSAR)

col.autCAR <- spautolm(CRIME~INC + HOVAL, data=as.data.frame(col.poly), listw=col.k5, family="CAR") 
summary(col.autCAR, Nagelkerke=T)

col.autSMA <- spautolm(CRIME~INC + HOVAL + X + Y, data=as.data.frame(col.poly), listw=col.k5, family="SMA") 
summary(col.autSMA, Nagelkerke=T)
                          
# Durbin spatial 
columlagsd<-lagsarlm(CRIME~INC + HOVAL + X + Y, data=as.data.frame(col.poly), listw=a.lw, type="mixed") 
summary(columlagsd, Nagelkerke=T)
predict.sarlm(columlagsd)
AIC(columlagsd)
deviance.sarlm(columlagsd)
residuals.sarlm(columlagsd)
coef.sarlm(columlagsd)
fitted.sarlm(columlagsd)
bptest.sarlm(columlagsd)
shapiro.test(columlagsd$residuals)

# Comparaci�n (Es mejor columlagsd)
anova(col.lag.sm,columlagsd)             
anova(columlagsd,col.error.sm)
anova(columlagsd,col.sarar)            
anova(col.lag.sm,col.error.sm)
            

# Se rezagan las variables, para ello se copia la informaci�n en col.poly1
col.poly1 <- as.data.frame(col.poly)
col.poly1$WX <- lag.listw(a.lw, col.poly1$X)
col.poly1$WY <- lag.listw(a.lw, col.poly1$Y)
col.poly1$WINC <- lag.listw(a.lw, col.poly1$INC)
col.poly1$WHOVAL <- lag.listw(a.lw, col.poly1$HOVAL)
a.lw1 <- nb2listw(col.lags[[1]], style="W")
a.lw2 <- nb2listw(col.lags[[2]], style="W")
a.lw3 <- nb2listw(col.lags[[3]], style="W")
a.lw4 <- nb2listw(col.lags[[4]], style="W")
a.lw5 <- nb2listw(col.lags[[5]], style="W")
a.lw6 <- nb2listw(col.lags[[6]], style="W",zero.policy =T)
a.lw7 <- nb2listw(col.lags[[7]], style="W",zero.policy =T)
col.poly1$W1CRIME <- lag.listw(a.lw1,col.poly1$CRIME)
col.poly1$W2CRIME <- lag.listw(a.lw2,col.poly1$CRIME)
col.poly1$W3CRIME <- lag.listw(a.lw3,col.poly1$CRIME)
col.poly1$W4CRIME <- lag.listw(a.lw4,col.poly1$CRIME)
col.poly1$W5CRIME <- lag.listw(a.lw5,col.poly1$CRIME)
col.poly1$W6CRIME <- lag.listw(a.lw6,col.poly1$CRIME, zero.policy =T)

# Las variables regresoras de abajo son s�lo las significativas estadisticamente   (col.k2)
col.lagx <- lagsarlm(CRIME ~ INC + HOVAL + X + WX, data=col.poly1, listw=a.lw) 
summary(col.lagx, Nagelkerke=T)
# Se prueban posibles transformaciones de acuerdo a los resultados de la transformaci�n Box-Cox
col.lagxt <- lagsarlm(CRIME~ log(INC) + I(1/sqrt(HOVAL)) + X + WX, data=col.poly1, listw=a.lw) 
summary(col.lagxt, Nagelkerke=T, correlation=TRUE)
bptest.sarlm(col.lagx)
hetero.plot(col.lagx)
shapiro.test(col.lagx$residuals)
predict.sarlm(col.lagx)
AIC(col.lagx)
deviance.sarlm(col.lagx)
residuals.sarlm(col.lagx)
coef.sarlm(col.lagx)
fitted.sarlm(col.lagx)
# Calculo pseudoR2:
pseudoR2 <- cor(fitted.sarlm(col.lagx),col.poly1$CRIME)^2

# Comparaci�n (Es mejor col.lagx "Durbin Espacial")
anova(col.lag.sm,col.lagx)            
anova(col.lagx,col.lagxt)

# Validaci�n supuestos Durbin Espacial:
bptest.sarlm(col.lagx)
shapiro.test(col.lagx$residuals)
moran.test(col.lagx$residuals, nb2listw(o.nb, style="W"),zero.policy=T,alternative="two.sided")

plot(col.lagx$residuals)                                  # Gr�fico de los residuales
dens.sd <- density(col.lagx$residuals)                    # Gr�fico de densidad
plot(dens.sd, main="Gr�fico de la funci�n de densidad\n de los residuos (Spatial Durbin)") 

# Modelo glm tasa criminalidad
col.glm<-glm(CRIME~W1CRIME + W5CRIME + INC + HOVAL + WHOVAL + X +Y + WX + WY, data=col.poly1, family=poisson(link="log")) 
summary(col.glm)

residuals(col.glm)
hetero.plot(col.glm)
moran.test(residuals(col.glm), nb2listw(o.nb, style="W"),zero.policy=T,alternative="two.sided")
moran.plot(residuals(col.glm), nb2listw(o.nb, style="W"), xlim=c(-6.5,6.5),ylim=c(-2,2))
# Calculo pseudoR2:
pseudoR2.glm <- cor(exp(predict(col.glm)),col.poly$CRIME)^2

bptest(col.glm)

# Supuesto de normalidad:
residuos.deviance <- sum(residuals(col.glm, type = "deviance")^2)
col.glm$deviance
res.dev.std <- rstandard(col.glm, type = "deviance")
hist(res.dev.std[-c(7,20)], breaks = 40)
abline(v = quantile(res.dev.std[-c(7,20)], probs = c(0.05/2, 1 - 0.05/2)), lty = 2)
shapiro.test(res.dev.std[-c(7,20)])

# Comportamiento residuales
par(mfrow = c(2, 2))
plot(col.glm)

# Heterocedasticidad con sphet:  PRUEBA
col.het <- gstslshet(CRIME~W1CRIME + W4CRIME + W6CRIME + INC + X, data=col.poly1, listw=nb2listw(o.nb, style="W"))
summary(col.het)

moran.test(residuals(col.het), nb2listw(o.nb, style="W"),zero.policy=T,alternative="two.sided")
moran.plot(residuals(col.het), nb2listw(o.nb, style="W"))

# Calculo pseudoR2:
pseudoR2.het <- cor(col.het$yhat,col.poly$CRIME)^2

data(coldis)
col.stls <- stslshac(CRIME~W1CRIME + W4CRIME + INC + X, data=col.poly1,listw=nb2listw(o.nb, style="W"),
                     distance=coldis, type= 'Triangular')
summary(col.stls)
# Calculo pseudoR2:
pseudoR2.stls <- cor(col.stls$yhat,col.poly$CRIME)^2

# Modelo Mixto Spatial Error (Este modelo funciona bien) SDEM:

col.errorM.sm <- errorsarlm(CRIME ~ INC + HOVAL + X + WX, data=col.poly1, listw=a.lw)
summary(col.errorM.sm, correlation=TRUE, Nagelkerke=T)
NK(col.errorM.sm,CRIME)
predict.sarlm(col.errorM.sm)
AIC(col.errorM.sm)
deviance.sarlm(col.errorM.sm)
residuals.sarlm(col.errorM.sm)
coef.sarlm(col.errorM.sm)
fitted.sarlm(col.errorM.sm)
bptest.sarlm(col.errorM.sm)
hetero.plot(col.errorM.sm)
moran.test(residuals(col.errorM.sm), nb2listw(o.nb, style="W"),zero.policy=T,alternative="two.sided")

# Calculo pseudoR2:
pseudoR2.errorM.sm <- cor(col.errorM.sm$fitted.values,col.poly$CRIME)^2

# Modelo Mixto Spatial Error Heteroscedastico (Chasco Yrigoyen, Anselin):
modelB <- lm(CRIME ~ INC + HOVAL, data=col.poly1)
CME <- anova(modelB)["Mean Sq"][3,]
res.est2 <- residuals(modelB)^2/CME
model.varE <- lm(res.est2 ~ I(HOVAL^2) , data=as.data.frame(col.poly))    # I(INC^2) no se considera ya uqe genera pesos negativos
summary(model.varE)
col.errorH.sm <- errorsarlm(CRIME ~ INC + HOVAL, data=col.poly1, listw=a.lw, weights=model.varE$fitted.values)
summary(col.errorH.sm, correlation=TRUE, Nagelkerke=T)
NK(col.errorH.sm,CRIME)
predict.sarlm(col.errorH.sm)
AIC(col.errorH.sm)
deviance.sarlm(col.errorH.sm)
residuals.sarlm(col.errorH.sm)
coef.sarlm(col.errorH.sm)
fitted.sarlm(col.errorH.sm)
bptest.sarlm(col.errorH.sm)                                               # A�n es heteroscedastico
hetero.plot(col.errorH.sm)
moran.test(residuals(col.errorH.sm), nb2listw(o.nb, style="W"),zero.policy=T,alternative="two.sided")           # Y libre de autocorrelaci�n

# Calculo pseudoR2:
pseudoR2.errorH.sm <- cor(col.errorH.sm$fitted.values,col.poly$CRIME)^2

############################
##  An�lisis de Impactos  ##
############################

col.new <- col.poly
col.new@data <- col.poly1

# Cambiando la tasa de criminalidad en u
col.new@data[col.new@data$NEIG == "38","CRIME"] <- 10.0

# Los valores de las predicciones originales
orig.pred <- as.data.frame(predict(col.lagx))

# Los valores predichos con la nueva tasa de criminalidad en el barrio 10
col_nbq1
new.pred <- as.data.frame(predict(col.lagx, newdata = col.new, listw = nb2listw(col_nbq1, style="W")))

# Las diferencias entre las predicciones
effect.10 <- new.pred$fit - orig.pred$fit                  
el <- data.frame(name = col.new$NEIG, dif_pred_CRIME = effect.10)
col.new$ef10 <- el$dif_pred_CRIME

# Ordenando los barrios por el valor absoluto del cambio en la predicci�n del CRIME
el <- el[rev(order(abs(el$dif_pred_CRIME))), ]
el[1:10, ]  #muestra los 10 primeros barrios

# Mapear estos cambios es tambien importante:

breaks <- c(min(col.new$ef10), -0.05, 0.05, max(col.new$ef10))
labels <- c("Efecto negativo (< -.05)", "Sin efecto (-.05 a .05)", 
    "Efecto positivo (> .05)")

# faltaba all.inside =T para evitar un poligono en blanco asociado al valor m�ximo  
np <- findInterval(col.new$ef10, breaks,all.inside =T) 
colors <- c("red", "yellow", "blue")

# Dibujando el mapa
plot(col.new, col = colors[np])
mtext("Efectos de un cambio en el barrio 10, (fijando CRIME=10)\n sobre los valores predichos en un modelo Durbin Espacial", 
    side = 3, line = 1)
legend("topleft", legend = labels, fill = colors, bty = "n")

# Tambi�n podr�amos mapear la magnitud de los cambios causados a ra�z de la disminuci�n de la criminalidad en el barrio 10.

pal5 <- brewer.pal(6, "Spectral")
cats5 <- classIntervals(col.new$ef10, n = 5, style = "jenks")
colors5 <- findColours(cats5, pal5)
plot(col.new, col = colors5)
legend("topleft", legend = round(cats5$brks, 2), fill = pal5, bty = "n")
mtext("Efectos de un cambio en el barrio 10, (fijando CRIME = 10)\n sobre los valores de las 
      predicciones en un modelo Durbin Espacial", side = 3, line = 1)
   
# While these maps and the effects for individual counties may be useful they only tell us about a 
# change in one place. The impacts() function gives us something like OLS regression coefficients for
# a spatial lag model. The logic of the impacts() function is similar to the code above, it tells you
# the direct (local), indirect (spill-over), and total effect of a unit change in each of the 
# predictor variables. The changes reported by impacts are the global average impact:


W <- as(as_dgRMatrix_listw(a.lw), "CsparseMatrix")
trMatc <- trW(W, type="mult")
trMC <- trW(W, type="MC")

impacts(col.lagx, listw=a.lw)

# Impacto Direct INC:
SrW.I <- solve(diag(49)-as.numeric(col.lagx["rho"])*as.matrix(W))%*%diag(49)*as.numeric(coefficients(col.lagx)["INC"])
sum(diag(SrW.I))/49

# Impacto Direct HOVAL:
SrW.H <- solve(diag(49)-as.numeric(col.lagx["rho"])*as.matrix(W))%*%diag(49)*as.numeric(coefficients(col.lagx)["HOVAL"])
sum(diag(SrW.H))/49

# Impacto Direct X:
SrW.X <- solve(diag(49)-as.numeric(col.lagx["rho"])*as.matrix(W))%*%diag(49)*as.numeric(coefficients(col.lagx)["X"])
sum(diag(SrW.X))/49

# Impacto Direct Y:
SrW.Y <- solve(diag(49)-as.numeric(col.lagx["rho"])*as.matrix(W))%*%diag(49)*as.numeric(coefficients(col.lagx)["Y"])
sum(diag(SrW.Y))/49

# Impacto Direct WX:
SrW.WX <- solve(diag(49)-as.numeric(col.lagx["rho"])*as.matrix(W))%*%diag(49)*as.numeric(coefficients(col.lagx)["WX"])
sum(diag(SrW.WX))/49

# Impacto Direct Y:
SrW.WY <- solve(diag(49)-as.numeric(col.lagx["rho"])*as.matrix(W))%*%diag(49)*as.numeric(coefficients(col.lagx)["WY"])
sum(diag(SrW.WY))/49

# Impacto Indirect INC:
(sum(SrW.I) - sum(diag(SrW.I)))/49

# Impacto Indirect INC:
(sum(SrW.H) - sum(diag(SrW.H)))/49

# Impacto Total INC:
rep(1,49)%*%(SrW.I)%*%matrix(rep(1,49),ncol=1)/49              # otra forma mas simple es: sum(SrW.I)/49

# Impacto Total HOVAL:
rep(1,49)%*%(SrW.H)%*%matrix(rep(1,49),ncol=1)/49              # otra forma mas simple es: sum(SrW.H)/49

# Nuestro modelo que trata las tasas de criminalidad CRIME (Robo de viviendas y veh�culos por cada 1000 hogares): en terminos de los impactos 
# nos dice que un aumento del 100% en HOVAL (valor de la vivienda en miles de dolares) conduce a una caida del  34.9% (en promedio) en la tasa 
# de CRIME (Robo de viviendas y veh�culos por cada 1000 hogares).

impacts(col.lagx, tr=trMatc)
impacts(col.lagx, tr=trMC)
               
# Procedimiento para evaluar los niveles de significancia
lobj <- lagsarlm(CRIME ~ INC + HOVAL + X +Y + WX + WY, col.poly1, a.lw)
summary(lobj)

lobj1 <- stsls(CRIME ~ INC + HOVAL + X +Y + WX + WY, col.poly1, a.lw)
loobj1 <- impacts(lobj1, tr=trMatc, R=200)
summary(loobj1, zstats=TRUE, short=TRUE)
lobj1r <- stsls(CRIME ~ INC + HOVAL + X +Y + WX + WY, col.poly1, a.lw, robust=TRUE)
loobj1r <- impacts(lobj1r, tr=trMatc, R=200)
summary(loobj1r, zstats=TRUE, short=TRUE)
lobjIQ5 <- impacts(lobj, tr=trMatc, R=200, Q=5)
summary(lobjIQ5, zstats=TRUE, short=TRUE)

summary(lobjIQ5, zstats=TRUE, short=TRUE, reportQ=TRUE)
impacts(mobj, listw=a.lw)
impacts(mobj, tr=trMatc)
impacts(mobj, tr=trMC)
summary(impacts(mobj, tr=trMatc, R=200), zstats=TRUE)

# http://rstudio-pubs-static.s3.amazonaws.com/5027_52298866e7924b18b54e5c9a0a21b450.html
# Interpreting Spatial Lag Models

########################################################################
##########       Regresi�n Geogr�ficamente Ponderada       #############
########################################################################

# https://gis.stackexchange.com/questions/241127/how-to-plot-output-from-gwr-in-r
# http://geokitchen.blogspot.com.co/2012/09/r-geographically-weighted-regression.html

adapt <- gwr.sel(CRIME~ INC + HOVAL + X + Y + WX + WY, data=col.poly1, coords=cbind(col.poly1$X,col.poly1$Y)) 
gwr_fit <- gwr(CRIME~ INC + HOVAL + X + Y + WX + WY, data=col.poly1, coords=cbind(col.poly1$X,col.poly1$Y), adapt = 1, hatmatrix = TRUE)
gwr_fit
gwr1 <- gwr.sel(CRIME~ INC + HOVAL + X + Y + WX + WY, data=col.poly1, coords=cbind(col.poly1$X,col.poly1$Y), gweight=gwr.bisquare)
# Creando un objecto con el valor de Quasi-global R2
globalR2 <- (1 - (gwr_fit$results$rss/gwr_fit$gTSS))

results<-as.data.frame(gwr_fit$SDF)
head(results)

head(gwr_fit$SDF)
# Create a colour palette 
lm.palette <- colorRampPalette(c("red","orange", "white"), space = "rgb")
colours = c("dark blue", "blue", "red", "dark red")
col.poly$Intercept_coef <- gwr_fit$SDF@data[,2]
col.poly$HOVAL_coef <- gwr_fit$SDF@data$HOVAL
col.poly$INC_coef <- gwr_fit$SDF@data$INC
col.poly$localR2 <- gwr_fit$SDF@data$localR2
col.poly$pred <- gwr_fit$SDF@data$pred
spplot(col.poly, "HOVAL_coef", col.regions=lm.palette(20), cex=0.6, main = "HOVAL Coef") 
spplot(col.poly, "localR2", key.space = "right", col.regions = lm.palette(20), cuts = 7, 
       main = "CRIME_localR2")
spplot(col.poly, "INC_coef", key.space = "right", col.regions = lm.palette(20), cuts = 7, 
       main = "INC_coef")
spplot(col.poly, "Intercept_coef", key.space = "right", col.regions = lm.palette(20), cuts = 7, 
       main = "Intercept_coef")
names(gwr_fit$SDF)

#Not unsurprisingly, the relationship is significant at a very high confidence level.
# But, what about the residuals ? Is there a geographical patterning to where the model over- or under-predicts? Let's look! 

# Makes the data into a "mappable" format 
map = col.poly

# Lists the "mappable" data 
names(map) 

#calculate t-value
t = gwr_fit$SDF$HOVAL / gwr_fit$SDF$HOVAL_se  
map@data$t = t 
colours2=c("green","red","green") 

#estimated GWR t- values, red indicates a relationship that is not significant
spplot(map, "t", col.regions=lm.palette(20), cex=c(0.6,1,0.6), main = "t - Value") 

## predicted and standard error
map@data$pred <- gwr_fit$SDF@data$pred
map@data$pred.se <- gwr_fit$SDF@data$pred.se
map@data$localR2 <- gwr_fit$SDF@data$localR2
spplot(map, "pred",col.regions=lm.palette(20), main = "Predicted Value")
spplot(map, "pred.se", col.regions=lm.palette(20), main = "Standard Error")

## bubble plot 
bubble(gwr_fit$SDF, "localR2",fill = F, main = "local R2 Coef")


#######Global tests of geographical weighted regressions ###########
##Four related test statistics for comparing OLS and GWR models based on papers by Brunsdon, 
#Fotheringham and Charlton (1999) and Leung et al (2000), and a development from the GWR book (2002).


##Brunsdon, Fotheringham & Charlton (1999) ANOVA
BFC99.gwr.test(gwr_fit)
#Brunsdon, Fotheringham & Charlton (2002, pp. 91-2)
BFC02.gwr.test(gwr_fit)

# anova
anova(gwr_fit)

col.poly2 <- col.poly1
coordinates(col.poly2) <- c("X", "Y")
xx <- gwr(CRIME~ INC + HOVAL + X + Y + WX + WY, col.poly2, bandwidth = 9.454624, hatmatrix=TRUE)
xx

# local R2
map@data$localR2b <- xx$SDF@data$localR2
spplot(map, "localR2b", key.space="right", col.regions = lm.palette(20), cuts=7, main = "Local R2")

