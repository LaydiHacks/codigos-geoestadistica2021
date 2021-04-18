library(lmtest)
library(RColorBrewer)
library(classInt)
library(spdep)
library(sphet)
library(RcmdrPlugin.epack)
library(rCarto)
library(pgirmess)
library(rgdal)


source("/cloud/project/funciones/northarrow.r")
source("/cloud/project/funciones/scalebar.R")
source("/cloud/project/funciones/moranbi.test.R")
source("/cloud/project/funciones/moran.cluster.R")
source("/cloud/project/funciones/moran.bi.R")
source("/cloud/project/funciones/moran.cluster.R")
source("/cloud/project/funciones/getis.cluster.R")
source("/cloud/project/funciones/localmoran.bi.R")
source("/cloud/project/funciones/moranbi.plot.R") 
source("/cloud/project/funciones/quantile.e.R")
source("/cloud/project/funciones/sp.na.omit.R")
source("/cloud/project/funciones/correlogram.d.R")


########################################
### Contiguidad Efecto Reina-Torre  ####
########################################

col.map <- readOGR("db/","Departamentos_IND")

## Ver los datos
#col.map@data

## Ver las Coordenadas
coords <- coordinates(col.map)


# Orden 1 
col_nbr <- poly2nb(col.map, queen=F)            # Torre
col_nbq1 <- poly2nb(col.map)                    # Reina Orden 1

col.lags <- nblag(col_nbq1, 2)                  # Orden 2
col.lags3 <- nblag(col_nbq1, 3)                 # Orden 3
col.lags10 <- nblag(col_nbq1, 10)


par(mar = c(0, 0, 0, 0), pty = "s")
plot(col.map)
# Efecto reina orden 1
plot(col_nbq1, coords, add=T, col=4, lwd=2)

plot(col.map)
# Efecto reina orden 2
plot(col.lags[[2]], coords, add=T, col=4, lwd=2)

plot(col.map)
# Efecto reina orden 3
plot(col.lags3[[3]], coords, add=T, col=4, lwd=2)

# Correlograma Moran a partir de matriz contiguidad espacial                                                  
sp.cr <- sp.correlogram(col_nbq1, col.map1@data$dbf.TAA04, order=6, method="corr", style="W", zero.policy=T) 
# Reina                                                                                                       
cor.s <- sp.correlogram(col_nbq1, col.map1@data$dbf.TAA04, order=6, method="I", style="W", zero.policy=T)    
plot(cor.s)                                                                                                   
# Torre                                                                                                       
cor.t <- sp.correlogram(col_nbr, col.map1@data$dbf.TAA04, order=6, method="I", style="W", zero.policy=T)     
plot(cor.t)                                                                                                   

################################
### K vecinos más cercanos  ####
################################

coords1 <- coordinates(col.map1)    
IDs <- row.names(as(col.map, "data.frame"))
#test <- knearneigh(coords, k=1, longlat = NULL, RANN=TRUE)
#knn2nb(test, row.names = NULL, sym = FALSE)

col_kn1<-knn2nb(knearneigh(coords, k=1), row.names=IDs)
col_kn2<-knn2nb(knearneigh(coords, k=2), row.names=IDs)
col_kn4<-knn2nb(knearneigh(coords, k=4), row.names=IDs)

par(mar = c(0, 0, 0, 0), pty = "s")
plot(col.map)
plot(col_kn1, coords, add=T, col=4, lwd=2)

plot(col.map)
plot(col_kn2, coords, add=T, col=4, lwd=2)

par(mar = c(0, 0, 0, 0), pty = "s")
plot(col.map)
plot(col_kn4, coords, add=T, col=4, lwd=2)

# Correlograma Moran k vecinos                                                                               
IDs <- row.names(col.map1@data)                                                                
col.kn5<-knn2nb(knearneigh(coords, k=5), row.names=IDs)                                                       
sp.cr <- sp.correlogram(col.kn5, col.map1$dbf.TAA04, order=5, method="corr", style="W", zero.policy=T) 
cor <- sp.correlogram(col.kn5, col.map1$dbf.TAA04, order=5, method="I", style="W", zero.policy=T)      
plot(cor)                                                                                                    

################################
##### Umbral (Distancias)  #####
################################

Dist <- unlist(nbdists(col.kn5, coords))
summary(Dist)
max_k1 <- max(Dist)

col_kd1<-dnearneigh(coords, d1=0, d2=100000, row.names=IDs)                   
col_kd2<-dnearneigh(coords, d1=0, d2=200000, row.names=IDs)                    
col_kd3<-dnearneigh(coords, d1=0, d2=300000, row.names=IDs)                    

par(mar = c(0, 0, 0, 0), pty = "s")
plot(col.map)
plot(col_kd1, coords, add=T, col=4, lwd=2)

plot(col.map)
plot(col_kd2, coords, add=T, col=4, lwd=2)

plot(col.map)
plot(col_kd3, coords, add=T, col=4, lwd=2)

# Correlograma Moran a partir de matriz distancias                                         
corD <- correlogram.d(coords,col.map1$dbf.TAA04,method="Moran",nbclass=10)           
corD$res                                                                                   
plot(corD$res,main="")  

# Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd
# BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral
par(mar=c(0,0,0,0))                            "YlOrBr"
quantile.e(x=col.map$A4, ic=5, digits=1, color="Spectral",style="quantile", border=col.map, north.arrow=c(70000,600000,1720000,0.7),scale.bar=list(470000,100000,500000,"m",0.7),leg.cex=0.8)

######################################################################################################

col.map1 <- readShapePoly("D:/CARLOS/Estadistica Espacial/Oscar/CEMM/Conflicto/deptos_colombia_magna_sirgas.shp", proj4string=CRS("+proj=longlat +datum=WGS84"))
plot(col.map1)
names(col.map1)
#col.map1@data[order(col.map1@data$NMG),]

AA <- read.dbf("D:/CARLOS/Estadistica Espacial/Oscar/CEMM/acciones armadas1.dbf")
#AA[order(AA$DEPARTAMEN),]

AVV <- read.dbf("D:/CARLOS/Estadistica Espacial/Oscar/CEMM/atencion a victimas.dbf")
# Área oficial KM2: AO
# Atención a Victimas de la Violencia 40 SMMLV-(año):  AVV

CE <- read.dbf("D:/CARLOS/Estadistica Espacial/Oscar/CEMM/cobertura en educacion.dbf")
# Área oficial KM2
# Cobertura Bruta en Educación por Nivel - Total: CBENT

CA <- read.dbf("D:/CARLOS/Estadistica Espacial/Oscar/CEMM/confrontaciones armadas.dbf")
# Área oficial KM2
# Confrontaciones Armadas por Año: CAA 

DCS <- read.dbf("D:/CARLOS/Estadistica Espacial/Oscar/CEMM/deficit cobertura salud.dbf")
# Área oficial KM2: AODCS
# Déficit de Cobertura en Salud: DCS2005

DFHR <- read.dbf("D:/CARLOS/Estadistica Espacial/Oscar/CEMM/despla forza hogares reci.dbf")
# Área oficial KM2: AODFHR
# Desplazamiento Forzado - Hogares Recibidos: DFHR

DFHE <- read.dbf("D:/CARLOS/Estadistica Espacial/Oscar/CEMM/despla forza hogares expu.dbf")
# Área oficial KM2: AODFHE
# Desplazamiento Forzado - Hogares Expulsados: DFHE

IDH <- read.dbf("D:/CARLOS/Estadistica Espacial/Oscar/CEMM/indice de desarrollo.dbf")
# Área oficial KM2: AOIDH
# Índice de Desarrollo Humano - IDH

IEGM <- read.dbf("D:/CARLOS/Estadistica Espacial/Oscar/CEMM/indice de gestion muni.dbf")
# Área oficial KM2: AOIE
# Índice de Eficacia de la Gestión Municipal: IEGM

PTU <- read.dbf("D:/CARLOS/Estadistica Espacial/Oscar/CEMM/pob total urb.dbf")
# Área oficial KM2

dfs <- list(col.map1@data[,-4], AA, matrix(AVV)[,-2], matrix(CE)[,-2], matrix(CA)[,-2], matrix(DCS)[,-2], matrix(DFHE)[,-2], matrix(DFHR)[,-2], matrix(PTU)[,-2], matrix(IDH)[,-2], matrix(IEGM)[,-2])

col.map1@data <- merge(col.map1@data, dfs, by = intersect("COD_DANE", "COD_DANE"), by.x="COD_DANE", all.x=TRUE, sort=F)

#Este comando permite guardar la información generar y el otro permite guardarla
#save.image("D:/CARLOS/Estadistica Espacial/Drive Estadistica Espacial/Datos Área/Ejercicios/Col")
#load("D:/CARLOS/Estadistica Espacial/Drive Estadistica Espacial/Datos Área/Ejercicios/Col")

