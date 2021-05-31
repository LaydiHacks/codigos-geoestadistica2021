# Transformacion de coordenadas
IDSTA <- read.table("stations_temp_xy_2008.csv", header=TRUE, sep=",", quote = "\"")
coordinates(IDSTA) <- ~LON+LAT
proj4string(IDSTA) <- CRS("+proj=longlat +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824")
# HRG2000 geoid [http://spatial-analyst.net/wiki/index.php?title=MGI_/_Balkans_coordinate_systems]
IDSTA.ll <- spTransform(IDSTA, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
writeOGR(IDSTA.ll, "gl_stations.kml", "IDSTA", "KML")
IDSTA.utm <- spTransform(IDSTA, CRS(utm33))
writeOGR(IDSTA.utm, "gl_stations.shp", "IDSTA", "ESRI Shapefile")
locs <- as.data.frame(IDSTA.utm)
names(locs)[c(5,6)] <- c("X", "Y")
str(locs)
# ST_ID     : Identification code of a station
# NAME      : Station name
# SP_CODE   : Station type "gl"-main, "kl"-climatological, "ks"-precipitation
# ELEV      : Station elevation in official database

# Lectura shapefile
library(maptools)
croatia.shp <- readShapePoly("D:/CARLOS/Doctorado/Tesis Doctoral/Kriging Universal Basado en Distancias/Articulos Universal Kriging/Articulo 1/CROACIA/croatia.shp")
croatia.shp <- readShapePoly("D:/........./croatia.shp")

# Anisotropia
library(intamap)
coordinates(Datos500A) <- ~x+y
AZ1_ILM <- estimateAnisotropy(Datos500A[Datos500A$zona==1,], "ILM")      # Para la especie LM zona1 hay anisotropia
AZ1_ILA <- estimateAnisotropy(Datos500A[Datos500A$zona==1,], "ILA")
AZ1_ISH <- estimateAnisotropy(Datos500A[Datos500A$zona==1,], "ISH")

library(geospt)
data(preci)
coordinates(preci) <- ~x+y
estimateAnisotropy(preci,prec)

#### Creación de la grilla con los puntos de predicción y ploteo
library(maptools)
pol_esp <- readShapePoly("D:/CARLOS/Doctorado/Tesis Isa/Borde.shp")
points(sp)
df.pts<-spsample(pol_esp, n=20000, type="regular")
plot(pol_esp)
points(df.pts)
