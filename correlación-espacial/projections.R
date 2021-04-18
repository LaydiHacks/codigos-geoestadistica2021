
library(rgdal)
library(plotKML)

# Indique debajo del directorio que contiene los archivos de mapa "powiaty"
setwd("db/")

mapaA <- readOGR("db/", "powiaty")
#?C?mo se define la ubicaci?n de los v?rtices? Veamos p.e. los bordes del mapa (caja delimitadora) ...
mapaA@bbox
#?Qu? significan estos n?meros?
getCRS(mapaA)
proj4string(mapaA)
#CRS = (Coordinates Reference System), Sistema de referencia de coordenadas; implica 3 par?metros clave
  #1.projection
    projInfo(type = "proj")
    #larga lista de opciones ... generalmente "longlat", tambi?n "tmerc" (Transverse Mercator) como aqu?
  #2.datum
    projInfo(type = "datum")
    #a menudo vac?o (0,0); en general: punto de anclaje del sistema de coordenadas
  #3.ellipsoid
    projInfo(type = "ellps")
    #explicando la forma elipsoide de la Tierra (eje ecuatorial> eje polar); 
    #aproximaci?n global popular: WGS84
    
#Para cambiar entre CRS individuales:
mapaB <- spTransform(mapaA, "+proj=longlat")
#?M?s intuitiva ahora?
mapaB@bbox
#Ahora se ve c?modamente como viejos buenos grados de longitud y latitud
getCRS(mapaB)
proj4string(mapaB)
#Y el elipsoide WGS84 vino en lugar de GRS80

#Se puede imponer un conjunto predefinido de opciones de CRS, denominado c?digo EPSG:
#Se puede examinar en: http://www.epsg-registry.org
mapaC <- spTransform(mapaA, CRS("+init=epsg:4326"))
getCRS(mapaC)
#Tambi?n se puede imponer la proyecci?n desde un mapa diferente.
mapaD <- mapaA
mapaD <- spTransform(mapaD, getCRS(mapaB))


#Ahora compare lo que implican las dos proyecciones en el plano 2D:
plot(mapaA)
plot(mapaB)
#UTM (Universal Transverse Mercator) distorsiona la distancia, pero conserva ?ngulos y direcciones
#(bueno por ejemplo en la navegaci?n).

#Para ver el ejemplo del impacto de la proyecci?n en el trazado del mapa, ejecute el siguiente c?digo:

#install.packages("mapproj")
#install.packages("maps")
#install.packages("ggplot2")
library(maps); library(ggplot2); library(mapproj)
states <- map_data("state")
usamap <- ggplot(states, aes(x=long, y=lat, group=group)) +
  geom_polygon(fill="white", colour="black")
mapproject(0,0)
usamap + coord_map("mercator")
usamap + coord_map("azequidistant") 

#Ver m?s:
#https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/OverviewCoordinateReferenceSystems.pdf

