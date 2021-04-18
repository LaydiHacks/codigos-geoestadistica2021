install.packages("rgdal")
install.packages("plotKML")
library(rgdal)
library(plotKML)
# Indique debajo del directorio que contiene los archivos de mapa "powiaty"
setwd("C:\Users\leonh\Desktop\GEOESTADISTICA\Nueva carpeta")

mapaA <- readOGR(".", "powiaty")
#¿Cómo se define la ubicación de los vértices? Veamos p.e. los bordes del mapa (caja delimitadora) ...
mapaA@bbox
#¿Qué significan estos números?
getCRS(mapaA)
proj4string(mapaA)
#CRS = (Coordinates Reference System), Sistema de referencia de coordenadas; implica 3 parámetros clave
  #1.projection
    projInfo(type = "proj")
    #larga lista de opciones ... generalmente "longlat", también "tmerc" (Transverse Mercator) como aquí
  #2.datum
    projInfo(type = "datum")
    #a menudo vacío (0,0); en general: punto de anclaje del sistema de coordenadas
  #3.ellipsoid
    projInfo(type = "ellps")
    #explicando la forma elipsoide de la Tierra (eje ecuatorial> eje polar); 
    #aproximación global popular: WGS84
    
#Para cambiar entre CRS individuales:
mapaB <- spTransform(mapaA, "+proj=longlat")
#¿Más intuitiva ahora?
mapaB@bbox
#Ahora se ve cómodamente como viejos buenos grados de longitud y latitud
getCRS(mapaB)
proj4string(mapaB)
#Y el elipsoide WGS84 vino en lugar de GRS80

#Se puede imponer un conjunto predefinido de opciones de CRS, denominado código EPSG:
#Se puede examinar en: http://www.epsg-registry.org
mapaC <- spTransform(mapaA, CRS("+init=epsg:4326"))
getCRS(mapaC)
#También se puede imponer la proyección desde un mapa diferente.
mapaD <- mapaA
mapaD <- spTransform(mapaD, getCRS(mapaB))


#Ahora compare lo que implican las dos proyecciones en el plano 2D:
plot(mapaA)
plot(mapaB)
#UTM (Universal Transverse Mercator) distorsiona la distancia, pero conserva ángulos y direcciones
#(bueno por ejemplo en la navegación).

#Para ver el ejemplo del impacto de la proyección en el trazado del mapa, ejecute el siguiente código:

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

#Ver más:
#https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/OverviewCoordinateReferenceSystems.pdf

