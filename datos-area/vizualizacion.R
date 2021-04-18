# https://mapshaper.org/ convierte a menor exactitud editor para mapa de datos
# http://web.sgh.waw.pl/~atoroj/
# Instale los paquetes que se requieren
#install.packages("rgdal")
#install.packages("sp")

# Paquetes para trabajar con mapas y visualizar datos en mapas
library(rgdal)
library(sp)

# Establecer directorio de trabajo: el mismo donde desempaquetamos los archivos descargados
#setwd("D:/CARLOS/Econometria Espacial/Doctorado Econom?a UN/Econometria Espacial Polonia")


# Mapa de importaci?n 1 - nivel de pobreza, proyecci?n correcta (source: Centro Central de Documentaci?n Geod?sica y Cartogr?fica, 
# http://www.codgik.gov.pl/index.php/darmowe-dane/prg.html)
mapa1 <- readOGR("db/", "powiaty")
# Este mapa es preciso y est? bien descrito, pero las coordenadas est?n codificadas de una manera diferente a la que necesitamos.
# Deber?amos recalcularlos en grados de longitud y latitud.
mapa1 <- spTransform(mapa1, "+proj=longlat")
plot(mapa1)

# Mapa de importaci?n 2: muchos niveles y pa?ses a la vez
#http://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/administrative-units-statistical-units/nuts#nuts13
mapa2 <- readOGR("db/", "NUTS_RG_01M_2013")
mapa2 <- spTransform(mapa2, "+proj=longlat")
plot(mapa2)
# Seleccionamos Polonia, Francia, Espa?a, Alemania, Reino Unido o el pa?s que deseen 
mapa2@data$NUTS_ID_char <- as.character(mapa2@data$NUTS_ID)
mapa2@data$country <- substr(mapa2@data$NUTS_ID_char, 1, 2) 
mapaPL <- mapa2[mapa2@data$country == "PL", ]
plot(mapaPL)
mapaFR <- mapa2[mapa2@data$country == "FR", ]
plot(mapaFR)
mapaES <- mapa2[mapa2@data$country == "ES", ]
plot(mapaES)
mapaDE <- mapa2[mapa2@data$country == "DE", ]
plot(mapaDE)
mapaUK <- mapa2[mapa2@data$country == "UK", ]
plot(mapaUK)
mapaRO <- mapa2[mapa2@data$country == "RO", ]
plot(mapaRO)


# Fronteras m?s amplias de ?reas administrativas
mapa2_NUTS1 <- mapa2[mapa2@data$STAT_LEVL_ == 1, ]
mapa2_NUTS2 <- mapa2[mapa2@data$STAT_LEVL_ == 2, ]
mapa2_NUTS3 <- mapa2[mapa2@data$STAT_LEVL_ == 3, ]
plot(mapa2_NUTS3, lwd = 0.3, border = rgb(0.7, 0.7, 0.7))
plot(mapa2_NUTS2, lwd = 2, add = TRUE)

mapaPLN2 <- mapa2_NUTS2[mapa2_NUTS2@data$country == "PL", ]
plot(mapaPLN2)

mapaESN2 <- mapa2_NUTS2[mapa2_NUTS2@data$country == "ES", ]
plot(mapaESN2)

mapaUKN2 <- mapa2_NUTS2[mapa2_NUTS2@data$country == "UK", ]
plot(mapaUKN2)

mapaRON2 <- mapa2_NUTS2[mapa2_NUTS2@data$country == "RO", ]
plot(mapaRON2)

mapaRON2@data$Infant_Mortality_Rates.2011 <- c(9.3,8.9,10.1,10.1,10.3,11.3,5.7,8.7)
writeOGR(mapaRON2,"D:/CARLOS/Econometria Espacial/Doctorado Econom?a UN/Parciales/Romania",layer="Romania.N2",driver="ESRI Shapefile")

Romania <- readOGR("D:/CARLOS/Econometria Espacial/Doctorado Econom?a UN/Parciales/Romania/Romania.N2.shp")
par(mai=c(0,0,0,0))
plot(Romania)
xy <- coordinates(Romania)
xy[5,1] <- xy[5,1]-0.5
xy[5,2] <- xy[5,2]+0.25
text(xy[,1],xy[,2],Romania$NUTS_ID, cex=1.2,col="red")

mapaUKN1 <- mapa2_NUTS1[mapa2_NUTS1@data$country == "UK", ]
plot(mapaUKN1)
names(UK.N1)[2] <- "NUTS_ID"
mapaUKN1@data <- merge(mapaUKN1@data, UK.N1@data, by = intersect("NUTS_ID","NUTS_ID"), by.x="NUTS_ID", all.x=TRUE, sort=F)
writeOGR(mapaUKN1,"D:/CARLOS/Econometria Espacial/Doctorado Econom?a UN/Parciales/UK",layer="UK.N1",driver="ESRI Shapefile")


UK.N1 <- readOGR("D:/CARLOS/Econometria Espacial/Doctorado Econom?a UN/Parciales/UK/United_Kingdom.N1.shp")
plot(UK.N1)
UK.N1@data$Labor.Productivity <- c(86.2,88.6,84.7,89.2,89.1,96.8,139.7,108.3,89.8,81.5,96.9,82.9)
UK.N1@data$Business.Birth.Rate <- c(11.2,11.1,10.5,10.3,10.5,10.5,14.6,10.8,9.6,9.3,10.9,6.5)



mapaUKN3 <- mapa2_NUTS3[mapa2_NUTS3@data$country == "UK", ]
plot(mapaUKN3)

#mapa de importaci?n 3 - formato R (rds; fuente: gadm.org)
mapa3 <- readRDS("POL_adm2.rds")
plot(mapa3)

#En ?ltima instancia, utilizamos el mapa de poviats proporcionado por CODGiK:
mapa <- mapa1
rm(mapa1, mapa2, mapa3, mapa2_NUTS2, mapa2_NUTS3)

#Importar otros datos
dane <- read.csv("BDL_dane.csv", header = TRUE, sep = ";", dec = ",")
mapa@data$kod <- as.numeric(as.character(mapa@data$jpt_kod_je))

#Ponga las bases de datos espaciales y econ?micas juntas, elimine las bases de datos parciales
spatial_data <- merge(y = dane, x = mapa, by.y = "kod", by.x = "kod")
rm(mapa)
rm(dane)

#Se ilustra con la variable "bezrobocie" desempleo
green_area <- rgb(24, 121, 104, 80, names = NULL, maxColorValue = 255)
pal <- colorRampPalette(c("white", green_area), bias = 1)
spplot(spatial_data, zcol = "bezrobocie", colorkey = TRUE, col.regions = pal(100), cuts = 99,
       par.settings = list(axis.line = list(col =  'transparent')),
       main = "Desempleo")
