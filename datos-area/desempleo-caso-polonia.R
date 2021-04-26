#INSTALLING PACKAGES (previously installed packages commented - please uncomment if necessary)
# install.packages("rgdal")
# install.packages("sp")
# install.packages("spdep")

# clear the workspace, plots and console
rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014")

library(rgdal)
library(sp)
library(spdep)
library(spatialreg)


mapa <- readOGR("../db/powiaty.shp")
mapa <- spTransform(mapa, "+proj=longlat")
dane <- read.csv("../db/BDL_dane.csv", header = TRUE, sep = ";", dec = ",")
mapa@data$kod <- as.numeric(as.character(mapa@data$jpt_kod_je))
spatial_data <- merge(y = dane, x = mapa, by.y = "kod", by.x = "kod")
centroids <- coordinates(spatial_data)
rm(mapa)
rm(dane)

#Modelo lineal (A partir de lo específico ...?)
spatial_data@data$salario <- spatial_data@data$wynagrodzenie
spatial_data@data$l_salario <- log(spatial_data@data$salario)
spatial_data@data$desempleo <- spatial_data@data$bezrobocie
mod.lin <- lm(l_salario ~ desempleo, data = spatial_data)
summary(mod.lin)

cont1 <- poly2nb(spatial_data, queen = T)
W1_list <- nb2listw(cont1, style = "W")
cont2 <- dnearneigh(centroids, 0.0001, 60, row.names = spatial_data@data$kod, longlat = TRUE)
W2_list <- nb2listw(cont2, style = "W")
lm.morantest(mod.lin, W1_list, alternative = "greater")
lm.morantest(mod.lin, W2_list, alternative = "greater")
#NO!

# Modelo general GNS (A partir de lo general ...?)
mod.GNS1 <- sacsarlm(l_salario ~ desempleo, listw = W1_list, data = spatial_data, type = "sacmixed")
summary(mod.GNS1,Nagelkerke=T)
mod.GNS2 <- sacsarlm(l_salario ~ desempleo, listw = W2_list, data = spatial_data, type = "sacmixed")
summary(mod.GNS2,Nagelkerke=T)
#Depende de la matriz W! (¿Cuál es el mejor en términos socioeconómicos?)
#...y en el nivel de significancia

#¿Qué dicen Monte Carlo y la literatura?


#De lo específico a lo general (modelos de fuente única)

#Corremos todas las pruebas 
LM_W1 <- lm.LMtests(mod.lin, listw = W1_list, test = "all")
# Aparece para el Spatial Error, Spatial Lag, las Robustas para el Spatil Error y Spatil Lag 
#y luego aparace el SARMA
t(sapply(LM_W1, function(x) unlist(x[1:3])))
LM_W2 <- lm.LMtests(mod.lin, listw = W2_list, test = "all")
t(sapply(LM_W2, function(x) unlist(x[1:3])))
#Aqui ya aparecen las estructuras
#...¿qué modelo es el indicado mediante las pruebas? ¿Cuál de las pruebas robustas?

#Modelo LagSar
mod.SLM1 <- lagsarlm(l_salario ~ desempleo, listw = W1_list, data = spatial_data)
summary(mod.SLM1,Nagelkerke=T) #Este modelo es significativo, por lo cual podemos trabajar con este 

mod.SLM2 <- lagsarlm(l_salario ~ desempleo, listw = W2_list, data = spatial_data)
summary(mod.SLM2,Nagelkerke=T) #Este modelo es significativo, pero los residuos son autocorrelacionados, por lo cual este no lo usamos

#Spatial Error
mod.SEM1 <- errorsarlm(l_salario ~ desempleo, listw = W1_list, data = spatial_data)
summary(mod.SEM1) #En este caso, también es signitivativo
mod.SEM2 <- errorsarlm(l_salario ~ desempleo, listw = W2_list, data = spatial_data)
summary(mod.SEM2)

#Modelo LmSLX
mod.SLX1 <- lmSLX(l_salario ~ desempleo, listw = W1_list, data = spatial_data)
summar(mod.SLX1) 
mod.SLX2 <- lmSLX(l_salario ~ desempleo, listw = W2_list, data = spatial_data)
summary(mod.SLX2)

#¿Qué dice la teoría? ¿Qué hay de las pruebas?



#Modelos de múltiples fuentes: ¿tiene sentido estimarlos individualmente?

#SARAR
mod.SARAR1 <- sacsarlm(l_salario ~ desempleo, listw = W1_list, data = spatial_data)
summary(mod.SARAR1,Nagelkerke=T)
mod.SARAR2 <- sacsarlm(l_salario ~ desempleo, listw = W2_list, data = spatial_data)
summary(mod.SARAR2,Nagelkerke=T)

#Durbin Espacial
mod.SDM1 <- lagsarlm(l_salario ~ desempleo, listw = W1_list, data = spatial_data, type = "Durbin")
summary(mod.SDM1,Nagelkerke=T) #Tambien se puede utilizar este modelo 
mod.SDM2 <- lagsarlm(l_salario ~ desempleo, listw = W2_list, data = spatial_data, type = "Durbin")
summary(mod.SDM2,Nagelkerke=T)

mod.SDEM1 <- errorsarlm(l_salario ~ desempleo, listw = W1_list, data = spatial_data, etype = "emixed")
summary(mod.SDEM1,Nagelkerke=T)
mod.SDEM2 <- errorsarlm(l_salario ~ desempleo, listw = W2_list, data = spatial_data, etype = "emixed")
summary(mod.SDEM2,Nagelkerke=T)



