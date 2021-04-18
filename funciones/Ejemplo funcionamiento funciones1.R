############################################################
### Funcionamiento funciones para dependencia espacial:  ###
############################################################

setwd("D:/CARLOS/Estadistica Espacial/Drive Estadistica Espacial/Datos Área/Funciones")
library(spatialreg)
library(spdep)
library(rgdal)
columbus <- readOGR(system.file("shapes/columbus.shp", package="spData")[1])
col_nbq <- poly2nb(columbus)
par.lags1 <- nblag(col_nbq, 2)                  # Orden 2
e.lw2 <- nb2listw(par.lags1[[2]], style="W",zero.policy=T)
a.lw <- nb2listw(col_nbq, style="W")
CRIME <- columbus$CRIME
INC <- columbus$INC

####################################
########     moran.bi      #########
####################################
source("moran.bi.R")
W <- as.matrix(as_dgRMatrix_listw(a.lw))
moran.bi(columbus$CRIME,columbus$INC,a.lw,zero.policy =T)

################################################
########  moranbi.test, localmoran.bi  #########
################################################

source("moranbi.test.R")
source("randomize_vector.R")
source("moran.bi.R")
source("localmoran.bi.R")
set.seed(123)
MBCrime <- moranbi.test(columbus$CRIME,columbus$INC,a.lw,999,graph=T,zero.policy =T,N=1000)
moranbi.test(columbus$INC,columbus$HOVAL,a.lw,999,graph=T,zero.policy =T,N=1000)
BiLisa.IC <- localmoran.bi(columbus$CRIME,columbus$INC,a.lw,zero.policy =T,mlvar=F,alternative = "greater")

########################################
########     moranbi.plot      #########
########################################

source("moranbi.plot.R")
# Editando las etiquetas de los ejes
CRIME <- as.vector(scale(columbus$CRIME))
INCOME <- as.vector(scale(columbus$INC))
moranbi.plot(CRIME,INCOME,quiet =F,zero.policy =F,listw=a.lw)
# Sin editar la etiqueta de los ejes
moranbi.plot(as.vector(scale(columbus$CRIME)),as.vector(scale(columbus$INC)),quiet =F,zero.policy =F,listw=a.lw)

##################################
########  moran.cluster  #########
##################################

source("moran.cluster.R")
# LISA Cluster Map: COLUMBUS
# Si desean el mapa de significancia deben escribir T, de lo contrario F
x11()
moran.cluster(CRIME, a.lw, zero.policy = FALSE, columbus, significant=T)
moran.cluster(CRIME, e.lw2, zero.policy = FALSE, columbus, significant=T)

###########################################
########     moranbi.cluster      #########
###########################################

source("moranbi.cluster.R")
# BiLISA Cluster Map: COLUMBUS
# Si desean el mapa de significancia deben escribir T, de lo contrario F
x11()
moranbi.cluster(CRIME, INC, a.lw, zero.policy = FALSE, columbus, significant=T)
moranbi.cluster(CRIME, INC, e.lw2, zero.policy = FALSE, columbus, significant=T)
moranbi.cluster(CRIME, INC, a.lw, zero.policy = FALSE, columbus, significant=T, alternative="two.sided")

#########################################
########     getis.cluster      #########
#########################################

source("getis.cluster.R")
# Getis Cluster Map: COLUMBUS
# Si desean el mapa de significancia deben escribir T, de lo contrario F
x11()
getis.cluster(CRIME, a.lw, zero.policy = FALSE, columbus, significant=T)
getis.cluster(CRIME, e.lw2, zero.policy = FALSE, columbus, significant=T)

#########################################
######## Correlograma Bivariado #########
#########################################

source("moranbi1.test.R")
source("geary.bi.R")
source("moran.bi1.R")
source("spcorrelogram.bi.R")
spc <- spcorrelogram.bi(col_nbq,columbus$CRIME,columbus$INC,a.lw,order=7,method="I",zero.policy = T,alternative="two.sided")
spc
x11()
plot(spc)

moranbi1.test(columbus$CRIME,columbus$INC,a.lw,zero.policy =T,alternative="less")
moranbi1.test(columbus$CRIME,columbus$INC,a.lw,zero.policy =T,alternative="two.sided")
moranbi1.test(columbus$CRIME,columbus$INC,a.lw,zero.policy =T,alternative="greater")

