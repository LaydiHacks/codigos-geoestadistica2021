#################################
######    Declustering     #######
#################################

library(geospt)
data(preci)
as.matrix(dist(preci[,2:3]))

library(FNN)
m <- get.knn(data.frame(x=preci$x, y=preci$y), 3)       # $nn.index
str(m) 
knns.ndx <- data.frame(m["nn.index"])
knns.ndi <- data.frame(m["nn.dist"])
head(knns.ndx)
tail(knns.ndx) 
Dist.knn <- apply(knns.ndi,1,sum)

preci$S.dist <- Dist.knn
preci$Wi <- preci$S.dist/sum(preci$S.dist)
z.prom <- preci$Wi%*%preci$prec
weighted.mean(preci$prec,preci$Wi,na.rm=FALSE)
var.c <- (preci$prec-z.prom)%*%(preci$prec-z.prom)/10
var.m <- (preci$Wi%*%(preci$prec-z.prom)^2)/sum(preci$Wi)


######################################
###### Anamorfosis Gaussiana   #######
######################################

library(RGeostats)
library(gstat)
library(sp)
data(meuse)

#01. Fit Gaussian Anamorphosis 
mdb = meuse[,c("x","y","zinc")] #must be 3 columns: x,y,z
mdb.rgdb = db.create(mdb,ndim=2,autoname=F)
mdb.herm = anam.fit(mdb.rgdb,name="zinc",type="gaus")
mdb.hermtrans = anam.z2y(mdb.rgdb,names="zinc",anam=mdb.herm)
zinc.trans = mdb.hermtrans@items$Gaussian.zinc 
#Sin transformar
shapiro.test(mdb$zinc)
#Transformación Box-Cox
library(RcmdrPlugin.epack)
bc2(mdb$zinc)
shapiro.test(mdb$zinc^-0.3)
#Con transformación: Anamorfosis Gaussiana
shapiro.test(zinc.trans)

#02. Variogram fit to gaussian transformed data using gstat 
meuse$zn.trans= zinc.trans
coordinates(meuse) = ~x+y
zn.svgm <- variogram(zn.trans~1,data=meuse)
vgm1 = vgm(1.17,"Sph",1022, 0.0925)
plot(zn.svgm,vgm1)
meuse.gstat <- gstat(id="zn.trans", formula=zn.trans~1, model=vgm1, data=meuse,nmin=2, nmax=30)

#03. Kriging predictions with grid (gstat)
data(meuse.grid)
coordinates(meuse.grid) = ~x+y
gridded(meuse.grid)=T
zntrans.ok = predict(meuse.gstat,meuse.grid) 

#04. Back transform using RGeostats
zn.bts = cbind(coordinates(zntrans.ok),zntrans.ok@data)
zn.bts.db = db.create(zn.bts,autoname = F)
tempdb = anam.y2z(zn.bts.db,names="zn.trans.pred",anam = mdb.herm)
#Prediction map
meuse.grid@data$pred.zinc <- tempdb@items$Raw.zn.trans.pred
spplot(meuse.grid,"pred.zinc")