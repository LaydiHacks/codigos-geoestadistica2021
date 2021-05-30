# Exploracion de datos geoestadisticos
setwd("D:/CARLOS/Estadistica Espacial/Estadistica Espacial/Geoestadística/Ejercicios")
source("aquifer.r")
attach(aquifer)
plot(easting,northing,pch=20,xlab="longitud",ylab="latitud")
plot(easting,northing,type="n",xlab="longitud",ylab="latitud")
text(easting,northing,head,cex=0.6)
library(scatterplot3d)
scatterplot3d(easting,northing,head,xlab="longitud",ylab="latitud",zlab="altura")
detach(aquifer)

# Estimacion del variograma
library(geoR)
aq.geo<-as.geodata(aquifer)
plot.geodata(aq.geo)

# estimador clásico
aq.cl1<-variog(aq.geo,option="cloud")
plot(aq.cl1,pch=20)
aq.v1<-variog(aq.geo, uvec=seq(0,260,20))                # uvec: inin dist, fin dist, tamaño rezago (lag)
plot(aq.v1,pch=20)
aq.vc1<-variog(aq.geo, uvec=seq(0,260,20), bin.cloud=T)  # las nubes de semivarianza son incluidas en cada bin, util para analizar la variabilidad local e inclusive detectar atipicos.
plot(aq.vc1,bin.cloud=T)

# estimador robusto
aq.cl2<-variog(aq.geo,option="cloud",estimator.type="modulus")
plot(aq.cl2,pch=20)
aq.v2<-variog(aq.geo,uvec=seq(0,260,20),estimator.type="modulus")
plot(aq.v2,pch=20)
aq.vc2<-variog(aq.geo,uvec=seq(0,260,20),estimator.type="modulus",bin.cloud=T)
plot(aq.vc2,bin.cloud=T)


# Pruebas con varios intervalos
plot(variog(aq.geo,uvec=seq(0,260,10)))
plot(variog(aq.geo,uvec=seq(0,260,5)))
plot(variog(aq.geo,uvec=seq(0,260,2)))
plot(variog(aq.geo,uvec=seq(0,260,1)))
plot(variog(aq.geo,uvec=seq(0,280,40)))
plot(variog(aq.geo,uvec=seq(0,300,50)))
plot(variog(aq.geo,uvec=seq(0,320,80)))

# Superposición del variograma empírico y de un modelo de variograma
aq.exp.ml<-likfit(geodata = aq.geo, ini = c(500000, 150))
aq.sph.ml<-likfit(geodata = aq.geo, ini = c(500000, 150),cov.model="sph")
aq.mat.ml<-likfit(geodata = aq.geo, ini = c(500000, 150),cov.model="mat",kappa=1.5)
aq.mat2.ml<-likfit(geodata = aq.geo, ini = c(500000, 150),cov.model="mat",kappa=1,fix.nugget=T)
aq.cir.ml<-likfit(geodata = aq.geo, ini = c(500000, 150),cov.model="cir")
aq.gau.ml<-likfit(geodata = aq.geo, ini = c(500000, 50),cov.model="gau")
aq.cub.ml<-likfit(geodata = aq.geo, ini = c(500000, 150),cov.model="cub")
aq.pow.ml<-likfit(geodata = aq.geo, ini = c(500000, 150),cov.model="powered.exponential",kappa=1.75)
aq.pow2.ml<-likfit(geodata = aq.geo, ini = c(500000, 150),cov.model="powered.exponential",kappa=1.75,fix.nugget=T)
plot(aq.v1)
lines(aq.pow2.ml,max.dist=250,lwd=3,col='red')
lines(aq.mat2.ml,max.dist=250,lwd=3,col='blue')
lines(aq.pow.ml,max.dist=250,lwd=3,col='green')
lines(aq.mat.ml,max.dist=250,lwd=3,col='yellow')
lines(aq.cub.ml,max.dist=250,lwd=3,col='orange')
lines(aq.gau.ml,max.dist=250,lwd=3,col='cyan')
lines(aq.cir.ml,max.dist=250,lwd=3,col='grey')
lines(aq.exp.ml,max.dist=250,lwd=3,col='magenta')
lines(aq.sph.ml,max.dist=250,lwd=3,col='pink')
# los mejores son el 'mat' y el 'powered.exponential'

# Diferentes metodos de ajuste del variograma
aq.mat.ml<-likfit(geodata = aq.geo, ini = c(500000, 150),cov.model="mat",kappa=1.5)
aq.mat.rml<-likfit(geodata = aq.geo, ini = c(500000, 150),cov.model="mat",kappa=1.5,method='RML')
aq.mat.ols<-variofit(vario = aq.v1, ini = c(300000, 100),cov.model="mat",kappa=1.5,weights="equal",minimisation.function="optim")
aq.mat.wls<-variofit(vario = aq.v1, ini = c(500000, 150),cov.model="mat",kappa=1.5,weights="npairs")
aq.pow.ml<-likfit(geodata = aq.geo, ini = c(500000, 150),cov.model="powered.exponential",kappa=1.75)
aq.pow.rml<-likfit(geodata = aq.geo, ini = c(500000, 150),cov.model="powered.exponential",kappa=1.75,method='RML')
aq.pow.ols<-variofit(vario = aq.v1, ini = c(500000, 150),cov.model="powered.exponential",kappa=1.75,weights="equal",minimisation.function="optim")
aq.pow.wls<-variofit(vario = aq.v1, ini = c(500000, 150),cov.model="powered.exponential",kappa=1.75,weights="npairs")
plot(aq.v1)
lines(aq.mat.ml,max.dist=250,lwd=3)
lines(aq.mat.rml,max.dist=250,lwd=3,lty=2)
lines(aq.mat.ols,max.dist=250,lwd=3,lty=3)
lines(aq.mat.wls,max.dist=250,lwd=3,lty=4)
legend(locator(1),legend=c('ML','RML','OLS','WLS'),,lwd=c(3,3,3,3),lty=c(1,2,3,4))
plot(aq.v1)
lines(aq.pow.ml,max.dist=250,lwd=3)
lines(aq.pow.rml,max.dist=250,lwd=3,lty=2)
lines(aq.pow.ols,max.dist=250,lwd=3,lty=3)
lines(aq.pow.wls,max.dist=250,lwd=3,lty=4)
legend(20,2000000,legend=c('ML','RML','OLS','WLS'),lwd=c(3,3,3,3),lty=c(1,2,3,4))

var4 <- variog4(aq.geo, max.dist=250)
plot(var4)

library(intamap)
coordinates(aquifer) <- ~easting+northing
estimateAnisotropy(aquifer,"head")

# Calculo de superficies de tendencia
library(spatial)
attach(aquifer)
aq.ls<-surf.ls(3,easting,northing,head)
aq.trsurf<-trmat(aq.ls, -150, 120, 0, 190, 100)

summary(lm(aquifer$head~aquifer$easting+aquifer$northing))

# Representación de superficies de tendencia
par(pty="s",mar=c(2,2,2,2))
contour(aq.trsurf)
points(easting,northing,pch=20)
par(mar=c(0,0,0,0))
image(aq.trsurf)
points(easting,northing,pch=20)
par(mfrow=c(1,1))
persp(aq.trsurf)
persp(aq.trsurf,theta=60,phi=30,col=2,ltheta=-20,shade=0.25,xlab="longitud",ylab="latitud")
detach(aquifer)

# Eliminacion de tendencias
aq.sin<-aquifer[,3]-predict(aq.ls,aquifer[,1],aquifer[,2])
aqs.geo<-aq.geo
aqs.geo$data<-aq.sin
plot.geodata(aqs.geo)
plot(variog(aqs.geo),pch=20)

aqs.v1<-variog(aqs.geo,uvec=seq(0,100,10),max.dist=120)
plot(aqs.v1)

# Estimacion del variograma sin tendencia
aqs.exp.ml<-likfit(geodata = aqs.geo, ini = c(30000, 50))
aqs.sph.ml<-likfit(geodata = aqs.geo, ini = c(30000, 50),cov.model="sph")
aqs.mat.ml<-likfit(geodata = aqs.geo, ini = c(30000, 50),cov.model="mat",kappa=1.5)
aqs.mat2.ml<-likfit(geodata = aqs.geo, ini = c(30000, 50),cov.model="mat",kappa=1,fix.nugget=T)
aqs.cir.ml<-likfit(geodata = aqs.geo, ini = c(30000, 50),cov.model="cir")
aqs.gau.ml<-likfit(geodata = aqs.geo, ini = c(30000, 50),cov.model="gau")
aqs.cub.ml<-likfit(geodata = aqs.geo, ini = c(30000, 50),cov.model="cub")
aqs.pow.ml<-likfit(geodata = aqs.geo, ini = c(30000, 50),cov.model="powered.exponential",kappa=1.75)
aqs.pow2.ml<-likfit(geodata = aqs.geo, ini = c(30000, 50),cov.model="powered.exponential",kappa=1.75,fix.nugget=T)
plot(aqs.v1)
lines(aqs.pow2.ml,max.dist=100,lwd=2,col='red')
lines(aqs.mat2.ml,max.dist=100,lwd=2,col='blue')
lines(aqs.pow.ml,max.dist=100,lwd=2,col='green')
lines(aqs.mat.ml,max.dist=100,lwd=2,col='yellow')
lines(aqs.cub.ml,max.dist=100,lwd=2,col='orange')
lines(aqs.gau.ml,max.dist=100,lwd=2,col='cyan')
lines(aqs.cir.ml,max.dist=100,lwd=2,col='grey')
lines(aqs.exp.ml,max.dist=100,lwd=2,col='magenta')
lines(aqs.sph.ml,max.dist=100,lwd=2,col='pink')

aqs.exp.wls<-variofit(vario = aqs.v1, ini = c(30000, 50),weights="npairs")
aqs.sph.wls<-variofit(vario = aqs.v1, ini = c(30000, 50),weights="npairs",cov.model="sph")
aqs.mat.wls<-variofit(vario = aqs.v1, ini = c(30000, 50),weights="npairs",cov.model="mat",kappa=1.5)
aqs.mat2.wls<-variofit(vario = aqs.v1, ini = c(30000, 50),weights="npairs",cov.model="mat",kappa=1,fix.nugget=T)
aqs.cir.wls<-variofit(vario = aqs.v1, ini = c(30000, 50),weights="npairs",cov.model="cir")
aqs.gau.wls<-variofit(vario = aqs.v1, ini = c(30000, 50),weights="npairs",cov.model="gau")
aqs.cub.wls<-variofit(vario = aqs.v1, ini = c(30000, 50),weights="npairs",cov.model="cub")
aqs.pow.wls<-variofit(vario = aqs.v1, ini = c(30000, 50),weights="npairs",cov.model="powered.exponential",kappa=1.75)
aqs.pow2.wls<-variofit(vario = aqs.v1, ini = c(30000, 50),weights="npairs",cov.model="powered.exponential",kappa=1.75,fix.nugget=T)
plot(aqs.v1)
lines(aqs.pow2.wls,max.dist=100,lwd=2,col='red')
lines(aqs.mat2.wls,max.dist=100,lwd=2,col='blue')
lines(aqs.pow.wls,max.dist=100,lwd=2,col='green')
lines(aqs.mat.wls,max.dist=100,lwd=2,col='yellow')
lines(aqs.cub.wls,max.dist=100,lwd=2,col='orange')
lines(aqs.gau.wls,max.dist=100,lwd=2,col='cyan')
lines(aqs.cir.wls,max.dist=100,lwd=2,col='grey')
lines(aqs.exp.wls,max.dist=100,lwd=2,col='magenta')
lines(aqs.sph.wls,max.dist=100,lwd=2,col='pink')

aqs.exp.ml<-likfit(geodata = aqs.geo, ini = c(30000, 50))
aqs.exp.rml<-likfit(geodata = aqs.geo, ini = c(30000, 50),method='RML')
aqs.exp.ols<-variofit(vario = aqs.v1, ini = c(30000, 50),weights="equal",minimisation.function="optim")
aqs.exp.wls<-variofit(vario = aqs.v1, ini = c(30000, 50),weights="npairs")
aqs.sph.ml<-likfit(geodata = aqs.geo, ini = c(30000, 50),cov.model="sph")
aqs.sph.rml<-likfit(geodata = aqs.geo, ini = c(30000, 50),method='RML',cov.model="sph")
aqs.sph.ols<-variofit(vario = aqs.v1, ini = c(30000, 50),weights="equal",minimisation.function="optim",cov.model="sph")
aqs.sph.wls<-variofit(vario = aqs.v1, ini = c(30000, 50),weights="npairs",cov.model="sph")

plot(aqs.v1)
lines(aqs.exp.ml,max.dist=100,lwd=2)
lines(aqs.exp.rml,max.dist=100,lwd=2,lty=2)
lines(aqs.exp.ols,max.dist=100,lwd=2,lty=3)
lines(aqs.exp.wls,max.dist=100,lwd=2,lty=4)
legend(60,10000,legend=c('ML','RML','OLS','WLS'),lty=c(1,2,3,4))

plot(aqs.v1)
lines(aqs.sph.ml,max.dist=100,lwd=2)
lines(aqs.sph.rml,max.dist=100,lwd=2,lty=2)
lines(aqs.sph.ols,max.dist=100,lwd=2,lty=3)
lines(aqs.sph.wls,max.dist=100,lwd=2,lty=4)
legend(60,10000,legend=c('ML','RML','OLS','WLS'),lty=c(1,2,3,4))


loci <- expand.grid(seq(-140,110,l=31), seq(10,180,l=31))
plot(aq.geo$coords)
points(loci,cex=0.3)
kc1<-krige.conv(aq.geo,locations=loci,krige=krige.control(
 cov.pars=aq.mat.ml$cov.pars,nugget=aq.mat.ml$nugget))
par(mfrow=c(2,2))
image.kriging(kc1,loc=loci,main='estimacion kriging')
image.kriging(kc1,loc=loci,val=sqrt(kc1$krige.var),
 main='error estandar')
persp.kriging(kc1,loc=loci,main='estimacion kriging',phi=30,theta=45)
persp.kriging(kc1,loc=loci,val=sqrt(kc1$krige.var),
 main='error estandar')

evaltend<-function(superf,puntos){
 predict(superf,puntos[,1],puntos[,2])}
kc2<-krige.conv(aqs.geo,locations=loci,krige=krige.control(
 cov.pars=aqs.exp.ols$cov.pars,nugget=aqs.exp.ols$nugget))
par(mfrow=c(2,2))
image.kriging(kc2,loc=loci,val=kc2$predict+evaltend(aq.ls,loci),
 main='estimacion kriging')
image.kriging(kc2,loc=loci,val=sqrt(kc2$krige.var),
 main='error estandar')
persp.kriging(kc2,loc=loci,val=kc2$predict+evaltend(aq.ls,loci),
 main='estimacion kriging',phi=30,theta=45)
persp.kriging(kc2,loc=loci,val=sqrt(kc2$krige.var),
 main='error estandar')

kc3<-krige.conv(aq.geo,locations=loci,krige=krige.control(trend.d=
 ~evaltend(aq.ls,aq.geo$coords),trend.l=~evaltend(aq.ls,loci),
 cov.pars=aqs.exp.ols$cov.pars,nugget=aqs.exp.ols$nugget))
par(mfrow=c(2,2))
image.kriging(kc3,loc=loci,main='estimacion kriging')
image.kriging(kc3,loc=loci,val=sqrt(kc3$krige.var),
 main='error estandar')
persp.kriging(kc3,loc=loci,main='estimacion kriging',phi=30,theta=45)
persp.kriging(kc3,loc=loci,val=sqrt(kc3$krige.var),
 main='error estandar')

source('carbon.r')
ca.geo<-as.geodata(carbon)
plot.geodata(ca.geo)
plot(carbon$x,carbon$y,type='n',lab=c(16,23,71),xlab='x',ylab='y')
text(carbon$x,carbon$y,format(round(carbon$carbon,1)),cex=0.5)
scatterplot3d(carbon$x,carbon$y,carbon$carbon,type='h',pch=16,
 xlab='x',ylab='y',zlab='carbon')

source('mpulido.r')
matca<-matrix(nrow=16,ncol=23)
for (i in 1:208)
   matca[carbon[i,1],carbon[i,2]]<- carbon[i,3]
ca.mp<-mpulido(matca)

sup.mp<-matrix(nrow=16,ncol=23)
for (i in 1:16)
  for (j in 1:23)
      sup.mp[i,j]<- ca.mp$overall+ca.mp$row[i]+ca.mp$col[j]
persp(sup.mp,theta=45)
image(sup.mp)

casmp.geo<-ca.geo
for (i in 1:208)
    casmp.geo$data[i]<-ca.geo$data[i]-
       sup.mp[ca.geo$coords[i,1],ca.geo$coords[i,2]]
plot.geodata(casmp.geo)

casmp.v<-variog(casmp.geo,uvec=seq(0,26,1),max.dist=20)
plot(casmp.v)
casmp.nug.ml<-likfit(geodata=casmp.geo,ini=c(1,0),
 cov.model='pure.nugget')
lines(casmp.nug.ml,max.dist=20,lwd=2)

kca<-krige.conv(casmp.geo,locations=c(8.2,10.6),krige=krige.control(
 cov.pars=casmp.nug.ml$cov.pars,nugget=casmp.nug.ml$nugget))
pred.ca<-kca$predict+ca.mp$overall+ca.mp$row[8]+ca.mp$col[10]+
 0.2*(ca.mp$row[9]-ca.mp$row[8])+0.6*(ca.mp$col[11]-ca.mp$col[10])
c(pred.ca-1.96*sqrt(kca$krige.var),pred.ca+1.96*sqrt(kca$krige.var))



# Otro ejemplo: Altitud de un terreno

source("altitud.r")
attach(altitud)
plot(x,y,type="n")
text(x,y,alt)
library(scatterplot3d)
scatterplot3d(x,y,alt)
detach(altitud)

alt.geo<-as.geodata(altitud)
plot(alt.geo)