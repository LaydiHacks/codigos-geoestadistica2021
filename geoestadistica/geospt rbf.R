#save.image("D:/CARLOS/Estadistica Espacial/Drive Estadistica Espacial/Geoestad?stica/Ejercicios/geospt")
load("D:/CARLOS/Estadistica Espacial/Drive Estadistica Espacial/Geoestad?stica/Ejercicios/geospt")

library(geospt)
data(preci)
coordinates(preci) <- ~x+y
# prediction case: a grid of points
puntos<-expand.grid(x=seq(min(preci$x),max(preci$x),0.05),y=seq(min(preci$y),max(preci$y),0.05))
coordinates(puntos) <- ~x+y

#Optimización de parametros
op.tps <- graph.rbf(prec~1, preci, eta.opt=TRUE, rho.opt=TRUE, n.neigh=9, func="TPS", 
          eta.dmax=2, rho.dmax=2, x0=c(0.1,0.1), iter=500)

op.tps1 <- graph.rbf(prec~x+y, preci, eta.opt=TRUE, rho.opt=TRUE, n.neigh=9, func="TPS", 
                    eta.dmax=2, rho.dmax=2, x0=c(0.1,0.1), iter=500)

#Establezco los parametros
pred.rbf <- rbf(prec~1, preci, eta=0.1503637, rho=0.01118794, newdata=puntos, n.neigh=10, func="TPS")
coordinates(pred.rbf) = c("x", "y")
gridded(pred.rbf) <- TRUE
# show prediction map
spplot(pred.rbf["var1.pred"], cuts=40, col.regions=bpy.colors(100),main = "", key.space=list(space="right", cex=0.8))

tabla <- rbf.tcv(prec~1, preci, eta=0.1503637, rho=0.01118794, n.neigh=9, func="TPS")
criteria.cv(tabla)

#Ejercicio con una base de datos real
library(geosptdb)
data(croatia)
data(croatiadb)

croatia.jan <- croatiadb[croatiadb$t==1,c(1:2,4)]
coordinates(croatia.jan) <- ~x+y
rbf.cv(MTEMP~1, croatia.jan, eta=1e-05, rho=0, n.neigh=10, func="M")
#1.23022
op.m <- graph.rbf(MTEMP~1, croatia.jan, eta.opt=T, rho.opt=T, n.neigh=10, func="M",
             eta.dmax=2, rho.dmax=2, iter=80)
# prediction case a grid of points
pts <- spsample(croatia, n=70000, type="regular")
pred.rbf <- rbf(MTEMP~1, croatia.jan, eta=1e-05, rho=0, newdata= pts,
                   n.neigh=10, func="M")
coordinates(pred.rbf) = c("x", "y")
gridded(pred.rbf) <- TRUE
spplot(pred.rbf["var1.pred"], cuts=40, scales = list(draw =T), col.regions
          =bpy.colors(100), key.space=list(space="right", cex=0.8))
