library(geospt) #Libreria de Melo
data(preci)
coordinates(preci) <- ~x+y
# prediction case: a grid of points
puntos<-expand.grid(x=seq(min(preci$x),max(preci$x),0.05),
                        y=seq(min(preci$y),max(preci$y),0.05))
coordinates(puntos) <- ~x+y
pred.rbf <- rbf(prec~x+y, preci, eta=0.1460814, rho=0, newdata=puntos,
                  n.neigh=10, func="TPS")
coordinates(pred.rbf) = c("x", "y")
gridded(pred.rbf) <- TRUE
# show prediction map
spplot(pred.rbf["var1.pred"], cuts=40, col.regions=bpy.colors(100),
           main = "", key.space=list(space="right", cex=0.8))

library(geospt)
data(preci)
coordinates(preci)<-~x+y
rbf.cv(prec~1, preci, eta=0.1460814, rho=0, n.neigh=9, func="TPS")
rbf.tcv(prec~x+y, preci, eta=0.1460814, rho=0, n.neigh=9, func="TPS")
graph.rbf(prec~x+y, preci, eta.opt=T, rho.opt=T, n.neigh=9, func=
            "TPS", np=80, eta.dmax=2, rho.dmax=2, P.T=TRUE, iter=100)


library(geosptdb)
data(croatia)
data(croatiadb)
croatia.jan <- croatiadb[croatiadb$t==1,c(1:2,4)]
coordinates(croatia.ene) <- ~x+y
rbf.cv(MTEMP~1, croatia.jan, eta=1e-05, rho=0, n.neigh=10, func="M")
1.23022
gp <- graph.rbf(MTEMP~1, croatia.jan, eta.opt=T, rho.opt=T, n.neigh=10, func="M",
             eta.dmax=2, rho.dmax=2, iter=80)
# prediction case a grid of points
pts <- spsample(croatia, n=70000, type="regular")
pred.rbf <- rbf(MTEMP~1, croatia.jan, eta=1e-05, rho=0, newdata= pts,
                   n.neigh=10, func="M")
coordinates(pred.rbf) = c("x", "y")
gridded(pred.rbf) <- TRUE
spplot(pred.rbf["var1.pred"], cuts=40, scales = list(draw =T), col.regions
          =bpy.colors(100), key.space=list(space="right", cex=0.8))



# Loading Croatia data
data(croatiadb)
coordinates(croatiadb) <- ~x+y

# prediction case: one point
point <- data.frame(670863,5043464,5)
names(point) <- c("x","y","t")

coordinates(point) <- ~x+y
idwST(MTEMP~1, data=croatiadb, newdata=point, n.neigh=60, C=1, factor.p=2)

## Not run: 
# prediction case: a grid of points Croatia (year 2008)
data(croatia)
points <- spsample(croatia, n=5000, type="regular")

data(croatiadb)
coordinates(croatiadb)<-~x+y

GridsT <- vector(mode = "list", length = 12)

for(i in 1:12){ 
  GridsT[[i]] <- data.frame(points@coords,i)
  names(GridsT[[i]]) <- c("x","y","t")
}

idw.croatia <- data.frame(matrix(NA, ncol = 14, nrow=nrow(GridsT[[1]])))
pb <- txtProgressBar(min = 0, max = 12, char = "=", style = 3)
for(i in 1:12){ 
  coordinates(GridsT[[i]]) <- c("x", "y")
  idw.croatia[,i+2] <- idwST(MTEMP~1, croatiadb, newdata=GridsT[[i]], n.neigh=10, C=1, 
                             factor.p=2, progress=FALSE)[,4]                  
  setTxtProgressBar(pb, i)
}
close(pb)

idw.croatia[,1:2] <- GridsT[[1]]@coords
nam <- paste(c("ENE","FEB","MAR","ABR","MAY","JUN","JUL","AGO","SEP","OCT","NOV","DIC"),
             2008,sep="")
names(idw.croatia) <- c("x","y",nam)

coordinates(idw.croatia) <- c("x", "y")
gridded(idw.croatia) <- TRUE

# show prediction map
pal2 <- colorRampPalette(c("blue3", "wheat1", "red3"))

p1 <- spplot(idw.croatia[,1:12], cuts=30, col.regions=pal2(35), colorkey=F, 
             scales = list(draw =T,cex=0.6, abbreviate=TRUE,minlength=1), pch=0.3, 
             cex.lab=0.3, cex.title=0.3, auto.key = F, main = "Earth's average 
             temperature IDW map 2008", key.space=list(space="right", cex=0.8))

split.screen( rbind(c(0, 1,0,1), c(1,1,0,1)))
split.screen(c(1,2), screen=1)-> ind
screen( ind[1])
p1
screen( ind[2])
image.plot(legend.only=TRUE, legend.width=0.5, col=pal2(100), 
           smallplot=c(0.7,0.75, 0.3,0.7), zlim=c(min(idw.croatia@data),
                                                  max(idw.croatia@data)), axis.args = list(cex.axis = 0.7))
close.screen( all=TRUE)

## End(Not run)

