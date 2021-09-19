library(geospt)
library(RANN)
library("plotrix")
data(preci)
So <- data.frame(t(c(4,2)))
names(So) <- c("x","y")
knn <- nn2(preci[,2:3],So,k=6)

plot(preci[,2:3],pch=19,ylim=c(-0.6,4.5))
draw.circle(4, 2, 2.3, nv = 1000, border = 2, col = NA, lty = 1, lwd = 1.5)
text(preci$x+0.12,preci$y+0.12,preci$Obs,cex=1)
points(So[1],So[2],col=2,pch=19)
text(So[1]+0.12,So[2]+0.12,"So",col=2,cex=1.1)

PP1 <- data.frame(preci[sort(knn$nn.idx),1:3],Dist=knn$nn.dists[order(knn$nn.idx)])
intensidad <- 6/(pi*2.3^2)
LHS <-  3/(pi*2.3^2)
PP1$RHS <- (1-(PP1$Dist/2.3)^2)^2
PP1$LHS.RHS <- LHS*PP1$RHS
PP1
KE <- sum(PP1$LHS.RHS)                          # Intensidad kernel cuártico

set.seed(123)
sample(10)


x = runif(100, -1, 1)
y = runif(100, -1, 1)
plot(x, y, asp = 1, xlim = c(-1, 1))
draw.circle(0, 0, 1, nv = 1000, border = NULL, col = NA, lty = 1, lwd = 1)

