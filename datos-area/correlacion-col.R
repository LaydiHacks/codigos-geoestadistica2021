#Resultados generados por el Profesor

load("correlación-espacial/Col")


par(mar = c(0, 0, 0, 0), pty = "s")
plot(col.map)

# Efecto reina orden 1
plot(col_nbq1, coords, add=T, col=4, lwd=2)
# Efecto reina orden 2
plot(col.lags[[2]], coords, add=T, col=4, lwd=2)
# Efecto reina orden 3
plot(col.lags3[[3]], coords, add=T, col=4, lwd=2)


# Correlograma I de Moran a partir de matriz contiguidad espacial                                                  
plot(cor.s) 

"En este caso estamos viendo el correlagrama para diferentes ordenes, 
lo que vamos a analizar aquí es si ese correlagrama de un determinado orden toca el cero, si toca el cero
la hipotesis nula va ser: hay ausencia de autocorrelación espacial. Por lo tanto si toca el cero no hay lio
pero en los casos donde toca el cero podemos decir que hay autocorrelacion negativa o positiva
"

# Vecinos más Cercanos
"Calcula cual es el vecino mas cercano"
plot(col.map)
plot(col_kn1, coords, add=T, col=4, lwd=2) 

plot(col.map)
plot(col_kn4, coords, add=T, col=4, lwd=2) #4 vecinos más cercanos


#Umbral de distancias
par(mar = c(0, 0, 0, 0), pty = "s")
plot(col.map)
plot(col_kd1, coords, add=T, col=4, lwd=2)


plot(corD$res,main="")
corD$res 
"En este caso podemos ver es significativa el 2 porque la probabilidad 0.009 es menor del 5%
es decir, en es variable hay autocorrelacion espacial"

par(mar=c(0,0,0,0))                           
quantile.e(x=col.map$A4, ic=5, digits=1, color="Spectral",style="quantile", border=col.map, north.arrow=c(70000,600000,1720000,0.7),scale.bar=list(470000,100000,500000,"m",0.7),leg.cex=0.8)



