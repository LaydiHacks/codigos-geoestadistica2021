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


