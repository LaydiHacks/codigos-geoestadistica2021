        ########################################################################################
        ###########                    Ejercicio Montecarlo Prevosti               #############
        ########################################################################################

# Instrucción generación matrices cuando se ingresan valores de una parte triángular
matriz <- function(vect.tri){
n.filas<-function(vect.tri){k<-1
                            while (sum(seq(1:k))<=length(vect.tri))
                            k<-k+1
                            k
                            }
n.filas1<-n.filas(vect.tri)
X <- vect.tri
m.i <- matrix(0, n.filas1, n.filas1)
m.i[lower.tri(m.i)] <- NA
m.i[is.na(m.i)] <- X
ma.i <- t(m.i)+m.i
print(ma.i)
}

# FUNCIÓN MANTEL p.prueba
p.prueba <- function(B,matriz.sim, matriz.dist){                      #system.time({
set.seed(127)
nf <- nrow(matriz.sim)
low.tri <- lower.tri(matriz.dist)
prod.datos <- sum(matriz.sim[low.tri]*matriz.dist[low.tri])
conteo <- replicate(B, {
                    filas.perm <- sample(1:nf)
                    ifelse(sum(matriz.sim[low.tri]*matriz.dist[filas.perm,filas.perm][low.tri])>=prod.datos,1,0)
                    }
                    )
pvalor <- (sum(conteo)+1)/(length(conteo)+1)
pvalor
                                                                      #})
}

# Función Mantel p.Mantel (Propuesto)
p.Mantel <- function(B,matriz.sim, matriz.dist, Gráfico=F){           #system.time({
set.seed(127)
nf <- nrow(matriz.sim)
low.tri <- lower.tri(matriz.dist)
t.tilde <- sum(matriz.sim[low.tri]*matriz.dist[low.tri])
t.estrella<-numeric(B)
for (i in 1:B){
filas.perm <- sample(1:nf)
t.estrella[i]<-sum(matriz.sim[low.tri]*matriz.dist[filas.perm,filas.perm][low.tri])
}
pvalor <- (sum(t.estrella>=t.tilde)+1)/(B+1)

if (Gráfico) {
       plot(density(t.estrella),type="l", col=2,
       main="Distribución Aproximada del Estadístico de Mantel",
       xlab="t.estrella",ylab="Número de Permutaciones",
       sub=paste("Est. Observado (t.tilde)=",round(t.tilde,3),
        ":", "Cola Superior P-value=",round(pvalor,4),":", B, "Permutaciones"),
       lwd=1.5, col.main="blue", col.lab="blue", cex.lab=1.2, cex.sub=0.7)
       abline(v=t.tilde, lty=2, col=3)
   }
list(T=t.tilde, T.Simetrico=t.tilde/2, Cola.Superior.P=pvalor)
                                                                      #})
}
# Esta función es más eficiente en tiempo, y además ofrece la posibilidad de gráficar la función de densidad del test y de mostrar algunos resultados asociados al mismo.


         #######################################################################################
         ########################             Ingreso Datos           ##########################
         #######################################################################################

# La matrices de barreras consideran 1 si hay de por medio una barrera (mar o montaña) entre los dos lugares y cero si no hay barreras.

# Distancias Geográficas, Geneticas y de Barreras en el Estrecho de Gibraltar (no se considera Valencia)
Tabla5.Gen.GibSt<-c(0.104, 0.285, 0.425, 0.392, 0.443, 0.296, 0.415, 0.383, 0.428, 0.209, 0.215, 0.228, 0.123, 0.064, 0.133)
Tabla5.Dist.GibSt<-c(228, 160, 720, 810, 900, 193, 680, 720, 825, 550, 650, 750, 189, 182, 135)
Tabla5.Barr.GibSt<-c(0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)

Tabla5.Gen.GS<-matriz(Tabla5.Gen.GibSt)
Tabla5.Dist.GS<-matriz(Tabla5.Dist.GibSt)
Tabla5.Barr.GS<-matriz(Tabla5.Barr.GibSt)


# Distancias Geográficas, Geneticas y de Barreras en el Área Tyrrhenian
Tabla6.Gen.Tyrr<-c(0.125, 0.465, 0.452, 0.456, 0.322, 0.451, 0.569, 0.573, 0.634, 0.465, 0.457, 0.463, 0.327, 0.448, 0.556, 0.560, 0.619, 0.080, 0.129, 0.297, 0.225, 0.516, 0.536, 0.544, 0.099, 0.262, 0.177, 0.478, 0.513, 0.527, 0.291, 0.203, 0.524, 0.539, 0.557, 0.184, 0.504, 0.509, 0.533, 0.515, 0.527, 0.550, 0.063, 0.152, 0.161)
Tabla6.Dist.Tyrr<-c(207, 445, 392, 413, 285, 407, 575, 653, 879, 598, 555, 508, 245, 265, 585, 693, 846, 63, 168, 458, 645, 443, 418, 762, 126, 400, 585, 420, 400, 678, 325, 503, 282, 282, 563, 192, 332, 443, 588, 430, 560, 601, 130, 340, 345)
Tabla6.Barr.Tyrr<-c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0)


Tabla6.Gen.Ty<-matriz(Tabla6.Gen.Tyrr)
Tabla6.Dist.Ty<-matriz(Tabla6.Dist.Tyrr)
Tabla6.Barr.Ty<-matriz(Tabla6.Barr.Tyrr)


# Distancias Geográficas, Geneticas y de Barreras en el Estrecho de Dover
Tabla7.Gen.SDover<-c(0.307, 0.307, 0.152, 0.271, 0.260, 0.235, 0.083, 0.290, 0.263, 0.399, 0.331, 0.276, 0.225, 0.370, 0.300, 0.150, 0.187, 0.112, 0.195, 0.120, 0.128)
Tabla7.Dist.SDover<-c(913, 913, 760, 1354, 1326, 1373, 16, 677, 906, 1558, 1222, 689, 921, 1571, 1237, 601, 884, 663, 1010, 450, 591)
Tabla7.Barr.SDover<-c(1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)

Tabla7.Gen.SD<-matriz(Tabla7.Gen.SDover)
Tabla7.Dist.SD<-matriz(Tabla7.Dist.SDover)
Tabla7.Barr.SD<-matriz(Tabla7.Barr.SDover)


# Distancias Geográficas, Geneticas y de Barreras en el Dardanelles y Bosphorus
Tabla8.Gen.DB<-c(0.119, 0.162, 0.322, 0.248, 0.305, 0.306, 0.200, 0.265, 0.232, 0.231, 0.268, 0.366, 0.288, 0.340, 0.312, 0.110, 0.120, 0.136, 0.138, 0.098, 0.154)
Tabla8.Dist.DB<-c(273, 162, 865, 647, 1061, 924, 272, 632, 510, 893, 878, 789, 519, 1066, 745, 395, 296, 514, 603, 267, 594)
Tabla8.Barr.DB<-c(0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)

Tabla8.Gen.DB<-matriz(Tabla8.Gen.DB)
Tabla8.Dist.DB<-matriz(Tabla8.Dist.DB)
Tabla8.Barr.DB<-matriz(Tabla8.Barr.DB)


# Distancias Geográficas, Geneticas y de Barreras en el Rango Pyrenean
Tabla9.Gen.Pyrenean<-c(0.130, 0.251, 0.262, 0.223, 0.229, 0.136)
Tabla9.Dist.Pyrenean<-c(117, 197, 281, 89, 189, 104)
Tabla9.Barr.Pyrenean<-c(0, 1, 1, 0, 1, 0)

Tabla9.Gen.P<-matriz(Tabla9.Gen.Pyrenean)
Tabla9.Dist.P<-matriz(Tabla9.Dist.Pyrenean)
Tabla9.Barr.P<-matriz(Tabla9.Barr.Pyrenean)



# Test de Mantel
p.prueba(10000, Tabla5.Gen.GS,Tabla5.Dist.GS)
p.prueba(10000, Tabla6.Gen.Ty, Tabla6.Dist.Ty)
p.prueba(10000, Tabla7.Gen.SD, Tabla7.Dist.SD)
p.prueba(10000, Tabla8.Gen.DB, Tabla8.Dist.DB)
p.prueba(10000, Tabla9.Gen.P, Tabla9.Dist.P)

# Dado que los valores p.prueba son menores del nivel de significancia 5% (en todos los casos),
# se tiene que las matrices de distancias genéticas y geográficas estan correlacionadas, dado
# que no hay evidencia estadística para no rechazar "Ho: No hay relación entre elementos en
# las matrices de distancias genéticas y geográficas".

p.prueba(10000, Tabla5.Gen.GS,Tabla5.Dist.GS)
p.prueba(10000, Tabla6.Gen.Ty, Tabla6.Barr.Ty)
p.prueba(10000, Tabla7.Gen.SD, Tabla7.Barr.SD)
p.prueba(10000, Tabla8.Gen.DB, Tabla8.Barr.DB)
p.prueba(10000, Tabla9.Gen.P, Tabla9.Barr.P)

# Los valores p.prueba son menores del nivel de significancia 5%, para los casos de las
# Tablas 6, 7 y 8. Y en los casos; Pirineos y Estrecho de Gibraltar el p.prueba es superior
# al 5%, estando muy cerca del mismo valor. Por lo cuál las matrices de distancias genéticas
# y de barreras estan correlacionadas en los casos de las Tablas 6, 7 y 8, encontrando
# evidencia estadística para rechazar "Ho: No hay relación entre elementos en las matrices de
# distancias genéticas y de barreras". Y en los casos del Estrecho de Gibraltar y de los
# Pirineos no se rechaza Ho.


# Solución a partir de la función para el Test de Mantel "propuesto"
par(mfrow = c(4,3))
p.Mantel(10000, Tabla5.Gen.GS,Tabla5.Dist.GS, Gráfico=T)
p.Mantel(10000, Tabla6.Gen.Ty,Tabla6.Dist.Ty, Gráfico=T)
p.Mantel(10000, Tabla7.Gen.SD,Tabla7.Dist.SD, Gráfico=T)
p.Mantel(10000, Tabla8.Gen.DB,Tabla8.Dist.DB, Gráfico=T)
p.Mantel(10000, Tabla9.Gen.P,Tabla9.Dist.P, Gráfico=T)

p.Mantel(10000, Tabla5.Gen.GS,Tabla5.Barr.GS, Gráfico=T)
p.Mantel(10000, Tabla6.Gen.Ty,Tabla6.Barr.Ty, Gráfico=T)
p.Mantel(10000, Tabla7.Gen.SD,Tabla7.Barr.SD, Gráfico=T)
p.Mantel(1000, Tabla8.Gen.DB,Tabla8.Barr.DB, Gráfico=T)
p.Mantel(10000, Tabla9.Gen.P,Tabla9.Barr.P, Gráfico=T)

# Lógicamente las conclusiones son las mismas ya que esta función considera tambien los productos, la diferencia es que el algoritmo es más eficiente en tiempo.



  #######################################################################################################
  #                                        Breves Comentarios:                                          #
  #######################################################################################################

# Si se encontro relación entre las matrices de distancias geográficas y las de distancias
# genéticas, dado que en todos los casos se encontró evidencia estadística para rechazar
# Ho: No hay relación entre los elementos de las matrices de distancias genéticas y
# geográficas.


# El test de Mantel efectivamente es un buen instrumento estadístico para estos casos de
# patrones espaciales. Encuentro que en el caso de las barreras geográficas, el test
# corrobora las afirmaciones del articulo. De hecho para el caso de las barreras geográficas
# el p-value del test de Mantel es superior al 5%, no rechazando la hipótesis nula en los
# casos de las tablas 5 y 9 (Estrecho de Gibraltar y Pirineos), esto confirma lo dicho al
# final del articulo en cuanto a la agudeza de las diferencias de las distancias:

#"En casos tales como los Pirineos y el Estrecho de Gibraltar, en los cuales la barrera que
# separa las dos areas continentales, la agudeza en las diferencias de las poblaciones de
# ambos lados es probablemente debido a una tasa asimétrica del flujo de genes. En estas
# poblaciones el flujo genético de un lado de la barrera es muy limitado considerando que
# éste esta abierto al otro lado" .


# En los Pirineos no se rechaza Ho: No hay relación entre elementos en las matrices de
# distancias genéticas y de barreras, al nivel de significancia del 5% el p-valor=0.08,
# es decir la barrera montañosa no incide en las relaciones de las distancias genéticas
# (es posible también que por ser la muestra muy pequeña el resultado no sea muy confiable)
# y en el Estrecho de Gibraltar tambien no se rechaza Ho al 5% dado que el p-valor=0.068.
# En cambio en los demas casos si se encuentra relación entre las distancias genéticas y las
# barreras geográficas, basicamente los diferentes estrechos que hay de por medio entre dos
# regiones geográficas o el mar mediterraneo, para los casos en que este éste de por medio.

# En algunos apartados del articulo se menciona la similitud de las distancias genéticas
# entre las poblaciones de Drosophilia suboscura cuando se encuentran en las islas o en el
# caso de encontrarse en los continentes y las diferencias entre las ubicadas en los
# continentes con respecto a las ubicadas en las islas. Esto se confirma en parte con lo
# realizado a partir del test de Mantel, pero para confirmarlo con mayor certeza seria
# necesario considerar solo los datos entre las islas sin considerar los datos de las
# poblaciones ubicadas en los continentes.

# Se puede agregar que el articulo es una evidencia más de la primera ley geográfica que dice
# que las cosas que estan más cerca son más similares y las que estan más distances más
# disimilares.



                             ##############################################
                             # Trabajo con library ecodist test de Mantel #
                             ##############################################


# Esta libreria tambien permite hacer el test de Mantel, trabaja internamente las dos posibilidades de ">=" y "<=", para relaciones directas e inversas respectivamente, y además ofrece alternativas para la hipótesis nula (en éste ejercicio sería r=0). Finalmente, la instrucción mgram permite obtener un mantelgrama, algo así como un correlograma con el test de Mantel, y este relaciona las diferentes distancias con los valores asociados de Mantel. La instrucción en éste caso es "mantel" y dentro se referencian las dos matrices a comparar, luego se establece el número de permutaciones, es importante que estas matrices esten como distancias (es decir, en éste caso se usaría: as.dist().)
# Los resultados obtenidos mediante esta instrucción corroboran los resultados obtenidos anteriormente. Y con los del ejercicio que ya habiamos realizado con los datos de la deriva continental, haciendo la observación que en dicho caso se trabajaba con "<=" (debido a la relación inversa) y con esta función "mantel", de la libreria "ecodist" es más simple.

library(ecodist)                                 # Disponible para versión: R.2.7.0
A5<-as.dist(Tabla5.Gen.GS)
B5<-as.dist(Tabla5.Dist.GS)
C5<-as.dist(Tabla5.Barr.GS)
t.est.T5<-mantel(A5~B5,nperm=10000)                # Ha: r<=0 (pval1), r>=0 (pval2), r=0 (pval3)
t.est1.T5<-mantel(A5~C5,nperm=10000)
t.mgram.T5<-mgram(A5, B5, nperm=0)
plot(t.mgram.T5)

A6<-as.dist(Tabla6.Gen.Ty)
B6<-as.dist(Tabla6.Dist.Ty)
C6<-as.dist(Tabla6.Barr.Ty)
t.est.T6<-mantel(A6~B6,nperm=10000)                # Ho: r<=0 (pval1), r>=0 (pval2), r=0 (pval2)
t.est1.T6<-mantel(A6~C6,nperm=10000)
t.mgram.T6<-mgram(A6, B6, nperm=0)
plot(t.mgram.T6)

A7<-as.dist(Tabla7.Gen.SD)
B7<-as.dist(Tabla7.Dist.SD)
C7<-as.dist(Tabla7.Barr.SD)
t.est.T7<-mantel(A7~B7,nperm=10000)               # Ho: r<=0 (pval1), r>=0 (pval2), r=0 (pval2)
t.est1.T7<-mantel(A7~C7,nperm=10000)
t.mgram.T7<-mgram(A7, B7, nperm=0)
plot(t.mgram.T7)

A8<-as.dist(Tabla8.Gen.DB)
B8<-as.dist(Tabla8.Dist.DB)
C8<-as.dist(Tabla8.Barr.DB)
t.est.T8<-mantel(A8~B8,nperm=10000)               # Ho: r<=0 (pval1), r>=0 (pval2), r=0 (pval2)
t.est1.T8<-mantel(A8~C8,nperm=10000)
t.mgram.T8<-mgram(A8, B8, nperm=0)
plot(t.mgram.T8)

A9<-as.dist(Tabla9.Gen.P)
B9<-as.dist(Tabla9.Dist.P)
C9<-as.dist(Tabla9.Barr.P)
t.est.T9<-mantel(A9~B9,nperm=10000)               # Ho: r<=0 (pval1), r>=0 (pval2), r=0 (pval2)
t.est1.T9<-mantel(A9~C9,nperm=10000)
t.mgram.T9<-mgram(A9, B9, nperm=0)
plot(t.mgram.T9)


# Situación Deriva Continental
A<-as.dist(m.earw)
B<-as.dist(m.dist1)
C<-as.dist(m.dist2)
t.est<-mantel(A~B,nperm=10000)                  # Ho: r<=0 (pval1), r>=0 (pval2), r=0 (pval2)
t.est1<-mantel(A~C,nperm=10000)
t.mgram<-mgram(A, B, nperm=0)
plot(t.mgram)

library(ape)
set.seed(127)
mantel.test(Tabla5.Gen.GS,Tabla5.Barr.GS, graph = TRUE)

