rm(list=ls())
library(spatstat) #para trabajar patrones puntuales

##########################################################################
# Tipos de patrones (aleatorio, uniforme, agregado), Patrones marcados,
# Covariables.
##########################################################################

## Aleatorio

aleatorio<-rpoispp(100)
x11()
plot(aleatorio, main="Aleatorio", pch=20)
summary(aleatorio)
plot(density(aleatorio), main="Estimación kernel de la intensidad")
contour(density(aleatorio), add=TRUE, main="Estimación kernel de la intensidad")
points(aleatorio, pch=20)
cx<-runif(600, 0,4)#coordenada x
cy<-runif(600,0,3) #coordenada y
#se generan 600 sitios aleatorios
patron<-ppp(cx,cy,c(0,4),c(0,3))
plot(patron, main="Aleatorio")
summary(patron) #lambda=50, porque el area=12 y 600/12=50
plot(density(patron),main="Intensidad Constante = 50")
contour(density(patron),add=TRUE)
plot(patron, add=TRUE, col=2)
plot(density(patron, bw=.25), main="Intensidad Constante = 50")
contour(density(patron, bw=.25),add=TRUE)
### No Homogéneo
lambda<-function(x,y)
{
100*(x+y)
}
nohomo<-rpoispp(lambda)
plot(nohomo)
### Proceso Matern Cluster
Mat<-rMatClust(kappa=10, r=0.1, mu=100)
Mat
plot(Mat, main="Agregado", pch=20)
summary(Mat)
plot(density(Mat), main="Agregado")
contour(density(Mat), add=TRUE, main="Estimación kernel de la intensidad")
points(Mat, pch=20)
### Patrón Regular
data(cells)
plot(cells, main="Regular", pch=20)


############################################################
## Simulando un patrón marcado (aleatorio y con tendencia)
############################################################

x<-runif(0,1,n=100)
y<-runif(0,1,n=100)
m1<-rnorm(100,10,2)
P<-ppp(x, y, c(-0.1,1.1),c(-.1,1.1), marks =m1)
plot(P, main="Patrón Marcado")
m2<-x-y+rnorm(1,0,1)
P2<-ppp(x, y, c(-0.1,1.1),c(-.1,1.1), marks =m2)
plot(P2, main="")

##########################################################################
# Estimación kernel de la densidad. Se usa un kernel Gaussiano
# K-Ripley
##########################################################################

############################
### Datos de pinos en Suecia
############################

demo(data)
data(swedishpines)
x<-(swedishpines)
plot(density(x,10), main="Estimación kernel de la intensidad")
contour(density(x,10), add=TRUE, main="Estimación kernel de la intensidad")
points(x, pch=20)
summary(x)

u<-Kest(x)
plot(u,, main="K-Ripley para los datos Swedishpines")
u

#####################################
### Datos de arboles en Nueva Zelanda
#####################################

data(nztrees)
y<-(nztrees)
plot(density(y,10), main="Estimación kernel de la intensidad")
contour(density(y,10), add=TRUE)
points(y, pch=20)
summary(y)

###########################################
### Datos de arboles en una selva tropical
############################################

data(bei)
z<-density(bei,sigma=70)
plot(z, main="Datos de árboles en una selva tropical")
plot(density(bei), main="Estimación de la intensidad")
#points(bei, pch=20)
contour(density(bei), add=TRUE)
contour(density(bei))

# Ejemplo de covariables
# Nota bei.extra$elev es formato im(pixel image)
elev<-bei.extra$elev
plot(elev, main="")
plot(bei, pch=20, main="Elevación como covariable", add=TRUE)
contour(density(bei), add=TRUE)
# Ejemplo de efecto de borde
plot(bei, add=TRUE)

##################################################################
# Métodos basados en cuadrantes
##################################################################

aleatorio<-rpoispp(100)
plot(aleatorio, main="Aleatorio", pch=20)
summary(aleatorio)
plot(density(aleatorio), main="Estimación kernel de la intensidad")
contour(density(aleatorio), add=TRUE, main="Estimación kernel de la intensidad")
points(aleatorio, pch=20)
Q<-quadratcount(aleatorio, nx=6,ny=3)
plot(aleatorio)
plot(Q, add=TRUE)
M<-quadrat.test(aleatorio, nx=6, ny=3)
plot(aleatorio)
plot(M, add=TRUE)
M

########################################################################
# Simulación de la estadistica X-squared para la prueba de aleatoriedad
# basado en el conteo de cuadrículas.
#########################################################################

aleatorio<-rpoispp(100)
Mobs<-quadrat.test(aleatorio, nx=3, ny=3)
plot(aleatorio)
plot(Mobs, add=TRUE)
Mobs
s<- 1000 # Número de simulaciones del proceso aleatorio
Xsq<-matrix(0, nrow=s, ncol=1)
for (i in 1:s)
{
aleatorio<-rpoispp(100)
M<-quadrat.test(aleatorio, nx=3, ny=3)
Xsq[i,]<-M$statistic
}
hist(Xsq, freq=FALSE) # Histograma de los datos simulados

x<-rchisq(10000, df=8)
curve(dchisq(x, df = 8), col = 2, lty = 2, lwd = 2, add = TRUE) # Curva teórica de la chi-cuadrado con 8 gl
lines(density(Xsq), col="blue") # Estimación kernel de la densidad
abline(v=Mobs$statistic, col=2, lty=2) # Valor de la estadística observada
# Lambda no homogéneo
lambda<-function(x,y)
{
100*(x+y)
}
nohomogeneo<-rpoispp(lambda)
nohomogeneo
MobsNohomoge<-quadrat.test(nohomogeneo, nx=3, ny=3)
abline(v=MobsNohomoge$statistic, col=3, lty=2)
# Nota: los grados de libertad de la chi cuadrado son [(3x3)-1]
# y entonces el percentil 95 (valor crítico de la prueba es
critico<-qchisq(0.95,8)
critico
# Estimación del alpha (Pr rechazar Ho siendo cierto) a partir de las simulaciones.
# Es decir en cuantas ocasiones de las 1000 la estadística observada es mayor que el crítico
suma<-0
for (i in 1:s)
{
if (Xsq[i,]> critico) suma<-suma+1 else suma<-suma
}
suma
alpha.estimado<-suma/s
alpha.estimado
# Probabilidad de rechazo con los datos simulados del proceso no homogéneo de media lambda
# Es decir en cuantas ocasiones de las 100 rechazo cuando el proceso es no homogéneo
sumanh<-0
XNoH<- matrix(0,nrow=s, ncol=1)
for (i in 1:s)
{
nohomogeneo<-rpoispp(lambda)

M2<-quadrat.test(nohomogeneo, nx=3, ny=3)
XNoH[i,]<-M2$statistic
if (XNoH[i,]> critico) sumanh<-sumanh+1 else sumanh<-sumanh
}
sumanh
probabilidad_rechazo<-sumanh/s
probabilidad_rechazo

#################################################
# Funciones G, F, K, L
#################################################

# Proceso completamente aleatorio
aleatorio<-rpoispp(100)
G1<-Gest(aleatorio)
plot(G1)
plot(envelope(aleatorio, Gest), main=" Función G con bandas de confianza")
F1<-Fest(aleatorio)
plot(F1)
plot(envelope(aleatorio, Fest), main=" Función F con bandas de confianza")
K1<-Kest(aleatorio)
plot(K1)
plot(envelope(aleatorio, Kest), main=" Funci´on K con bandas de confianza")
E <- envelope(aleatorio, Kest)
plot(E, sqrt(./pi) ~ r, main="Función L")
# Proceso de inhibición (regular)
data(cells)
plot(cells)
G2<-Gest(cells)
plot(G2)
plot(envelope(cells, Gest), main= "Función G con bandas de confianza")
F2<-Fest(cells)
plot(F2)
plot(envelope(cells, Fest), main= "Función F con bandas de confianza")
x11()

plot(envelope(cells, Fest), main=" Función F con bandas de confianza")
K2<-Kest(cells)
plot(K2)
plot(envelope(cells, Kest), main=" Función K con bandas de confianza")
E <- envelope(cells, Kest)
plot(E, sqrt(./pi) ~ r, main="Función L")

# Proceso de agregado
data(redwood)
plot(redwood)
G3<-Gest(redwood)
plot(G3)
plot(envelope(redwood, Gest), main= "Función G con bandas de confianza")
F3<-Fest(redwood)
plot(F3)
plot(envelope(redwood, Fest), main= "Función F con bandas de confianza")
K3<-Kest(redwood)
plot(K3)
plot(envelope(redwood, Kest), main=" Función K con bandas de confianza")
E <- envelope(redwood, Kest)
plot(E, sqrt(./pi) ~ r, main="Función L")

#################################################
# Simulacion de procesos
# Homogeneo, No homogeneo, Cluster, Cox, Matern
##################################################

#######################################
## Poisson Homogéneo
#######################################

aleatorio<-rpoispp(100)
plot(aleatorio, main="Aleatorio", pch=20)

summary(aleatorio)
plot(density(aleatorio), main="Estimación kernel de la intensidad")
contour(density(aleatorio), add=TRUE, main="Estimación kernel de la intensidad")
points(aleatorio, pch=20)
points(0.5, 0.8, col=2, pch=19)
fit_homo<-ppm(aleatorio, ~1)
names(fit_homo)
fit_homo$coef[[1]]
exp(fit_homo$coef[[1]]) # Estimación de la intensidad basada en el beta estimado.

##################################
# Proceso Poisson No Homogéneo
#################################

lambda<-function(x,y)
{
100*(x+y)
} #en la medida que xy y aumenten la
#funcion de intensidad tmb aumenta
nohomo<-rpoispp(lambda) #va ir cambiando dependiendo de las coordenadas x y
summary(nohomo)
windows()
plot(nohomo)
plot(density(nohomo), main ="Patrón No Homogéneo")
contour(density(nohomo), add=TRUE)

# Ajuste de un modelo lineal a la intensidad
fit_nohomo<-ppm(nohomo, ~x+y) #lambda(x,y)= exp(beta_0+beta1*x+beta2*y)
names(fit_nohomo)
fit_nohomo$coef
# Función para calcular la intensidad en un punto (x,y)
lambda_lineal<-function(x,y,fit)
{
lambda<-exp(fit[[1]]+fit[[2]]*x+fit[[3]]*y)
return(lambda)
}

x<-0.7
y<-0.8
lambda_lineal(x,y,fit_nohomo$coef) # Estimaci´on de la función de intensidad en x=0.5 y y=0.5

# Comparar contra el gráfico
par(mfrow=c(1,2))
plot(density(nohomo), main ="Patrón No Homogéneo. Ajuste No Paramétrico")
contour(density(nohomo), add=TRUE)
points(x,y, col=2, pch=19)
plot(fit_nohomo, se=FALSE)
points(x,y, col=2, pch=19)
# Ajuste de un modelo polinomial
fit_nohomo2<-ppm(nohomo, ~polynom(x,2)+y)
fit_nohomo2
plot(fit_nohomo2, se=FALSE)
fit_nohomo3<-ppm(nohomo, ~polynom(x,y,2))
fit_nohomo3
plot(fit_nohomo3, se=FALSE)

# Comparación de modelos

#(usando la tabla de Desviación basada en la razón de verosimilitud)
fit_nohomo<-ppm(nohomo, ~x+y)
fitnull<-update(fit_nohomo,~1)
anova(fitnull, fit_nohomo, test="Chi")
fit_nohomo2<-ppm(nohomo, ~polynom(x,2)+y)
fitnull2<-update(fit_nohomo2, ~x+y)
anova(fitnull2, fit_nohomo2, test="Chi")

# Interpretación: Si el p-valor es pequeño (p.ej menor del 5%)
# se rechaza la hipótesis nula (es decir que el modelo de
# menos parámetros es mejor)

lambda<-function(x,y)
{
200*(x+y+x^2+y^2+x*y)
}
nohomo<-rpoispp(lambda)
summary(nohomo)
plot(nohomo)

plot(density(nohomo), main ="Patrón No Homogéneo")
contour(density(nohomo), add=TRUE)
fit_nohomo3<-ppm(nohomo, ~polynom(x,y,2))
plot(fit_nohomo3, se=FALSE)
fitnull3<-update(fit_nohomo3, ~x+y)
anova(fitnull3, fit_nohomo3, test="Chi")

###############################################################
# Realizaciones de un proceso de Cox con intendidad constante
# lambda (lambda cambia en cada iteración)
##################################################################

par(mfrow=c(2,2))
vectorl<-c(rep(0,10000))
vectorn<-c(rep(0,10000))
for (i in 1:10000)
{
alpha<-10
beta<-5
lambda<-rgamma(1, shape=alpha, scale=beta)
vectorl[i]<-lambda
Cox<-rpoispp(lambda)
#plot(Cox)
#print(Cox)
vectorn[i]<-Cox$n
}
vectorl
vectorn
mean(vectorn)
mean(vectorl)
# Nota: Cuando se simula un proceso Poisson homogéneo,
# en cada iteración el lambda es el mismo,
# en el de Cox el lambda cambia!!!

################################
#Poisson no homogeneo aleatorio
################################

par(mfrow=c(2,2))
vectorl<-c(rep(0,100))
vectorn<-c(rep(0,100))

for (i in 1:2)
{
alpha<-10
beta<-5
lambda1<-rgamma(1, shape=alpha, scale=beta)
lambda<-function(x,y)
{
lambda1*(x+y)
}
Cox<-rpoispp(lambda)
plot(Cox)
plot(density(Cox), main ="Patrón No Homogéneo")
contour(density(Cox), add=TRUE)
print(Cox)
vectorn[i]<-Cox$n
}
mean(vectorn)
mean(vectorl)
# Cox Log_normal
par(mfrow=c(3,2))
for (i in 1:3)
{
lambda<-function(x,y)
{
exp(rnorm(1,5,1)*(x+y))
}
Cox_log<-rpoispp(lambda)
plot(Cox_log)
plot(density(Cox_log), main ="Patrón No Homogéneo")
contour(density(Cox_log), add=TRUE)
#print(Cox_log)
}
hist(e)

###############################################################
# Realizaciones de un proceso de Cox con intendidad NO constante
# lambda (lambda cambia en cada iteraci´on y depende de x,y)
##################################################################

# Poisson/Gamma
par(mfrow=c(1,3))
for (i in 1:3)
{
alpha<-10
beta<-5
lambda<-rgamma(1, shape=alpha, scale=beta)
lambda1<-function(x,y)
{
lambda*(x^2+y^2)
}
Cox_nohomo<-rpoispp(lambda1)
plot(Cox_nohomo)
}
# Cox Log_normal
par(mfrow=c(1,3))
for (i in 1:3)
{
lambda<-function(x,y)
{
exp(rnorm(1,5,1)*(x+y))
}
Cox_log<-rpoispp(lambda)
plot(Cox_log)
}

#############################################
# Proceso Poisson Agregado (cluster)
#############################################

# Proceso de Thomas
rThomas
pp_Thomas<-rThomas(10, 2, 10)
plot(pp_Thomas)
fit<-kppm(pp_Thomas, ~x+y, "Thomas")

fit # Estima los parámetros de la función K(h) usando el método del mínimo contraste
#ensayo para ver los clusters
pp_Thomas<-rThomas(2, 0.02, 100)
plot(pp_Thomas)
fit<-kppm(pp_Thomas, ~1, "Thomas")
fit

#funcion k de ripley para proceso cluster
K1<-Kest(pp_Thomas,correction="none")
plot(K1)
pi

#teorica
clase<-function(h)
{
k<-(pi*h^2)+(1-exp(-(h^2))/4)*(0.01^2)))/2
retorn(k)
}
clase(0.05)
dist=seq(0,12,by=0.01)
k_repley<-clase(dist)
plot(dist,k_ripley, type="l")


