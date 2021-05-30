#########################################################
############    EJERCICIO GWR "COLUMBUS"     ############
#########################################################
# https://rpubs.com/bogotan/AMESP12ModelosDatosareas
# http://homepage.ntu.edu.tw/~wenthung/R_Spatial/RLab_9.html
# clear the workspace, plots and console
rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014")


library(rgdal)
library (spdep)
library(spgwr)
library (GWmodel) 


require(maptools)
columbus <- readOGR(system.file("shapes/columbus.shp", package="spData")[1])
col.poly <- as.data.frame(columbus)
o.nb <- read.gal("db/columbus.gal")    
a.lw <- nb2listw(o.nb, style="W")
col.poly$WX <- lag.listw(a.lw, col.poly$X)
col.poly$WY <- lag.listw(a.lw, col.poly$Y)

#  Mapping INCOME
lm.palette <- colorRampPalette(c("white","orange", "red"), space = "rgb")
spplot(columbus, zcol="INC", col.regions=rev(lm.palette(20)), main="INCOME")

#  Mapping CRIME
spplot(columbus["CRIME"], col.regions=rev(lm.palette(20)), main="CRIME")
#  Mapping HOVAL
spplot(columbus["HOVAL"], col.regions=rev(lm.palette(20)), main="HOVAL")

# GW Summary: Geographically weighted summary statistics (GWSS)
# This function also calculates basic geographically weighted
# covariances together with basic and robust geographically weighted correlations.
# 1. Geographically-Weighted (GW) summary Statistics
# Mean, Variance and Correlation

GWS<-gwss(columbus, vars = c("HOVAL", "INC","CRIME"), adaptive = TRUE, bw = 10)
head(data.frame(GWS$SDF))
# GW Summary: Monte Carlo (randomization) tests
GWS.Sim<-montecarlo.gwss(columbus, vars = c("HOVAL", "INC","CRIME"), adaptive = TRUE, bw = 10)

# Mapping GW CRIME Mean
lm.palette <- colorRampPalette(c("white","orange", "red"), space = "rgb")
spplot(GWS$SDF, zcol="CRIME_LM", col.regions=rev(lm.palette(20)), main="GW CRIME")
# Mapping GW HOVAL Mean
spplot(col.poly, zcol="HOVAL", col.regions=rev(lm.palette(20)), main="HOVAL")
spplot(GWS$SDF, zcol="HOVAL_LM", col.regions=rev(lm.palette(20)), main="GW HOVAL")
# Mapping GW HOVAL p-value
col.poly$GW.Hoval.p<-GWS.Sim[,1]
spplot(col.poly, zcol="GW.Hoval.p", col.regions=lm.palette(20), main="p-value for GW Hoval")
# Mapping GW INCOME Mean
spplot(GWS$SDF, zcol="INC_LM", col.regions=rev(lm.palette(20)), main="GW INCOME")
# Mapping GW INCOME p-value
col.poly$GW.income.p<-GWS.Sim[,2]
lm.palette <- colorRampPalette(c("red","orange", "white"), space = "rgb")
spplot(col.poly, zcol="GW.income.p", col.regions=lm.palette(20), main="p-value for GW Income")

#  Mapping CRIME
spplot(col.poly, zcol="CRIME", col.regions=rev(lm.palette(20)), main="CRIME")
spplot(GWS$SDF, zcol="CRIME_LM", col.regions=rev(lm.palette(20)), main="GW CRIME")

# Mapping GW CRIME p-value
col.poly$GW.crime.p <- GWS.Sim[,3]
spplot(col.poly, zcol="GW.crime.p", col.regions=lm.palette(20), main="p-value for GW Crime")

# Spatial Non-stationary
spplot(GWS$SDF, zcol="Spearman_rho_HOVAL.CRIME", col.regions=lm.palette(20), main="GW Correlation")
spplot(GWS$SDF, zcol="Spearman_rho_INC.CRIME", col.regions=lm.palette(20), main="GW Correlation")
spplot(GWS$SDF, zcol="Spearman_rho_HOVAL.INC", col.regions=lm.palette(20), main="GW Correlation")

col.poly$GW.Corr.HI.p <- GWS.Sim[,16]
lm.palette <- colorRampPalette(c("red","orange", "white"), space = "rgb")
spplot(col.poly, zcol="GW.Corr.HI.p", col.regions=lm.palette(20), main="p-value for Corr. between GW Hoval and GW Income")

col.poly$GW.Corr.HC.p <- GWS.Sim[,17]
spplot(col.poly, zcol="GW.Corr.HC.p", col.regions=lm.palette(20), main="p-value for Corr. between GW Hoval and GW Crime")

col.poly$GW.Corr.IC.p <- GWS.Sim[,18]
spplot(col.poly, zcol="GW.Corr.IC.p", col.regions=lm.palette(20), main="p-value for Corr. between GW Income and GW Crime")

# 2. GWR Model Specification
# OLS Model

lm.global <- lm(CRIME ~ INC + HOVAL, data=col.poly)
summary(lm.global)

#GW Model
# Crossvalidation of bandwidth for geographically weighted regression
adapt <- gwr.sel(CRIME~ INC + HOVAL + X + Y + WX + WY, data=columbus, coords=cbind(col.poly1$X,col.poly1$Y)) 


#Define Variables DeV: Var Dep  e  InV: Var Indep
DeV<-"CRIME"
InV<-c("INC","HOVAL","X","Y","WX","WY")

#Model Selection
coordinates(columbus) <- c("X","Y")
model.sel <- model.selection.gwr(DeV,InV, data = col.poly1, kernel = "bisquare", adaptive = TRUE, bw = 10)

data.frame(model.sel[2])

sorted.models <- model.sort.gwr(model.sel, numVars = length(InV), ruler.vector = model.sel[[2]][,2])
model.list <- sorted.models[[1]]
x11()
str(model.view.gwr(DeV, InV, model.list = model.list))
model.view.gwr(DeV, InV, model.list = model.list)

plot(sorted.models[[2]][,2], col = "black", pch = 20, lty = 5, main = "Alternative view of GWR model selection procedure", ylab = "AICc", xlab = "Model number", type = "b")

# Optimal Bandwidth
bw.gwr.1 <- bw.gwr(CRIME~HOVAL+INC, data = columbus, approach = "AICc", kernel = "bisquare", adaptive = TRUE)
bw.gwr.1

# 3. Fitting a GWR model
# Understanding GWR outputs

gwr.result<-gwr.basic(CRIME~HOVAL+INC, data=columbus, kernel="bisquare", adaptive=TRUE, bw=bw.gwr.1)
print(gwr.result)

head(data.frame(gwr.result$SDF),n=49)

Table<-data.frame(gwr.result$SDF)
hist(Table$yhat, breaks=20)
hist(Table$HOVAL, breaks=20)

# Mapping the model parameters
lm.palette <- colorRampPalette(c("blue","orange", "red"), space = "rgb")
spplot(gwr.result$SDF, zcol="INC", col.regions=lm.palette(20), main="GWR Results INC_coef")

# Exporting to Shape files
ogrListLayers(dsn)
drv="ESRI Shapefile"
getwd()
writeOGR(gwr.result$SDF, dsn="ESRI Shapefile", layer="CRIME_GWR", driver=drv)
## Warning in writeOGR(gwr.result$SDF, dsn = "shapefiles/GWR", layer =
## "TB_GWR", : Field names abbreviated for ESRI Shapefile driver

#Mapping p-values from the pseudo t-tests of GWR outputs
#Monte Carlo test for significance of GWR parameter variability

gwr.sim.result<-montecarlo.gwr(CRIME~HOVAL+INC, data=col.poly1, kernel="bisquare", adaptive=TRUE, bw=bw.gwr.1)

# p-values from the pseudo t-tests of GWR outputs
pvalue<-gwr.t.adjust(gwr.result) 
pvalueTable<-pvalue$SDF@data
names(pvalueTable)

# Mapping pseudo t-statistic and p-value

# Mapping t-statistic
lm.palette <- colorRampPalette(c("white","orange", "red"), space = "rgb")
spplot(pvalue$SDF, zcol="HOVAL_t", col.regions=rev(lm.palette(20)), main="t-statistic of Hoval")
spplot(pvalue$SDF, zcol="INC_t", col.regions=rev(lm.palette(20)), main="t-statistic of Income")

# Mapping p-value
lm.palette <- colorRampPalette(c("green","yellow", "white"), space = "rgb")
spplot(pvalue$SDF, zcol="HOVAL_p", col.regions=lm.palette(20), main="p-value of Hoval")
spplot(pvalue$SDF, zcol="INC_p", col.regions=lm.palette(20), main="p-value of Income")

#spplot(pvalue$SDF, zcol=c("Intercept_t","HOVAL_t","INC_t"), col.regions=lm.palette(20), main="")
#spplot(pvalue$SDF, zcol=c("Intercept_p","HOVAL_p","INC_p"), col.regions=lm.palette(20), main="")

# Mapping p-value (Bonferroni)
lm.palette <- colorRampPalette(c("red","orange", "white"), space = "rgb")
spplot(pvalue$SDF, zcol="HOVAL_p_bo", col.regions=lm.palette(20), main="p-value (Bonferroni) of Hoval")
spplot(pvalue$SDF, zcol="INC_p_bo", col.regions=lm.palette(20), main="p-value (Bonferroni) of Hoval")



