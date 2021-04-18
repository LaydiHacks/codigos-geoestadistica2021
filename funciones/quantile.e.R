assign("quantile.e",  
function(x,ic,digits=1,style="quantile",border=border,size,xn,yn,xsb,ysb,length,unit){ 
plotvar <- round(x,digits=digits)    
nclr <- ic
plotclr <- brewer.pal(nclr,"YlOrBr")
class <- classIntervals(plotvar, nclr, style=style)
colcode <- findColours(class, plotclr)
plot(border)   
plot(border, col=colcode, add=T)
legend(locator(1), legend=names(attr(colcode, "table")), fill=attr(colcode, "palette"), cex=0.8, bty="n")
northarrow(c(xn,yn),size=size,cex=0.8)
scalebar(c(xsb,ysb),length=length,unit=unit,division.cex=0.6)
}
)
  
  
 