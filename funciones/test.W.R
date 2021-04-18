#' Function to compute and test eigenvectors of spatial weighting matrices
#' 
#' This function is a user-friendly way to compute and test eigenvectors for 
#' various definitions of spatial weighting matrices. It combines calls to the 
#' functions \code{scores.listw} and \code{ortho.AIC}. It allows to test various
#' definitions of the spatial weighting matrix and return results of 
#' \code{scores.listw} for the best one.
#' 
#' @details This functions allows to test one binary spatial weighting matrix 
#'   (if only Y and nb are provided). It allows also to test a weighting 
#'   function based on distances (if f is provided) and a weighting function 
#'   with different values of parameters if other arguments of \code{f} are 
#'   provided.
#'   
#' @param Y A matrix with response variables (univariate or multivariate 
#'   response)
#' @param nb An object of the class \code{nb} created by functions of the 
#'   \code{spdep} package
#' @param xy Coordinates of the samples, this argument is optional and is 
#'   required only if the argument \code{f} is not null.
#' @param MEM.autocor A string indicating if all MEM must be returned or only 
#'   those corresponding to positive or negative autocorrelation
#' @param f A function of the distance that can be used as a weighting spatial 
#'   function. This argument is optional
#' @param \dots Others arguments for the function \code{f}. It defines the range
#'   of parameters which will be tested
#' @return A list with the following elements: \item{all }{A data.frame where 
#'   each row correspond to one spatial weighint matrix tested. It contains 
#'   value of parameteres tested and corrected AIC and number of orthogonal 
#'   vectors for the best model.} \item{best }{A list containing results of 
#'   scores.listw and ortho.AIC of the best spatial weighting matrix according 
#'   to corrected AIC.}
#' @author Stéphane Dray \email{stephane.dray@@univ-lyon1.fr}
#' @seealso   \code{\link{ortho.AIC}}, \code{\link{scores.listw}}
#' @references Dray S., Legendre P. and Peres-Neto P. R. (2006) Spatial 
#'   modeling: a comprehensive framework for principal coordinate analysis of 
#'   neighbor matrices (PCNM). Ecological Modelling, 196, 483--493
#' @keywords spatial
#' @examples
#' 
#' if(require(ade4) & require(spdep)){
#' 
#' data(oribatid)
#' # Hellinger transformation
#' fau <- sqrt(oribatid$fau / outer(apply(oribatid$fau, 1, sum), rep(1, ncol(oribatid$fau)), "*"))
#' # remove gradient effect
#' faudt <- resid(lm(as.matrix(fau) ~ as.matrix(oribatid$xy)))
#' 
#' # test a binary spatial weighting matrix
#' nbtri <- tri2nb(as.matrix(oribatid$xy))
#' tri.res <- test.W(faudt, nbtri)
#' 
#' maxi <- max(unlist(nbdists(nbtri, as.matrix(oribatid$xy))))
#' 
#' # test a simple spatial weighting function of the distance
#' f1 <- function(x) {1-(x)/(maxi)}
#' tri.f1 <- test.W(faudt, nbtri, f = f1, xy = as.matrix(oribatid$xy))
#' 
#' # test a spatial weighting function with various values of parameters
#' f2 <- function(x,dmax,y) {1-(x^y)/(dmax)^y}
#' tri.f2 <- test.W(faudt,nbtri, f = f2, y = 2:10, dmax = maxi, xy = as.matrix(oribatid$xy))
#' }
#' 
#' @importFrom spdep nb2listw
#' @export
#' 
assign("test.W", 
       function(Y,nb,xy,MEM.autocor = c("all", "positive", "negative"),f = NULL,...) {
         mycall <- pairlist(...)
         res <- list()
         MEM.autocor <- match.arg(MEM.autocor)
         if (!(is.null(f))) {
           nbdist <- nbdists(nb, as.matrix(xy))
           if (!(is.null(mycall))) {
             param <- expand.grid(as.list(mycall))
             m1 <- match(names(param), names(formals(f)))
             for (i in 1:nrow(param)) {
               formals(f)[m1] <- unclass(param[i,])
               res[[i]] <-
                 scores.listw(nb2listw(nb,style = "B",glist = lapply(nbdist, f),zero.policy = TRUE),MEM.autocor = MEM.autocor)
             }
           }
           else {
             res[[1]] <-scores.listw(nb2listw(nb, style = "B", glist = lapply(nbdist, f)),MEM.autocor = MEM.autocor)
           }
         }
         else {
           res[[1]] <-
             scores.listw(nb2listw(nb, style = "B"), MEM.autocor = MEM.autocor)
         }
         res2 <-
           lapply(res, function(x)
             ortho.AIC(Y = Y,X = x,ord.var = TRUE)
           )
         if (!(is.null(mycall))) {
           res3 <-
             data.frame(AICc = unlist(lapply(res2, function(x)
               min(x[[1]], na.rm = TRUE))), NbVar = unlist(lapply(res2, function(x)
                 which.min(x[[1]]))))
           res3 <- cbind(param, res3)
         }
         else{
           res3 <-
             data.frame(AICc = unlist(lapply(res2, function(x)
               min(x[[1]], na.rm = TRUE))), NbVar = unlist(lapply(res2, function(x)
                 which.min(x[[1]]))))
         }
         
         thebest <- which.min(res3$AICc)
         cat (paste("\n\nAICc for the null model:", res2[[thebest]]$AICc0, "\n"))
         cat ("\nBest spatial model:\n")
         print(res3[thebest,])
         
         return(list(all = res3, best = list(MEM = res[[thebest]], AIC = res2[[thebest]])))
       }
)

#' Compute AIC for models with orthonormal explanatory variables
#' 
#' This function compute corrected AIC for models with orthonormal and centered 
#' explanatory variables such as MEM spatial eigenfunctions. Variables are
#' sorted by their contribution to R2.
#' 
#' It ensures that a model with k variables is the best one that can be 
#' obtained. By default, response variables are centered (model with intercept).
#' 
#' 
#' @param Y A matrix with response variables (univariate or multivariate 
#'   response)
#' @param X A set of orthonormal and centered vectors
#' @param ord.var A logical value indicating if the order of variables and 
#'   cumulative R2 must be returned
#' @return A vector with corrected AIC if \code{ord.var=FALSE}. A list if
#'   \code{ord.var=TRUE} with: \item{AICc }{Values of corrected AIC.}
#'   \item{AICc0 }{Values of corrected AIC for the null model (only intercept).}
#'   \item{ord }{Order of variables to be enter in the model} \item{R2
#'   }{Cumulative R2}
#' @author Stéphane Dray \email{stephane.dray@@univ-lyon1.fr}
#' @references Godinez-Dominguez E. and Freire J. (2003) Information-theoretic 
#'   approach for selection of spatial and temporal models of community 
#'   organization. Marine Ecology - Progress Series. 253, 17--24
#' @keywords models
#' @examples
#' 
#' y <- matrix(rnorm(50),50,1)
#' x <- svd(scale(y \%*\% c(0.1,0.5,2,0,0.7)+matrix(rnorm(250),50,5)))$u
#' res <- ortho.AIC(y,x,ord.var=TRUE)
#' minAIC <- which.min(res$AICc)
#' nvar <- length(1:minAIC)+1 # number of orthogonal vectors + 1 for intercept
#' lm1 <- lm(y~x[,res$ord[1:minAIC]])
#' summary(lm1)$r.squared # R2
#' res$R2[minAIC] # the same

#' min(res$AICc) # corrected AIC
#' extractAIC(lm1) # classical AIC
#' min(res$AICc)-2*(nvar*(nvar+1))/(nrow(x)-nvar-1) # the same
#'
#' lm2 <- lm(y~1)
#'
#' res$AICc0 # corrected AIC for the null model
#' extractAIC(lm2) # classical AIC
#' res$AICc0-2*(1*(1+1))/(nrow(x)-1-1) # the same
#'
#' @export

assign("ortho.AIC", 
       function(Y, X, ord.var = FALSE) {
         # Fast Forward Selection AIC if X is orthonormal (XtX=I)
         # return a vector of AICc
         # if ord.var=TRUE, a list containing also order of variables is returned
         if (sum(apply(as.matrix(apply(X, 2, mean)), 1, function(x)
           identical(all.equal(x, 0), TRUE))) != ncol(X))
           stop("X variables are not centered")
         X <- sweep(X, 2, sqrt(colSums(X^2)), "/")
         if (!(sum(identical(all.equal(
           sum(crossprod(as.matrix(X)) - diag(ncol(X))), 0
         ), TRUE))))
           stop("X variables are not orthonormalized")
         
         f1 <- function(resp, X) {
           R2 <- t(as.matrix(X)) %*% as.matrix(Y)
           R2 <- t(R2) %*% R2
           return(sum(diag(R2)))
         }
         
         Y <- scale(Y, scale = FALSE)
         Y <- as.matrix(Y)
         R2 <- apply(X, 2, f1, resp = Y)
         SSTot <- sum(diag(t(Y) %*% Y))
         RSS <- SSTot - R2
         ordre <- order(RSS)
         RSS <- sort(RSS)
         R2 <- R2[ordre]
         RSScum <- cumsum(c(RSS[1], -R2[-1]))
         RSScum <- c(SSTot, RSScum)
         # By default, Y is centered
         # K is the number of othogonal vectors + 1 (for intercept)
         
         K <- (1 + (0:ncol(X)))
         AICtri <-
           nrow(X) * log(ifelse(RSScum <= 0, NA, RSScum / nrow(X))) + 2 * K
         correct <- 2 * (K * (K + 1)) / (nrow(X) - K - 1)
         correct <- ifelse(is.finite(correct) & (correct > 0), correct, NA)
         AICc <- AICtri + correct
         if (ord.var) {
           AICc <-
             list(
               AICc = AICc[-1],
               AICc0 = AICc[1],
               ord = ordre,
               R2 = cumsum(R2 / SSTot)
             )
         }
         if (!ord.var)
           AICc <- AICc[-1]
         return(AICc)
       }
)