geary.bi <- function(x, y, listw, zero.policy=NULL,adjust.n=TRUE,alternative="greater") {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
        stopifnot(is.vector(x))
        stopifnot(is.vector(y))
  wc <- spweights.constants(listw, zero.policy = zero.policy, adjust.n = adjust.n)
  n <- wc$n
  n1 <- wc$n1
  S0 <- wc$S0
  n2 <- length(x)
  res <- as.matrix(as_dgRMatrix_listw(listw))*(kronecker(matrix(x,nrow=n2),t(rep(1,n2)))-
         t(kronecker(matrix(y,nrow=n2),t(rep(1,n2)))))^2
	z <- scale(x, scale=FALSE)
	zz <- sum(z^2)
	K <- (n*sum(z^4))/(zz^2)
	C <- (n1 / (2*S0)) * (sum(res) / zz)
	res <- list(C=C, K=K)
	res
}

gearybi.test <- function(x, y, listw, randomisation=TRUE, zero.policy=NULL,
    alternative="greater", spChk=NULL, adjust.n=TRUE) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	alternative <- match.arg(alternative, c("less", "greater", "two.sided"))
	if(!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.numeric(x)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
	if (any(is.na(x))) stop("NA in X")
	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw))
		stop("Check of data and weights ID integrity failed")
	
	wc <- spweights.constants(listw, zero.policy, adjust.n=adjust.n)
	S02 <- wc$S0*wc$S0
	res <- geary.bi(x, y, listw, zero.policy,adjust.n)
	C <- res$C
	if (is.na(C)) stop("NAs generated in geary - check zero.policy")
	K <- res$K
	EC <- 1
	if(randomisation) {
		VC <- (wc$n1*wc$S1*(wc$nn - 3*n + 3 - K*wc$n1))
		VC <- VC - ((1/4) * (wc$n1*wc$S2*(wc$nn + 3*n - 6 - 
			K*(wc$nn - n + 2))))
		VC <- VC + (S02*(wc$nn - 3 - K*(wc$n1^2)))
		VC <- VC / (n*wc$n2*wc$n3*S02)
	} else {
		VC <- ((2*wc$S1 + wc$S2)*wc$n1 - 4*S02) / (2*(n + 1)*S02)
	}
#	ZC <- (C - EC) / sqrt(VC)
# order changed 090609 RSB (C&O 1973, p. 21)
	ZC <- (EC - C) / sqrt(VC)
	statistic <- ZC
	names(statistic) <- "Geary C statistic ZC"
	PrC <- NA
	if (is.finite(ZC)) {
        	if (alternative == "two.sided") PrC <- 2 * pnorm(abs(ZC), 
			lower.tail=FALSE)
        	else if (alternative == "greater")
            	PrC <- pnorm(ZC, lower.tail=FALSE)
        	else PrC <- pnorm(ZC)
		if (!is.finite(PrC) || PrC < 0 || PrC > 1) 
		    warning("Out-of-range p-value: reconsider test arguments")
	}
	vec <- c(C, EC, VC)
	names(vec) <- c("Geary C statistic", "Expectation", "Variance")
	method <- paste("Geary C test under", ifelse(randomisation,
	    "randomisation", "normality"))
	data.name <- paste(deparse(substitute(x)), "\nweights:",
	    deparse(substitute(listw)), "\n")
	res <- list(statistic=statistic, p.value=PrC, estimate=vec, 
	    alternative=ifelse(alternative == "two.sided", alternative, 
	    paste("Expectation", alternative, "than statistic")), 
	    method=method, data.name=data.name)
	class(res) <- "htest"
	res
}

                                                   
