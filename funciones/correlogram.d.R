assign("correlogram.d",
function (coords, z, method = "Moran", nbclass = NULL, zero.policy=NULL, ...) 
{
    coords <- as.matrix(coords)
    matdist <- dist(coords)
    if (is.null(nbclass)) 
        nbclass <- nclass.Sturges(matdist)
    etendue <- range(matdist)
    breaks1 <- seq(etendue[1], etendue[2], l = nbclass + 1)
    breaks2 <- breaks1 + 1e-06
    breaks <- cbind(breaks1[1:length(breaks1) - 1], breaks2[2:length(breaks2)])
    breaks[1, 1] <- breaks[1, 1] - 1e-06
    lst.nb1 <- rep(list(NA), nbclass)
    lst.z1 <- rep(list(NA), nbclass)
    for (i in 1:length(breaks[, 1])) {
        lst.z1[[i]] <- z
        lst.nb1[[i]] <- dnearneigh(coords, breaks[i, 1], breaks[i, 
            2])
        zero <- which(card(lst.nb1[[i]]) == 0)
        if (length(zero) > 0) {
            lst.nb1[[i]] <- dnearneigh(coords[-zero, ], breaks[i, 
                1], breaks[i, 2])
            lst.z1[[i]] <- z[-zero]
        }
    }
    lst.res1 <- rep(list(NA), nbclass)
    for (i in 1:length(breaks[, 1])) {
        xt <- switch(pmatch(method, c("Moran", "Geary"), nomatch = 3), 
            try(moran.test(lst.z1[[i]], nb2listw(lst.nb1[[i]], 
                style = "W"),zero.policy=zero.policy, ...), silent = TRUE), try(geary.test(lst.z1[[i]], 
                nb2listw(lst.nb1[[i]], style = "W"),zero.policy=zero.policy, ...), silent = TRUE), 
            stop("Method must be 'Moran' or 'Geary'"))
        if (inherits(xt, "try-error")) {
            stop("Bad selection of class breaks, try another one...")
        }
else {
            x <- xt$estimate[1]
            p <- xt$p.value
            N <- sum(card(lst.nb1[[i]]))
        }
        lst.res1[[i]] <- c(x = x, p = p, N = N)
    }
    meth <- names(xt[[3]][1])
    mat <- matrix(unlist(lst.res1), ncol = 3, byrow = TRUE)
    res <- cbind(dist.class = rowMeans(breaks), coef = mat[, 
        1], p.value = mat[, 2], n = mat[, 3], low.l = breaks[,1], up.l = breaks[,2])
    attributes(res) <- c(attributes(res), list(Method = meth))
    class(res) <- c("correlog", "matrix")
    res1 <- list(res=res,z1=lst.z1, dnn=lst.nb1)
    res1
}
)