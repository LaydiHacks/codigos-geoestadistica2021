

W <- as.matrix(as_dgRMatrix_listw(a.lw))   # wc$n=49
S0 <- t(rep(1,49))%*%(W)%*%rep(1,49)
S02 <- S0^2
S3 <- t(rep(1,49))%*%(W*t(W))%*%rep(1,49)
S4 <- t(rep(1,49))%*%(W*W)%*%rep(1,49)
S5 <- t(rep(1,49))%*%(W%*%W)%*%rep(1,49)
S6 <- t(rep(1,49))%*%(t(W)%*%W+W%*%t(W))%*%rep(1,49)
S1 <- S3+S4 # t(rep(1,49))%*%(W*W+W*t(W))%*%rep(1,49)
S2 <- 2*S5+S6
rxy <- cor(x,y)
K <- sum((as.numeric(scale(x,center=TRUE, scale=F))^2)*as.numeric(scale(y,center=TRUE, scale=F))^2)/(((var(x)*wc$n1)/wc$n)*((var(y)*wc$n1)/wc$n))
EI <- -rxy/(wc$n1)
EI2 <- ((rxy^2)*wc$n*(2*(S02-S2+S1)+(2*S3-2*S5)*wc$n3+S3*wc$n2*wc$n3)-K*(6*(S02-S2+S1)+(4*S1-2*S2)*wc$n3+S1*wc$n2*wc$n3)+
       wc$n*((S02-S2+S1)+(2*S4-S6)*wc$n3+S4*wc$n2*wc$n3))/(wc$n1*wc$n2*wc$n3*S02)
VI <- EI2-EI^2
Z(Ixy) <- (I-EI)/sqrt(VI)
