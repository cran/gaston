
require(gaston)

if(!require("HGDP.CEPH")) install.packages("HGDP.CEPH", repos="https://genostats.github.io/R/")
filepath <-system.file("extdata", "hgdp_ceph.bed", package="HGDP.CEPH")
x <- read.bed.matrix(filepath)
x <- set.stats(x)

# two bed matrices with the same genetic map
I <- sample.int(nrow(x), 500)
x1 <- x[ I, ]
x2 <- x[ -I, ]

# PCA of the first matrix
K <- GRM(x1)
eiK <- eigen(K)
plot( sqrt(eiK$values[1]) * eiK$vectors[,1], sqrt(eiK$values[2]) * eiK$vectors[,2] , col = x1@ped$region7)

# loadings of the first two PCs
L <- bed.loadings(x1, eiK$vectors[,1:2])

# recompute the first two PCs (rounding errors may lead to slightly different values
standardize(x1) <- "p"
PC1 <- (x1 %*% L) / sqrt(nrow(L)-1)
plot(PC1, col = x1@ped$region7)

# projection of x2 usinh the loadings
x2@p <- x1@p # use same allelic frequences should slightly improve the result
standardize(x2) <- "p"
PC2 <- (x2 %*% L) / sqrt(nrow(L)-1)
points(PC2, pch = 5, col = x2@ped$region7)

