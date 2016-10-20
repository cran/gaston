## Basic Statistics
set.stats <- function(x, set.p = TRUE, set.mu_sigma = TRUE, verbose = getOption("gaston.verbose",TRUE)) {
  if( is(x)!='bed.matrix' ) stop('x must be an object of class bed.matrix')
  if(!is.logical(set.p) | !is.logical(set.mu_sigma)) 
    stop('set.* arguments must be logical')

  w.a <- x@snps$chr %in% getOption("gaston.autosomes")
  w.x <- x@snps$chr %in% getOption("gaston.chr.x")
  w.y  <- x@snps$chr %in% getOption("gaston.chr.y")
  w.mt <- x@snps$chr %in% getOption("gaston.chr.mt")
  w.f  <- x@ped$sex  == 2
  st <- .Call('gg_geno_stats', PACKAGE = 'gaston', x@bed, w.x, w.y, w.mt, w.f) 

  ############  completer snps
  st$snps$callrate <- 1-st$snps$NAs/nrow(x)
  # correction pour chr y 
  st$snps$callrate[w.y] <- 1-(st$snps$NAs[w.y]-st$snps$NAs.f[w.y])/sum(x@ped$sex == 1)

  # freq allele alt
  n <- nrow(x) - st$snps$NAs;
  pp <- (2*st$snps$N2 + st$snps$N1)/(2*n);
  # correction pour chr x
  a <- st$snps$N2.f[w.x] + st$snps$N1.f[w.x] + st$snps$N2[w.x]
  b <- st$snps$N0.f[w.x] + st$snps$N1.f[w.x] + st$snps$N0[w.x]
  pp[w.x] <- a/(a+b)

  st$snps$maf <- pmin(pp,1-pp)
  st$snps$hz <- st$snps$N1/n

  # correction pour chr x
  nb.fe <- sum(x@ped$sex == 2)
  st$snps$hz[w.x] <- st$snps$N1.f[w.x]/(nb.fe - st$snps$NAs.f[w.x])

  ############ completer inds/ped
  n.a <- sum(w.a)
  st$inds$callrate <- 1-st$inds$NAs/n.a
  st$inds$hz <- st$inds$N1/(n.a-st$inds$NAs)

  n.x <- sum(w.x)
  st$inds$callrate.x <- 1-st$inds$NAs.x/n.x
  st$inds$hz.x <- st$inds$N1.x/(n.x-st$inds$NAs.x)

  n.y <- sum(w.y)
  st$inds$callrate.y <- 1-st$inds$NAs.y/n.y
  st$inds$hz.y <- st$inds$N1.mt/(n.y-st$inds$NAs.y)

  n.mt <- sum(w.mt)
  st$inds$callrate.mt <- 1-st$inds$NAs.mt/n.mt
  st$inds$hz.mt <- st$inds$N1.mt/(n.mt-st$inds$NAs.mt)

  ########## insérer dans x
  x@snps[, names(st$snps)] <- st$snps
  x@ped[ , names(st$inds)] <- st$inds

  if(verbose) cat("ped stats and snps stats have been set. \n") 

  if(set.p) {
    x@p <- pp
    if(verbose) cat("'p' has been set. \n")
  }  

  if(set.mu_sigma) { # calcul brutal
    n <- nrow(x) - x@snps$NAs;
    mu <- (2*x@snps$N2 + x@snps$N1)/n;
    N <- nrow(x)
    s <- sqrt( (x@snps$N1 + 4*x@snps$N2 + mu**2*x@snps$NAs)/(N-1) - N/(N-1)*mu**2 )
    x@mu <- mu;
    x@sigma <- s
    if(verbose) cat("'mu' and 'sigma' have been set.\n");
  }
  x
}

## Basic Statistics, only individuals ---------------------------------------------------------
set.stats.ped <- function(x, verbose = getOption("gaston.verbose",TRUE)) {
  if( is(x)!='bed.matrix' ) stop('x must be an object of class bed.matrix')

  w.a <- x@snps$chr %in% getOption("gaston.autosomes")
  w.x <- x@snps$chr %in% getOption("gaston.chr.x")
  w.y  <- x@snps$chr %in% getOption("gaston.chr.y")
  w.mt <- x@snps$chr %in% getOption("gaston.chr.mt")
  st <- .Call('gg_geno_stats_inds', PACKAGE = 'gaston', x@bed, w.x, w.y, w.mt) 

  ############ completer inds/ped
  n.a <- sum(w.a)
  st$inds$callrate <- 1-st$inds$NAs/n.a
  st$inds$hz <- st$inds$N1/(n.a-st$inds$NAs)

  n.x <- sum(w.x)
  st$inds$callrate.x <- 1-st$inds$NAs.x/n.x
  st$inds$hz.x <- st$inds$N1.x/(n.x-st$inds$NAs.x)

  n.y <- sum(w.y)
  st$inds$callrate.y <- 1-st$inds$NAs.y/n.y
  st$inds$hz.y <- st$inds$N1.mt/(n.y-st$inds$NAs.y)

  n.mt <- sum(w.mt)
  st$inds$callrate.mt <- 1-st$inds$NAs.mt/n.mt
  st$inds$hz.mt <- st$inds$N1.mt/(n.mt-st$inds$NAs.mt)

  ########## insérer dans x
  x@ped[ , names(st$inds)] <- st$inds

  if(verbose) cat("ped stats have been set. \n") 
  x
}

## Basic Statistics, only snps -----------------------------------------------------------------------------------
set.stats.snps <- function(x, set.p = TRUE, set.mu_sigma = TRUE, verbose = getOption("gaston.verbose",TRUE)) {
  if( is(x)!='bed.matrix' ) stop('x must be an object of class bed.matrix')
  if(!is.logical(set.p) | !is.logical(set.mu_sigma)) 
    stop('set.* arguments must be logical')

  w.x <- x@snps$chr %in% getOption("gaston.chr.x")
  w.y <- x@snps$chr %in% getOption("gaston.chr.y")
  w.f <- x@ped$sex  == 2
  st <- .Call('gg_geno_stats_snps', PACKAGE = 'gaston', x@bed, w.x, w.f) 

  ############  completer snps
  st$snps$callrate <- 1-st$snps$NAs/nrow(x)
  # correction pour chr y 
  st$snps$callrate[w.y] <- 1-(st$snps$NAs[w.y]-st$snps$NAs.f[w.y])/sum(x@ped$sex == 1)

  # freq allele alt
  n <- nrow(x) - st$snps$NAs;
  pp <- (2*st$snps$N2 + st$snps$N1)/(2*n);
  # correction pour chr x
  a <- st$snps$N2.f[w.x] + st$snps$N1.f[w.x] + st$snps$N2[w.x]
  b <- st$snps$N0.f[w.x] + st$snps$N1.f[w.x] + st$snps$N0[w.x]
  pp[w.x] <- a/(a+b)

  st$snps$maf <- pmin(pp,1-pp)
  st$snps$hz <- st$snps$N1/n

  # correction pour chr x
  nb.fe <- sum(x@ped$sex == 2)
  st$snps$hz[w.x] <- st$snps$N1.f[w.x]/(nb.fe - st$snps$NAs.f[w.x])

  ########## insérer dans x
  x@snps[, names(st$snps)] <- st$snps

  if(verbose) cat("snps stats have been set. \n") 

  if(set.p) {
    x@p <- pp
    if(verbose) cat("'p' has been set. \n")
  }  

  if(set.mu_sigma) { # calcul brutal
    n <- nrow(x) - x@snps$NAs;
    mu <- (2*x@snps$N2 + x@snps$N1)/n;
    N <- nrow(x)
    s <- sqrt( (x@snps$N1 + 4*x@snps$N2 + mu**2*x@snps$NAs)/(N-1) - N/(N-1)*mu**2 )
    x@mu <- mu;
    x@sigma <- s
    if(verbose) cat("'mu' and 'sigma' have been set.\n");
  }
  x
}




set.hwe <- function(x, method = c("chisquare", "exact"), verbose = getOption("gaston.verbose",TRUE)) {
  if( !all(c("N0", "N1", "N2") %in% names(x@snps) )) {
    if(verbose) cat("Computing basic stats\n")
    x <- set.stats(x)
  }
  method <- match.arg(method)
  if(method == 'chisquare') {
    if(verbose) cat("Computing HW chi-square p-values\n")
    hwe_ <- .Call('gg_hwe_chi', PACKAGE = 'gaston', x@snps$N0, x@snps$N1, x@snps$N2)
  } else {
    if(verbose) cat("Computing HW exact test p-values\n")
    hwe_ <- .Call('gg_hwe', PACKAGE = 'gaston', x@snps$N0, x@snps$N1, x@snps$N2)
  }
  x@snps$hwe <- hwe_
  x
}



