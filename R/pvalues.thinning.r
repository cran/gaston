pvalues.thinning <- function(p, max.bin.size = 100, bins = 100) {
  p <- p[ !is.nan(p) & p > 0 ]
  lp <- -log10(p) 
  s <- seq( min(lp), max(lp), length=bins+1)
  r <- NULL
  w1 <- TRUE
  for(i in 1:bins) {
    w2 <- (lp <= s[i+1])
    w <- which(w1 & w2)
    if(length(w) <= max.bin.size)
      r <- c(r, w)
    else
      r <- c(r, sample(w, max.bin.size))
    w1 <- !w2
  }
  r
}


