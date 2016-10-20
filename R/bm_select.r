select.snps <- function(x, condition) {
  if(!is(x, "bed.matrix")) stop("x is not a bed.matrix")
  w <- eval(substitute(condition), x@snps, parent.frame())
  miss <- is.na(w)
  if(sum(miss)>0) {
    warning(paste(sum(miss), 'SNP(s) with undefined condition are removed from bed.matrix'))
    w <- w & !miss
  }
  x[,w]
}

select.inds <- function(x, condition) {
  if(!is(x, "bed.matrix")) stop("x is not a bed.matrix")
  w <- eval(substitute(condition), x@ped, parent.frame())
  miss <- is.na(w)
  if(sum(miss)>0) {
    warning(paste(sum(miss), 'individual(s) with undefined condition are removed from bed.matrix'))
    w <- w & !miss
  }
  x[w,]
}
