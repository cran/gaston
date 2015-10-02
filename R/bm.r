setClassUnion("data.frameOrNULL",members=c("data.frame", "NULL"))
setClassUnion("numericOrNULL",members=c("numeric", "NULL"))


pednames <- c("famid", "id", "father", "mother", "sex", "pheno")
snpnames <- c("chr", "id", "dist", "pos", "A1", "A2")

snpstatnames <- c("N0", "N1", "N2", "NAs", "callrate", "maf", "hz", "hwe")
pedstatnames <- c("N0", "N1", "N2", "NAs", "callrate")

is.null.df <- function(x) is.data.frame(x) & nrow(x) == 0 & ncol(x) == 0

## Class bed.matrix
setClass("bed.matrix", representation( 
                   ped = 'data.frame',
                   snps = 'data.frame', 
                   bed = 'externalptr',
                   p = 'numericOrNULL',
                   mu = 'numericOrNULL',
                   sigma = 'numericOrNULL',
                   standardize_p = 'logical',
                   standardize_mu_sigma = 'logical' ))

setValidity('bed.matrix',
           function(object) {
             errors <- character()
             if ( object@standardize_p & object@standardize_mu_sigma ) 
                errors <- c(errors, "Only one center scale parameter can be TRUE.")
             if ( object@standardize_p & is.null(object@p) ) 
                errors <- c(errors, "If 'standardize_p' is TRUE, 'p' must be defined.")
             if ( object@standardize_mu_sigma & ( is.null(object@mu) | is.null(object@sigma) ) ) 
                errors <- c(errors, "If 'standardize_mu_sigma' is TRUE, 'mu' and 'sigma' must be defined.")
             if ( !is.null(object@p) & length(object@p) != ncol(object) ) 
                errors <- c(errors, "The length of 'p' must be equal to the number of markers.")
             if ( !is.null(object@mu) & length(object@mu) != ncol(object) ) 
                errors <- c(errors, "The length of 'mu' must be equal to the number of markers.")
             if ( !is.null(object@sigma) & length(object@sigma) != ncol(object) ) 
                errors <- c(errors, "The length of 'sigma' must be equal to the number of markers.")
             if ( length(errors)==0 ) return(TRUE) else return(errors)
           } );


setAs("bed.matrix", "matrix",
  function(from) {
    validObject(from)
    to <- if(from@standardize_p) 
      .Call('gg_m4_as_scaled_matrix_p', PACKAGE = 'gaston', from@bed, from@p)
    else if(from@standardize_mu_sigma)
      .Call('gg_m4_as_scaled_matrix_mu_sigma', PACKAGE = 'gaston', from@bed, from@mu, from@sigma)
    else
      .Call('gg_m4_as012', PACKAGE = 'gaston', from@bed)
    colnames(to) <- from@snps$id
    rownames(to) <- if(any(duplicated(from@ped$id))) paste(from@ped$fam, from@ped$id, sep="_")
                    else from@ped$id
    to
  } );

setGeneric('as.matrix')
setMethod("as.matrix", signature="bed.matrix", definition = function(x) as(x,"matrix") )

setAs("matrix", "bed.matrix", 
  function(from) {
    bed <- .Call('gg_as_matrix4', PACKAGE = 'gaston', from)

    ped <- if(is.null(rownames(from))) 
             structure(list(), row.names = c(NA, -nrow(from)), class = "data.frame") # empty data frame with right number of lines
           else 
             data.frame(famid = rownames(from), id = rownames(from), father = 0, mother = 0, sex = 0, pheno = 0, stringsAsFactors = FALSE)

    snp <- if(is.null(colnames(from)))
             structure(list(), row.names = c(NA, -ncol(from)), class = "data.frame") #idem
           else 
             data.frame(chr = NA, id = colnames(from), dist = NA, pos = NA, A1 = NA, A2 = NA, stringsAsFactors = FALSE)

    new("bed.matrix", bed = bed, snps = snp, ped = ped, p = NULL, mu = NULL,
        sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )
  } );

setGeneric('as.bed.matrix',function(x, fam = NULL, bim = NULL) standardGeneric("as.bed.matrix"),package='gaston',valueClass="bed.matrix")
setMethod('as.bed.matrix', signature=c(x='matrix', fam = 'data.frameOrNULL', bim = 'data.frameOrNULL'),
          def=function(x, fam = NULL, bim = NULL) {
   if ( !is.null(fam) & !is.data.frame(fam) ) stop('fam must be a data.frame or NULL.')
   if ( !is.null(bim) & !is.data.frame(bim) ) stop('bim must be a data.frame or NULL.')
   to <- as(x,"bed.matrix")
   if (!is.null(fam))
   {
     to@ped <- if(all(pednames %in% names(fam))) 
                 data.frame(famid = fam$famid, id = fam$id, father = fam$father, mother = fam$mother, sex = fam$sex, pheno = fam$pheno, stringsAsFactors = FALSE)
               else 
                 stop('"fam" data frame must contains "famid", "id", "father", "mother", "sex" and "pheno" variables')
   }
   if (!is.null(bim)) {
     to@snps <- if(all(snpnames %in% names(bim))) 
                  data.frame(chr = bim$chr, id = bim$id, dist = bim$dist, pos = bim$pos, A1 = bim$A1, A2 = bim$A2, stringsAsFactors = FALSE)
               else 
                  stop('"bim" data frame must contains "chr", "id", "dist", "pos", "A1" and "A2" variables')
   }
   to } )

setGeneric('dim')
setMethod("dim", signature = "bed.matrix", 
    function(x) c(.Call('gg_ninds', PACKAGE = 'gaston', x@bed), .Call('gg_nsnps', PACKAGE = 'gaston', x@bed)))

setGeneric('head')
setMethod( 'head', signature(x='bed.matrix'), function(x, nrow=10, ncol=10) print( as.matrix( x[1:min( nrow, nrow(x) ),1:min( ncol, ncol(x) )] ) ) )

setMethod(show, signature("bed.matrix"), 
      function(object) { 
        cat('A bed.matrix with ', nrow(object), ' individuals and ', ncol(object), ' markers.\n', sep='')
      } )

	  
	  
