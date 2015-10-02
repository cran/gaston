\name{association.test}
\alias{association.test}

\title{ Association Test }

\description{ Association tests between phenotype and SNPs. }
\usage{ 
association.test(x, Y = x@ped$pheno, X = matrix(1, nrow(x)), 
    eigenK, beg = 1, end = ncol(x), p = 0, test = c("wald", "lrt"),
    tol = .Machine$double.eps^0.25, multithreaded = FALSE)
}

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{Y}{ The phenotype vector. Default is the column (\code{pheno}) of \code{x@ped} }
  \item{X}{ A covariable matrix. The default is a column vector of ones, to include an intercept in the model }
  \item{eigenK}{ Eigen decomposition of the Genetic Relationship Matrix (as produced by the function \code{eigen}) }
  \item{beg}{ Index of the first SNP tested for association }
  \item{end}{ Index of the last SNP tested for association }
  \item{p}{ Number of Principal Components to include in the mixed model with fixed effect }
  \item{test}{ Which association test to use: \code{"wald"} for the Wald test or \code{"lrt"} for the likelihood ratio test }
  \item{tol}{ Parameter for the likelihood maximization (as in \code{optimize}) }
  \item{multithreaded}{ \code{Logical}. If \code{TRUE}, several threads are used. }
}

\details{
  Tests the association between the phenotype and requested SNPs in \code{x}.
  
  For each SNP, the following model in considered
  \deqn{ Y = (X|PC)\alpha + G\beta + \omega + \varepsilon }{ Y = (X|PC) alpha + G beta + omega + epsilon }
  with \eqn{ \omega \sim N(0,\tau K) }{omega ~ N(0, tau K)} and \eqn{ \varepsilon \sim N(0,\sigma^2 I_n) }{epsilon ~ N(0, sigma^2 I_n)}.
  \eqn{G} is the genotype vector of the SNP, \eqn{K} is the Genetic Relationship Matrix (GRM) computed on the whole genome 
  (specified via the parameter \code{eigenK}), 
  \eqn{X} the covariable matrix, and \eqn{PC} the matrix of the first \eqn{p} principal components.
  All parameters are estimated with the same procedure as in \code{lmm.diago}.

  Note: in the present implementation, using several threads is likely to be slower than
  using only one, hence the \code{multithreaded} parameter. This might change in the near future.
}

\value{
  A named list of vectors, giving for each considered SNP:
  \item{h2}{ Estimated value of \eqn{\tau \over {\tau + \sigma^2}}{tau/(tau + sigma^2)} }
  \item{beta}{ (only whith \code{test = "wald"}) Estimated value of \eqn{\beta}{beta} }
  \item{sd}{ (only whith \code{test = "wald"}) Estimated standard deviation of the \eqn{\beta}{beta} estimation }
  \item{LRT}{ (only whith \code{test = "lrt"}) Value of the Likelihood Ratio Test}
  \item{p}{ The corresponding p-value }
}

\seealso{ \code{\link{lmm.diago}}, \code{\link[stats:optimize]{optimize}} }

\author{ Hervé Perdry and Claire Dandine-Roulland }

\examples{
\donttest{
# Load data
data(AGT)
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
x <- set.stats(x)
standardize(x) <- 'mu'

# generate a random positive matrix (to play the role of the GRM)
set.seed(1)
R <- random.pm(nrow(x))

# simulate phenotype with effect of SNP #351 and a polygenic component 
y <- x \%*\% c(rep(0,350),0.25,rep(0,ncol(x)-351)) + lmm.simu(0.3,1,eigenK=R$eigen)$y

# association test
t <- association.test(x, y, eigenK = R$eigen)

# (mini) Manhattan plot
plot(-log10(t$p), xlab="SNP index", ylab = "-log(p)")

# link between p-values and LD with SNP #351
lds <- LD(x, 351, c(1,ncol(x)))
plot(lds, -log10(t$p), xlab="r^2", ylab="-log(p)")
}}

\keyword{ Association Test }
