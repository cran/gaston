\name{set.stats}
\alias{set.stats}

\title{ Basic statistics for a \code{\link{bed.matrix}} }

\description{
  Return an updated \code{\link{bed.matrix}} with new variables for 
  several basic statistics.
}
\usage{
set.stats(x, set.ped_stats=TRUE, set.snps_stats=TRUE, set.p = TRUE, 
          set.mu_sigma = TRUE, verbose = getOption("gaston.verbose", TRUE))
}

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
  \item{set.ped_stats}{ If \code{TRUE}, \code{x@ped} is updated}
  \item{set.snps_stats}{ If \code{TRUE}, \code{x@snps} is updated}
  \item{set.p}{ If \code{TRUE}, \code{x@p} is updated }
  \item{set.mu_sigma}{ If \code{TRUE}, \code{x@mu} and \code{x@sigma} are updated }
 \item{verbose}{ If \code{TRUE}, display information on the function actions }
}

\details{
  If \code{set.ped_stats = TRUE}, the following variables are updated in \code{x@ped} :
\itemize{
  \item \code{N0}, \code{N1}, \code{N2} and \code{NAs} give for each individual the number of SNPs
  with a genotype equal to 0, 1, 2 and missing, respectively
  \item \code{callrate} is the individual callrate.
  }

  If \code{set.snps_stats = TRUE}, the following variables are updated in \code{x@snps} :
\itemize{
  \item \code{N0}, \code{N1}, \code{N2} and \code{NAs} give for each SNP the number of individuals
  with a genotype equal to 0, 1, 2 and missing, respectively
  \item \code{callrate} is the SNP callrate
  \item \code{maf} is the Minor Allele Frequency 
  \item \code{hz} is the heterozygosity 
  }

  If \code{set.p = TRUE}, \code{x@p} is updated with the alternate allele frequency.

  If \code{set.mu_sigma = TRUE}, \code{x@mu} is updated with the genotype mean (equal to \code{2*x@p})
  and \code{x@sigma} with the genotype standard deviation (should be approximately \code{sqrt(2*x@p*(1-x@p))}
  under Hardy-Weinberg Equilibrium).
}

\value{ 
  A \code{\link{bed.matrix}} similar to \code{x}, with slots updated as described above.
}

\author{ Hervé Perdry and Claire Dandine-Roulland }

\seealso{ \code{\link{set.hwe}} }

\examples{
# Load data
data(TTN)
x <- as.bed.matrix(TTN.gen, TTN.fam, TTN.bim)

# Compute statistics
x <- set.stats(x)
str(x@ped)
str(x@snps)
}

\keyword{ Genetic }
\keyword{ Statistics }
