// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <iostream>
#include "matrix4.h"
#include "loubar.h"
#include "m4_ld.h"

#define debug false

using namespace Rcpp;

// un matrice de génotypes, les deux vecteurs nécessaires à centrer/réduire
// threshold : le seuil de r² accepté
// pos : les positions des SNPs 
// chr : les chromosomes 
// max_dist : la distance au dessus de laquelle on ne calcule pas le LD
// beg, end : les indices des SNPs à considérer (inclusifs)
// [[Rcpp::export]]
LogicalVector ld_thin_right(XPtr<matrix4> pA, NumericVector mu, NumericVector sd, double threshold, 
                      IntegerVector pos, IntegerVector chr, int max_dist, int beg, int end) {
  int n = end - beg + 1;
  LogicalVector w(n);
  for(int i = 0; i < n; i++) w(i) = true;
  int i = beg;
  threshold = sqrt(threshold);
  if(debug) std::cout << "threshold = " << threshold;
  while(i <= end) {
    if(debug) std::cout << "\ni = " << i << "\n";
    int j = i + 1;
    int max_pos = pos(i) + max_dist;
    int chr_i = chr(i);
    double mu_i = mu(i);
    double sd_i = sd(i);
    int next_i = 0;
    bool gotnexti = false;
    while(j <= end && pos(j) < max_pos && chr(j) == chr_i) {
      if(debug) std::cout << "\ni = " << i << " j = " << j;
      double ld = LD_colxx(*pA, mu_i, mu(j), sd_i*sd(j), i, j);
      if(debug) std::cout << " abs(ld) = " << fabs(ld);
      if(!gotnexti) {
        next_i = j;
        gotnexti = true;
      }
      if(fabs(ld) > threshold) {
        if(debug) std::cout << " removing snp " << i;
        w(i-beg) = false; 
        break;
      }
      j++;
    }
    if(debug) std::cout << " out of j loop with j = " << j ;
    if(gotnexti) i = next_i; else i = j; // seulement si la boucle j était vide ! (gap)
  }
  return w; 
}


RcppExport SEXP gg_ld_thin_right(SEXP pASEXP, SEXP muSEXP, SEXP sdSEXP, SEXP thresholdSEXP, SEXP posSEXP, SEXP chrSEXP, SEXP max_distSEXP, SEXP begSEXP, SEXP endSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP );
        Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type pos(posSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type chr(chrSEXP );
        Rcpp::traits::input_parameter< int >::type max_dist(max_distSEXP );
        Rcpp::traits::input_parameter< int >::type beg(begSEXP );
        Rcpp::traits::input_parameter< int >::type end(endSEXP );
        LogicalVector __result = ld_thin_right(pA, mu, sd, threshold, pos, chr, max_dist, beg, end);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

