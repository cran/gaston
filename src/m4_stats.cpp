// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "matrix4.h"
#include "loubar.h"

using namespace Rcpp;
using namespace RcppParallel;

uint8_t N0[256] = {
4, 3, 3, 3, 3, 2, 2, 2, 3, 2, 2, 2, 3, 2, 2, 2, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};

uint8_t N1[256] = {
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
2, 3, 2, 2, 3, 4, 3, 3, 2, 3, 2, 2, 2, 3, 2, 2, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0};

List geno_stats(matrix4 & A) {
  // Stats SNP
  IntegerMatrix SN(4,A.nrow);

  // Stats individus
  IntegerMatrix sIN(4, A.true_ncol*4);
  int * pIN = sIN.begin();

  for(size_t i = 0; i < A.nrow; i++) {
    for(size_t j = 0; j < A.true_ncol; j++) {
      uint8_t d = A.data[i][j];
      SN(0,i) += N0[d];
      SN(3,i) += N0[255-d];
      SN(1,i) += N1[d];
      SN(2,i) += N1[255-d];

      pIN[16*j + ((int) d&3)]++;
      pIN[16*j + 4  + ((int) (d>>2)&3)]++;
      pIN[16*j + 8  + ((int) (d>>4)&3)]++;
      pIN[16*j + 12 + ((int) (d>>6)&3)]++;
    }
    SN(3,i) -= (4*A.true_ncol - A.ncol);
  }

  sIN = sIN(_, Range(0,A.ncol-1));
  List L;

  L["snps"] = DataFrame::create(Named("N0") = SN(0,_), Named("N1")  = SN(1,_), 
                                Named("N2") = SN(2,_), Named("NAs") = SN(3,_) );
  L["inds"] = DataFrame::create(Named("N0") = sIN(0,_), Named("N1")  = sIN(1,_), 
                                Named("N2") = sIN(2,_), Named("NAs") = sIN(3,_) );
  return L;
}


DataFrame geno_stats_snp(matrix4 & A) {
  // Stats SNP
  IntegerMatrix SN(4,A.nrow);

  for(size_t i = 0; i < A.nrow; i++) {
    for(size_t j = 0; j < A.true_ncol; j++) {
      uint8_t d = A.data[i][j];
      SN(0,i) += N0[d];
      SN(3,i) += N0[255-d];
      SN(1,i) += N1[d];
      SN(2,i) += N1[255-d];
    }
    SN(3,i) -= (4*A.true_ncol - A.ncol);
  }

  return DataFrame::create(Named("N0") = SN(0,_), Named("N1")  = SN(1,_), 
                           Named("N2") = SN(2,_), Named("NAs") = SN(3,_) );
}

DataFrame geno_stats_ind(matrix4 & A) {
  // Stats individus
  IntegerMatrix sIN(4, A.true_ncol*4);
  int * pIN = sIN.begin();

  for(size_t i = 0; i < A.nrow; i++) {
    for(size_t j = 0; j < A.true_ncol; j++) {
      uint8_t d = A.data[i][j];
      pIN[16*j + ((int) d&3)]++;
      pIN[16*j + 4  + ((int) (d>>2)&3)]++;
      pIN[16*j + 8  + ((int) (d>>4)&3)]++;
      pIN[16*j + 12 + ((int) (d>>6)&3)]++;
    }
  }

  sIN = sIN(_, Range(0,A.ncol-1));
  return DataFrame::create(Named("N0") = sIN(0,_), Named("N1")  = sIN(1,_), 
                           Named("N2") = sIN(2,_), Named("NAs") = sIN(3,_) );
}

//[[Rcpp::export]]
List geno_stats(XPtr<matrix4> p_A) {
  return geno_stats(*p_A);
}

//[[Rcpp::export]]
DataFrame geno_stats_snp(XPtr<matrix4> p_A) {
  return geno_stats_snp(*p_A);
}

//[[Rcpp::export]]
DataFrame geno_stats_ind(XPtr<matrix4> p_A) {
  return geno_stats_ind(*p_A);
}

RcppExport SEXP gg_geno_stats(SEXP p_ASEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP );
        List __result = geno_stats(p_A);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

RcppExport SEXP gg_geno_stats_snp(SEXP p_ASEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP );
        DataFrame __result = geno_stats_snp(p_A);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

RcppExport SEXP gg_geno_stats_ind(SEXP p_ASEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP );
        DataFrame __result = geno_stats_ind(p_A);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}


