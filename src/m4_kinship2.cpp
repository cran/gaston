// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "matrix4.h"
#include "loubar.h"

using namespace Rcpp;
using namespace RcppParallel;


// ************* [on ne symmétrise pas] ***********

struct paraKin2 : public Worker {
  // input and others
  char * data;
  const size_t ncol;
  double * gg;
  // output
  double * K;
  
  // constructeur
  paraKin2(char * data, size_t ncol, double * gg, double * K) 
         : data(data), ncol(ncol), gg(gg), K(K) {}

  // worker !
  void operator()(size_t beg, size_t end) {
    for(size_t j1 = beg; j1 < end; j1++) {
      char x1 = data[j1];
      for(int ss1 = 0; ss1 < 4 && 4*j1 + ss1 < ncol; ss1++) {
        double * ggg = gg + ((x1&3)<<2);
        size_t k = (4*j1+ss1)*ncol;
        for(size_t j2 = 0; j2 < j1; j2++) {
          char x2 = data[j2];
          for(int ss2 = 0; ss2 < 4; ss2++) {
            K[k++] += ggg[x2&3];
            x2 >>= 2;
          }
        }
        size_t j2 = j1;
        char x2 = data[j2];
        for(int ss2 = 0; ss2 <= ss1; ss2++) {
          K[k++] += ggg[x2&3];
          x2 >>= 2;
        }
        x1 >>= 2;
      }
    }
  }
};




// [[Rcpp::export]]
NumericMatrix Kinship2(XPtr<matrix4> p_A, const std::vector<double> & mu, const std::vector<double> & sd, int chunk) {

  NumericMatrix Y(p_A->ncol,p_A->ncol);

  for(size_t i = 0; i < p_A->nrow; i++) {
    double gg[16];
    gg[3] = gg[7] = gg[11] = gg[12] = gg[13] = gg[14] = gg[15] = 0;
    double sd_ = (sd[i]==0)?1:sd[i];
    double mu_ = mu[i];
    double v0 = -mu_/sd_;
    double v1 = (1-mu_)/sd_;
    double v2 = (2-mu_)/sd_;

    gg[0] = v0*v0; 
    gg[1] = gg[4] = v1*v0;
    gg[2] = gg[8] = v2*v0;
    gg[5] = v1*v1;
    gg[6] = gg[9] = v2*v1;
    gg[10] = v2*v2;

    paraKin2 X(p_A->data[i], p_A->ncol, gg, Y.begin());
    parallelFor(0, p_A->true_ncol, X, chunk);
  }

  // symmetriser
  for(size_t i = 0; i < p_A->ncol; i++) {
    for(size_t j = 0; j < i; j++) {
      Y(i,j) = Y(j,i);
    }
  }

  return Y;
}

RcppExport SEXP gg_Kinship2(SEXP p_ASEXP, SEXP muSEXP, SEXP sdSEXP, SEXP chunkSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type mu(muSEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type sd(sdSEXP );
        Rcpp::traits::input_parameter< int >::type chunk(chunkSEXP );
        NumericMatrix __result = Kinship2(p_A, mu, sd, chunk);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}


