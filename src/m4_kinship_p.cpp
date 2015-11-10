// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "matrix4.h"
#include "loubar.h"
#include "m4_kinship_type.h"

using namespace Rcpp;
using namespace RcppParallel;


// ************* [on ne symmétrise pas] ***********

struct paraKin_p : public Worker {
  // input and others
  const matrix4 & A;
  const std::vector<double> p;
  const size_t ncol;
  const size_t true_ncol;
  const Ktype n;
  // output
  Ktype * K __attribute__ ((aligned(16)));
  
  // constructeurs
  paraKin_p(matrix4 & A, std::vector<double> p) : A(A), p(p), n(A.nrow-1),  // ******* BEWARE THE nrow-1 !!! ***********
        ncol(A.ncol), true_ncol(A.true_ncol) { 
          K = new Ktype[(ncol*(ncol+1))/2];
          std::fill(K, K+(ncol*(ncol+1))/2, 0);
        }
  paraKin_p(paraKin_p & Q, Split) : A(Q.A), p(Q.p), n(Q.n), ncol(Q.ncol), true_ncol(Q.true_ncol) {
          K = new Ktype[(ncol*(ncol+1))/2];
          std::fill(K, K+(ncol*(ncol+1))/2, 0);
        }

  // destructeur
  ~paraKin_p() { 
          delete [] K; 
  }

  // worker !
  void operator()(size_t beg, size_t end) {
    Ktype gg[16] __attribute__ ((aligned(16)));
    gg[3] = gg[7] = gg[11] = gg[12] = gg[13] = gg[14] = gg[15] = 0;
    for(size_t i = beg; i < end; i++) {
      Ktype p_ = (Ktype) p[i]; 
      if(p_ == 0 || p_ == 1) continue;
      Ktype w_ = 1/sqrt(2*p_*(1-p_)*n);
      Ktype v0 = -2*p_*w_;
      Ktype v1 = (1-2*p_)*w_;
      Ktype v2 = (2-2*p_)*w_;
      size_t k = 0;

      gg[0] = v0*v0; 
      gg[1] = gg[4] = v1*v0;
      gg[2] = gg[8] = v2*v0;
      gg[5] = v1*v1;
      gg[6] = gg[9] = v2*v1;
      gg[10] = v2*v2;
      char * dd = A.data[i];
      for(size_t j1 = 0; j1 < true_ncol; j1++) {
        char x1 = dd[j1];
        for(unsigned int ss1 = 0; (ss1 < 4) && (4*j1 + ss1 < ncol); ss1++) {
          Ktype * ggg __attribute__ ((aligned(16))) = gg + ((x1&3)<<2);
          for(size_t j2 = 0; j2 < j1; j2++) {
            char x2 = dd[j2];
            for(unsigned int ss2 = 0; ss2 < 4; ss2++) {
              K[k++] += ggg[x2&3];
              x2 >>= 2;
            }
          }
          size_t j2 = j1;
          char x2 = dd[j2];
          for(unsigned int ss2 = 0; ss2 <= ss1; ss2++) {
            K[k++] += ggg[x2&3];
            x2 >>= 2;
          }
          x1 >>= 2;
        }
      }  
    }
  }

  // recoller
  void join(const paraKin_p & Q) {
    std::transform(K, K + (ncol*(ncol+1))/2, Q.K, K, std::plus<Ktype>());
    // autrement dit : K += Q.K;
  }

};




// [[Rcpp::export]]
NumericMatrix Kinship_p(XPtr<matrix4> p_A, const std::vector<double> & p, int chunk) {
  paraKin_p X(*p_A, p);
  parallelReduce(0, p_A->nrow, X, chunk);

  NumericMatrix Y(p_A->ncol,p_A->ncol);
  size_t k = 0;
  for(size_t i = 0; i < p_A->ncol; i++) {
    for(size_t j = 0; j <= i; j++) {
      Y(j,i) = (double) X.K[k++];
    }
  }

  // symmetriser
  k = 0;
  for(size_t i = 0; i < p_A->ncol; i++) {
    for(size_t j = 0; j <= i; j++) {
      Y(i,j) = (double) X.K[k++]; // ou Y(j,i)
    }
  }
  
  return Y;
}


RcppExport SEXP gg_Kinship_p(SEXP p_ASEXP, SEXP pSEXP, SEXP chunkSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type p(pSEXP );
        Rcpp::traits::input_parameter< int >::type chunk(chunkSEXP );
        NumericMatrix __result = Kinship_p(p_A, p, chunk);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

