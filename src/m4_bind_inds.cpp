#include <Rcpp.h>
#include <iostream>
#include "matrix4.h"

using namespace Rcpp;

// [[Rcpp::export]]
XPtr<matrix4> bind_inds(List L) {
  int s = L.size();
  if(s < 2) Rf_error("Can't bind less than two matrices!");
  XPtr<matrix4> first = as<XPtr<matrix4> >(L[0]);
  int n = first->ncol;
  int m = first->nrow;
  for(int i = 1; i < s; i++) {
    XPtr<matrix4> nxt = as<XPtr<matrix4> >(L[i]);
    if(m != nxt->nrow) Rf_error("Dimensions mismatch");
    n += nxt->ncol;
  }
  XPtr<matrix4> r(new matrix4(m,n));
  for(int i = 0; i < m; i++) {
    int k = 0;
    for(int j = 0; j < s; j++) {
      XPtr<matrix4> nxt = as<XPtr<matrix4> >(L[j]);
      for(int jj = 0; jj < nxt->ncol; jj++) {
        (*r)(i,k++) = (*nxt)(i,jj); }
    }
  }
  return r;
}

RcppExport SEXP gg_bind_inds(SEXP LSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type L(LSEXP );
        XPtr<matrix4> __result = bind_inds(L);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

