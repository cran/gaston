// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
using namespace Rcpp;
using namespace Eigen;

typedef Map<MatrixXd> Map_MatrixXd;

NumericMatrix crossprod(NumericMatrix X) {
  Map_MatrixXd xx(as<Map<MatrixXd> >(X));
  int n(xx.cols());
  MatrixXd aa( MatrixXd(n,n).setZero().selfadjointView<Lower>().rankUpdate(xx.adjoint()));
  return wrap(aa);
}

NumericMatrix tcrossprod(NumericMatrix X) {
  Map_MatrixXd xx(as<Map<MatrixXd> >(X));
  int n(xx.rows());
  MatrixXd aa( MatrixXd(n,n).setZero().selfadjointView<Lower>().rankUpdate(xx));
  return wrap(aa);
}

RcppExport SEXP gg_crossprod(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP );
        NumericMatrix __result = crossprod(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

RcppExport SEXP gg_tcrossprod(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP );
        NumericMatrix __result = tcrossprod(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

