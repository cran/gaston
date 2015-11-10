#include <RcppEigen.h>
#include <iostream>

using namespace Rcpp;
using namespace Eigen;

typedef Map<MatrixXd> Map_MatrixXd;

template<typename T1, typename T2, typename T3>
void AIREML1(const Eigen::MatrixBase<T1> & y, const Eigen::MatrixBase<T3> & x, const Eigen::MatrixBase<T2> & K, 
              int EMsteps, int EMsteps_fail, double EM_alpha, bool constraint, double min_s2, double min_tau, int max_iter, 
              double eps, bool verbose, Vector2d & theta, double & logL, double & logL0, int & niter, 
              double & gr_norm, VectorXd & Py, VectorXd & KPy, VectorXd & beta, MatrixXd & XViX_i, double & varXbeta, bool start_theta) {

  int n(y.rows()), p(x.cols());
 
  MatrixXd V(n,n);
  MatrixXd Vi(n,n), P(n,n);
  MatrixXd XViX(p,p);
  MatrixXd ViX(n,p);
  VectorXd PPy(n), PKPy(n);

  // X'X
  MatrixXd xtx( MatrixXd(p,p).setZero().selfadjointView<Lower>().rankUpdate( x.transpose() ));
  MatrixXd xtxi(p,p); // et son inverse
  double det_xtx, ldet_xtx;
  MatrixXd xtx0(xtx);
  sym_inverse(xtx0, xtxi, ldet_xtx, det_xtx, 1e-5); // détruit xtx0

  // Calcul de log L0 et
  // choix paramètres initiaux
  VectorXd xty = x.transpose() * y.col(0);
  double s2_0 = (y.col(0).dot(y.col(0)) - xty.dot( xtxi*xty ))/(n-p);
  logL0 = -0.5*((n-p)*log(s2_0) + ldet_xtx + (n-p));
 
  Vector2d theta0, gr;
  Matrix2d AI;
  double log_detV, detV, old_logL, d, log_d;

  double mean_diag = K.diagonal().mean();
  if(!start_theta) {
    theta(0) =  s2_0/2; // s2
    theta(1) =  s2_0/2/mean_diag; // tau
  }
  //------------


  bool bloc_tau = false, bloc_s2 = false;
  bool EM = false;

  gr_norm = eps+1;
  int i;
  for(i = 0; i < max_iter; i++) {
    if(verbose) Rcout << "[Iteration " << i+1 << "] theta = " << theta.transpose() << "\n";    

    V = theta(0)*MatrixXd::Identity(n,n) + theta(1)*K;

    // Calcul de Vi = inverse(V)
    sym_inverse(V,Vi,log_detV,detV,1e-7);

    // Calcul de P
    ViX.noalias() = Vi * x;
    XViX.noalias() = x.transpose() * ViX;
    sym_inverse(XViX, XViX_i, log_d, d, 1e-5);
    P.noalias() = Vi - ViX * XViX_i * ViX.transpose();
 
    // Optimiser les produits pour tenir 
    // compte de la symmétrie de P [ça change peanuts]
    Py.noalias()   =  P.selfadjointView<Lower>() * y;
    old_logL = logL;
    logL = -0.5*(log_detV + log_d + Py.dot(y.col(0)));
    if(verbose) Rcout << "[Iteration " << i+1 << "] log L = " << logL << "\n";


    // Is new value of likelihood OK ?
    if(i > 0 &&  logL < old_logL) {
      if(EM) {
        Rcout << "EM step failed to improve likelihood (this should not happen)\n";
      }
      else {
        EMsteps = EMsteps_fail;
        if(verbose) Rcout << "[Iteration " << i+1 << "] AI algorithm failed to improve likelihood, trying " 
                          << EMsteps << "EM step\n";
        theta = theta0;
        continue;
      }
    }

    // gradient
    KPy.noalias()  = K * Py; // le .selfadjointView ne compile pas avec le template !!
    PPy.noalias()  = P.selfadjointView<Lower>() * Py;
    PKPy.noalias() = P.selfadjointView<Lower>() * KPy;

    gr(0) = -0.5*(P.trace() - Py.squaredNorm());
    gr(1) = -0.5*(trace_of_product(K,P) - Py.dot(KPy));

    //update theta
    theta0 = theta;
    if(EMsteps > 0) {
      theta(0) = theta0(0) + 2*EM_alpha*theta0(0)*theta0(0)/n*gr(0);
      theta(1) = theta0(1) + 2*EM_alpha*theta0(1)*theta0(1)/n*gr(1);
      // logL = old_logL;
      if(verbose) Rcout << "[Iteration " << i+1 << "] EM update" << "\n";
      EM = true;
      EMsteps--;
    } else {

      if(constraint) { // on débloque les paramètres si le gradient ne pointe plus hors de la boîte
        if(bloc_s2  && gr(0) > 0) bloc_s2  = false; 
        if(bloc_tau && gr(1) > 0) bloc_tau = false;
      }

      // Average Information
      AI(0,0) = 0.5*PPy.dot(Py);
      AI(1,0) = AI(0,1) = 0.5*PPy.dot(KPy);
      AI(1,1) = 0.5*PKPy.dot(KPy);

      // Rcout << "\n unconstrained AI =\n" << AI << "\n------\n";

      if(constraint && bloc_s2) {
        theta(1) += gr(1)/AI(1,1);
        gr(0) = 0;
      }
      else if(constraint && bloc_tau) {
        theta(0) += gr(0)/AI(0,0);
        gr(1) = 0;
      }
      else {
        theta += AI.inverse()*gr;
      }

      if(constraint) {
        if(theta(0)  < min_s2)  { theta(0) = min_s2;  bloc_s2  = true; }
        if(theta(1)  < min_tau) { theta(1) = min_tau; bloc_tau = true; }
      }

      if(verbose) Rcout << "[Iteration " << i+1 << "] AI-REML update" << "\n";
      EM = false;
    }

    gr_norm = gr.norm();
    if(verbose) Rcout << "[Iteration " << i+1 << "] ||gradient|| = " << gr_norm << "\n";
    // Rcout << "theta avant contrainte = " << theta.transpose() << "\n";
    // Rcout << "nouveau theta = " << theta.transpose() << "\n";
    if(gr_norm < eps) {
      logL += gr.dot(theta-theta0);  // update linéaire de logL avant de sortir...
      break; 
    }
    checkUserInterrupt();
  }
  niter = i+1;
  // Rque : l'utilisateur récupère Py qui est utile dans calcul des BLUP 
  // tau Py correspond à v dans la formulation Y = K v + e avec v ~ N(0, tau K+)
  // [l'utilisateur récupère celui qui est calculé avec le Py qui correspond à theta0 !! tant pis pour lui]
  // l'utilisateur récupère aussi logL

  // BLUP POUR beta
  // X beta = Y - omega - sigma2 * Py
  // beta = (X'X)^{-1} X'(Y - omega - sigma² Py)
  beta = x.transpose() * (y - theta(1)*KPy - theta(0)*Py);
  beta = xtxi * beta;

  // Calcul débiaisé de var (X \hat beta)
  double psi = trace_of_product(xtx, XViX_i)/(n-1);
  VectorXd gg(x.transpose()*VectorXd::Ones(n));
  psi -= gg.dot(XViX_i*gg)/n/(n-1);

  VectorXd Xbeta = x*beta;
  double SXb = Xbeta.sum();
  varXbeta = (Xbeta.squaredNorm() - SXb*SXb/n)/(n-1) - psi;
}


