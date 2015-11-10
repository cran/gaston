// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>

using namespace Rcpp;
using namespace Eigen;

typedef Map<MatrixXd> Map_MatrixXd;

template<typename T1, typename T2>
void AIREML_nofix(const Eigen::MatrixBase<T1> & y, const Eigen::MatrixBase<T2> & K, int EMsteps, int EMsteps_fail,
                  double EM_alpha, bool constraint, double min_s2, double min_tau, int max_iter, double eps,
                  bool verbose, Vector2d & theta, double & logL, double & logL0, int & niter, double & gr_norm, 
                  VectorXd & Py, VectorXd & KPy, bool start_theta) {
  int n(y.rows());
 
  double var_y = y.squaredNorm()/n; // eh oui, pas d'effet fixe : Y est supposé centré !!

  // if(verbose) Rcout << "var(Y) = " << var_y << "\n";
  MatrixXd V(n,n);
  MatrixXd P(n,n), Delta(n,n);
  VectorXd PPy(n), PKPy(n);
  
  Vector2d theta0, gr;
  Matrix2d AI;
  double log_detV, detV, old_logL;

  // Choix paramètres initiaux
  double mean_diag = K.diagonal().mean();
  if(!start_theta) {
    theta(0) =  var_y/2; // s2
    theta(1) =  var_y/2/mean_diag; // tau
  }

  bool bloc_tau = false, bloc_s2 = false;
  bool EM = false;

  gr_norm = eps+1;
  int i;
  for(i = 0; i < max_iter; i++) {
    if(verbose) Rcout << "[Iteration " << i+1 << "] theta = " << theta.transpose() << "\n";    

    V = theta(0)*MatrixXd::Identity(n,n) + theta(1)*K;

    // Calcul de P = inverse(V)
    sym_inverse(V,P,log_detV,detV,1e-7);

    // Optimiser les produits pour tenir 
    // compte de la symmétrie de P [ça change peanuts]
    Py.noalias()   =  P.selfadjointView<Lower>() * y;
    old_logL = logL;
    logL = -0.5*(log_detV + Py.dot(y.col(0)));
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
     // if(verbose) Rcout << "[Iteration " << i+1 << "] EM step ok, going back to AI algorithm\n"; 
   
    // computing gradient
    KPy.noalias()  = K * Py; // le .selfadjointView ne compile pas avec le template !!
    PPy.noalias()  = P.selfadjointView<Lower>() * Py;
    PKPy.noalias() = P.selfadjointView<Lower>() * KPy;

    gr(0) = -0.5*(P.trace() - Py.squaredNorm());
    gr(1) = -0.5*(trace_of_product(K,P) - Py.dot(KPy));
    // Rcout << "\n unconstrained gr = " << gr.transpose() << "\n";
    // updating theta

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
  // Rque : l'utilisateur récupère Py qui est utile dans calcul des BLUP 
  // tau Py correspond à v dans la formulation Y = K v + e avec v ~ N(0, tau K+)
  // [l'utilisateur récupère celui qui est calculé avec le Py qui correspond à theta0 !! tant pis pour lui]
  // l'utilisateur récupère aussi logL
  logL0 = -0.5*n*(log(var_y)+1);
  niter = i+1;
}


