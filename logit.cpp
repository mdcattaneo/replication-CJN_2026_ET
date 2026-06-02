// logit.cpp
//
// C++ (Rcpp) functions for the pairwise logit simulation.
// Exports to R:
//   generate_data           — simulate data from DGP models 1–3
//   precompute_active_pairs — collect kernel-weighted discordant pairs
//   fast_wlogit             — Newton-Raphson weighted logistic regression (no intercept)
//   gradient_fnc            — gradient of the pairwise logistic objective
//   objective_fnc           — value of the pairwise logistic objective

#include <Rcpp.h>
#include <algorithm>
#include <numeric>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame generate_data(int n,
                   NumericVector beta = NumericVector::create(1.0, 1.0),
                   int model = 1) {

  NumericVector v(n), x1(n), x2(n), g(n), eps(n), y(n);

  if (model == 1) {
    // w ~ N(0,1), x1 = v + w^2, g(w) = w^2 - 2
    NumericVector w(n);
    for (int i = 0; i < n; i++) {
      w[i]   = R::rnorm(0.0, 1.0);
      v[i]   = R::rnorm(0.0, 1.0);
      x2[i]  = (R::runif(0.0, 1.0) < 0.5) ? -1.0 : 1.0;
      x1[i]  = v[i] + w[i] * w[i];
      g[i]   = w[i] * w[i] - 2.0;
      eps[i] = R::rlogis(0.0, 1.0);
      double y_star = beta[0]*x1[i] + beta[1]*x2[i] + g[i] + eps[i];
      y[i] = (y_star >= 0.0);
    }
    return DataFrame::create(
      Named("y")  = y,
      Named("x1") = x1,
      Named("x2") = x2,
      Named("w")  = w
    );

  } else if (model == 2) {
    // (w1,w2) ~ BivNormal(0, Sigma), Sigma = [[1,0.2],[0.2,1]]
    // Cholesky: w1=z1, w2=0.2*z1 + sqrt(0.96)*z2
    // x1 = w1^2 + w2^2 + v, g(w) = w1^2 + w2^2 - 3
    const double L21 = 0.2, L22 = 0.9797958971;  // sqrt(0.96)
    NumericVector w1(n), w2(n);
    for (int i = 0; i < n; i++) {
      double z1 = R::rnorm(0.0, 1.0), z2 = R::rnorm(0.0, 1.0);
      w1[i]  = z1;
      w2[i]  = L21*z1 + L22*z2;
      v[i]   = R::rnorm(0.0, 1.0);
      x2[i]  = (R::runif(0.0, 1.0) < 0.5) ? -1.0 : 1.0;
      double ss = w1[i]*w1[i] + w2[i]*w2[i];
      x1[i]  = ss + v[i];
      g[i]   = ss - 3.0;
      eps[i] = R::rlogis(0.0, 1.0);
      double y_star = beta[0]*x1[i] + beta[1]*x2[i] + g[i] + eps[i];
      y[i] = (y_star >= 0.0);
    }
    return DataFrame::create(
      Named("y")  = y,
      Named("x1") = x1,
      Named("x2") = x2,
      Named("w1") = w1,
      Named("w2") = w2
    );

  } else {
    // model == 3
    // (w1,w2,w3) ~ TriNormal(0, Sigma), Sigma_{lm} = 0.2 (l!=m), 1 (l==m)
    // Cholesky of Sigma:
    //   L = [[1,0,0],[0.2,sqrt(0.96),0],[0.2, 0.16/sqrt(0.96), sqrt(1-0.04-(0.16^2/0.96))]]
    const double L21 = 0.2,  L22 = 0.9797958971;             // sqrt(0.96)
    const double L31 = 0.2,  L32 = 0.16 / 0.9797958971;      // 0.16/sqrt(0.96)
    const double L33 = 0.9660917831;                          // sqrt(1-0.04-L32^2)
    NumericVector w1(n), w2(n), w3(n);
    for (int i = 0; i < n; i++) {
      double z1 = R::rnorm(0.0,1.0), z2 = R::rnorm(0.0,1.0), z3 = R::rnorm(0.0,1.0);
      w1[i]  = z1;
      w2[i]  = L21*z1 + L22*z2;
      w3[i]  = L31*z1 + L32*z2 + L33*z3;
      v[i]   = R::rnorm(0.0, 1.0);
      x2[i]  = (R::runif(0.0, 1.0) < 0.5) ? -1.0 : 1.0;
      double ss = w1[i]*w1[i] + w2[i]*w2[i] + w3[i]*w3[i];
      x1[i]  = ss + v[i];
      g[i]   = ss - 4.0;
      eps[i] = R::rlogis(0.0, 1.0);
      double y_star = beta[0]*x1[i] + beta[1]*x2[i] + g[i] + eps[i];
      y[i] = (y_star >= 0.0);
    }
    return DataFrame::create(
      Named("y")  = y,
      Named("x1") = x1,
      Named("x2") = x2,
      Named("w1") = w1,
      Named("w2") = w2,
      Named("w3") = w3
    );
  }
}



// biweight (quartic) kernel
inline double biweight(double u) {
  if (std::abs(u) <= 1.0) {
    double t = 1.0 - u * u;
    return 15.0 / 16.0 * t * t;
  } else {
    return 0.0;
  }
}

// logistic function
inline double logistic(double x) {
  return 1.0 / (1.0 + std::exp(-x));
}


// Returns list with: kern_vals, dX (matrix of X[i,]-X[j,]), yi, scale, n_active
// Uses sorted w[,0] + binary search on first dimension, then checks remaining dimensions.
// Product kernel across all w columns.
// [[Rcpp::export]]
List precompute_active_pairs(NumericVector y, NumericMatrix X, NumericMatrix w, double h) {

  int n = y.size();
  int k = X.ncol();
  int d = w.ncol();
  double scale = (double)n * (n - 1) / 2.0 * std::pow(h, d);

  // Sort observation indices by w[,0] so we can use binary search on first dimension
  std::vector<int> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b){ return w(a, 0) < w(b, 0); });

  std::vector<double> ws0(n);
  for (int s = 0; s < n; s++) ws0[s] = w(ord[s], 0);

  // Pass 1: count active pairs
  int n_active = 0;
  for (int s = 0; s < n - 1; s++) {
    int hi = (int)(std::lower_bound(ws0.begin() + s + 1, ws0.end(), ws0[s] + h) - ws0.begin());
    int i  = ord[s];
    for (int t = s + 1; t < hi; t++) {
      int j = ord[t];
      if (y[i] == y[j]) continue;
      bool active = true;
      for (int l = 1; l < d; l++) {
        if (std::abs(w(i, l) - w(j, l)) >= h) { active = false; break; }
      }
      if (active) n_active++;
    }
  }

  // Allocate output storage
  NumericVector kern_vals(n_active);
  NumericMatrix dX(n_active, k);
  NumericVector yi(n_active);

  // Pass 2: fill active pairs
  int pos = 0;
  for (int s = 0; s < n - 1; s++) {
    int hi = (int)(std::lower_bound(ws0.begin() + s + 1, ws0.end(), ws0[s] + h) - ws0.begin());
    int i  = ord[s];
    for (int t = s + 1; t < hi; t++) {
      int j = ord[t];
      if (y[i] == y[j]) continue;
      double kern = biweight((w(i, 0) - w(j, 0)) / h);
      bool active = (kern > 0.0);
      for (int l = 1; l < d && active; l++) {
        double kl = biweight((w(i, l) - w(j, l)) / h);
        if (kl <= 0.0) active = false; else kern *= kl;
      }
      if (active) {
        kern_vals[pos] = kern;
        yi[pos]        = y[i];
        for (int l = 0; l < k; l++) dX(pos, l) = X(i, l) - X(j, l);
        pos++;
      }
    }
  }

  return List::create(
    Named("kern_vals") = kern_vals,
    Named("dX")        = dX,
    Named("scale")     = scale,
    Named("n_active")  = n_active,
    Named("yi")        = yi
  );
}

// Newton-Raphson weighted logistic regression (no intercept).
// start: warm-start beta (length k); maxit/tol: convergence controls.
// [[Rcpp::export]]
NumericVector fast_wlogit(NumericMatrix dX, NumericVector yi,
                          NumericVector kern_vals, NumericVector start,
                          int maxit = 25, double tol = 1e-8) {
  int n = dX.nrow(), k = dX.ncol();
  NumericVector beta = clone(start);
  if (beta.size() != k) { beta = NumericVector(k, 0.0); }
  if (n == 0) return beta;

  for (int iter = 0; iter < maxit; iter++) {
    // Accumulate gradient and upper-triangle of Hessian
    double g0 = 0.0, g1 = 0.0;
    double H00 = 0.0, H01 = 0.0, H11 = 0.0;

    for (int i = 0; i < n; i++) {
      double eta = dX(i, 0) * beta[0] + dX(i, 1) * beta[1];
      double q   = 1.0 / (1.0 + std::exp(-eta));
      double r   = kern_vals[i] * (yi[i] - q);   // score contribution
      double s   = kern_vals[i] * q * (1.0 - q); // info contribution
      g0  -= r * dX(i, 0);
      g1  -= r * dX(i, 1);
      H00 += s * dX(i, 0) * dX(i, 0);
      H01 += s * dX(i, 0) * dX(i, 1);
      H11 += s * dX(i, 1) * dX(i, 1);
    }

    double det = H00 * H11 - H01 * H01;
    if (std::abs(det) < 1e-15) break;

    double step0 = ( H11 * g0 - H01 * g1) / det;
    double step1 = (-H01 * g0 + H00 * g1) / det;
    beta[0] -= step0;
    beta[1] -= step1;

    if (std::sqrt(step0 * step0 + step1 * step1) < tol) break;
  }
  return beta;
}


// Gradient of the pairwise logistic objective (normalised by scale).
// Convention: eta = -dX * beta, where dX = X[i,] - X[j,] for the discordant
// pair (i,j) with y[i] = 1.
// [[Rcpp::export]]
NumericVector gradient_fnc(NumericVector beta,
                           NumericMatrix dX,
                           NumericVector kern,
                           NumericVector yi,
                           double scale){
  int n = yi.size();
  int k = dX.ncol();

  NumericVector grad(k);
  for (int i =0; i < n; i++){
    double eta = 0;
    for (int l = 0; l< k; l++){
      eta += -1*dX(i,l)*beta[l];
    }
    double p = 1.0 / (1.0 + std::exp(-eta));
    double coeff = kern[i] * (p - 1 + yi[i]);
    
    for (int l = 0; l < k; l++) {
      grad[l] += -1*coeff * dX(i,l);
    }
  }
  return grad/scale;
}

// Value of the pairwise logistic objective (normalised by scale).
// [[Rcpp::export]]
double objective_fnc(NumericVector beta,
                     NumericMatrix dX,
                     NumericVector kern,
                     NumericVector yi,
                     double scale) {
  int n = dX.nrow();
  int k = dX.ncol();
  double val = 0.0;
  
  for (int i = 0; i < n; i++){
    double eta = 0.0;
    for (int l = 0; l< k; l++){
      eta += -1*dX(i,l)*beta[l];
    }
    val += kern[i]*std::log1p(std::exp( eta*(2*yi[i]-1) ));
  }
  return val/scale;
}
