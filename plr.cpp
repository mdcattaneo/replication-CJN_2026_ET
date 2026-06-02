#include <Rcpp.h>
using namespace Rcpp;

// ============================================================
// PLR pairwise difference estimator
// Robinson (1988) DGP, scalar x (k=1).
//
// Model:  y_i = x_i * theta_0 + gamma_0(w_i) + eps_i
//
// Closed-form estimator:
//   theta_hat = SXy / SXX
//   SXX = sum_{i<j} K_h(w_i - w_j) * (x_i - x_j)^2
//   SXy = sum_{i<j} K_h(w_i - w_j) * (x_i - x_j) * (y_i - y_j)
//
// biweight kernel K(u) = (15/16)(1-u^2)^2 * 1(|u|<=1)
// For d>1, we use the product kernel
// ============================================================


// [[Rcpp::export]]
List generate_data_plr(int n, double theta0 = 1.0, int model = 1) {

  NumericVector y(n), x(n);
  NumericMatrix w;

  const double sqrt2 = std::sqrt(2.0);

  if (model == 1 ) {
    // ---- Model (6.1): d = 1 ----
    // Cholesky of Sigma = [[4,2],[2,3]]:  L = [[2,0],[1,sqrt(2)]]
    // x = 2*z1,  w = z1 + sqrt(2)*z2
    // => Var(x)=4, Var(w)=3, Cov(x,w)=2  ✓
    w = NumericMatrix(n, 1);

    for (int i = 0; i < n; i++) {
      double z1 = R::rnorm(0.0, 1.0);
      double z2 = R::rnorm(0.0, 1.0);

      x[i]    = 2.0 * z1;
      w(i, 0) = z1 + sqrt2 * z2;

      double gamma = (1.0 + w(i, 0) * w(i, 0));
      y[i] = theta0 * x[i] + gamma + R::rnorm(0.0, 1.0);
    }

  } else if (model == 2 || model == 3) {
    // ---- Models 2&3 modified (6.2), d=2, d=3 ----
    //
    // Both use the same equicorrelated factorisation:
    //   x   = 1 + sqrt(2)*z0 + eta_x
    //   z_j = 1 + sqrt(2)*z0 + eta_j,  j = 1,...,q
    // z0, eta_x, eta_1,...,eta_q iid N(0,1)
    // => E[.]=1, Var[.]=2+1=3, Cov(any pair)=2, Corr=2/3  ✓
    // gamma_0(w) = 1 + sum_{j=1}^q z_j^2   (alpha=gamma_j=1)
    //
    // The DGP mirrors in (6.2) Robinson (1988)..
    const int q = (model == 2) ? 2 : 3;
    w = NumericMatrix(n, q);

    for (int i = 0; i < n; i++) {
      double z0   = R::rnorm(0.0, 1.0);           // common factor
      x[i]        = 1.0 + sqrt2 * z0 + R::rnorm(0.0, 1.0);

      double gamma = 1.0;                           // alpha = 1
      for (int j = 0; j < q; j++) {
        double zj = 1.0 + sqrt2 * z0 + R::rnorm(0.0, 1.0);
        w(i, j)   = zj;
        gamma     += zj * zj;                       // gamma_j = 1
      }
      y[i] = theta0 * x[i] + gamma + R::rnorm(0.0, 1.0);
    }
  } else {
    Rcpp::stop("model must be 1, 2, or 3.");
  }

  return List::create(Named("y") = y, Named("x") = x, Named("w") = w);
}


// Biweight (quartic) kernel
inline double biweight(double u) {
  if (std::abs(u) <= 1.0) {
    double t = 1.0 - u * u;
    return 15.0 / 16.0 * t * t;
  }
  return 0.0;
}


// [[Rcpp::export]]
double estimate_plr(NumericVector y,
                    NumericVector x,
                    NumericMatrix w,   // n x d  (d=1 or d=5)
                    double        h) {
  //
  // Closed-form WLS (CJN Section 2.1):
  //   theta_hat = SXy / SXX
  //
  // Kernel: product over d dimensions.
  //   K_h(w_i - w_j) = prod_{l=1}^d biweight((w_{il} - w_{jl}) / h)
  //
  // The h^{-d} normalisation cancels in the ratio SXy / SXX and is omitted.
  //
  int    n   = y.size();
  int    d   = w.ncol();
  double SXX = 0.0;
  double SXy = 0.0;

  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {

      // Product kernel: exit early as soon as any factor is zero
      double kern = 1.0;
      for (int l = 0; l < d; l++) {
        kern *= biweight((w(i, l) - w(j, l)) / h);
        if (kern == 0.0) break;
      }
      if (kern <= 0.0) continue;

      double dx  = x[i] - x[j];
      double dy  = y[i] - y[j];

      SXX += kern * dx * dx;
      SXy += kern * dx * dy;
    }
  }

  if (SXX < 1e-300) return R_NaN;
  return SXy / SXX;
}


// [[Rcpp::export]]
NumericVector estimate_plr_grid(NumericVector y,
                                NumericVector x,
                                NumericMatrix w,
                                NumericVector hvals) {
  //
  // Evaluate the PLR estimator at every bandwidth in hvals in a single
  // pass over the O(n^2) pairs.  hvals must be sorted ascending.
  //
  // For each pair (i < j):
  //   1. Compute max_l |w_il - w_jl|.  If >= h_max, the pair contributes
  //      to no bandwidth — skip immediately.
  //   2. Binary-search for the smallest h index where the pair first
  //      contributes (hvals[hh] > max_l |w_il - w_jl|).
  //   3. For all hh from that index onward, compute the product kernel
  //      and accumulate into SXX[hh] and SXy[hh].
  //
  // This avoids repeating the O(n^2) scan once per bandwidth, cutting
  // total work from O(n^2 * H) independent scans to a single O(n^2)
  // scan with O(H) accumulation per contributing pair.
  //
  int n    = y.size();
  int d    = w.ncol();
  int H    = hvals.size();

  NumericVector SXX(H, 0.0);
  NumericVector SXy(H, 0.0);

  double h_max = hvals[H - 1];

  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {

      // Pass 1: max component distance, early exit at h_max
      double max_dist = 0.0;
      for (int l = 0; l < d; l++) {
        double adw = std::abs(w(i, l) - w(j, l));
        if (adw >= h_max) { max_dist = h_max; break; }
        if (adw > max_dist) max_dist = adw;
      }
      if (max_dist >= h_max) continue;

      // First bandwidth index where the pair contributes
      // (linear scan; H is small, typically < 100)
      int hh0 = 0;
      while (hh0 < H && hvals[hh0] <= max_dist) hh0++;
      if (hh0 >= H) continue;

      double dx   = x[i] - x[j];
      double dy   = y[i] - y[j];
      double dxdx = dx * dx;
      double dxdy = dx * dy;

      // Pass 2: accumulate for each contributing bandwidth
      for (int hh = hh0; hh < H; hh++) {
        double h    = hvals[hh];
        double kern = 1.0;
        for (int l = 0; l < d; l++) {
          double u = std::abs(w(i, l) - w(j, l)) / h;
          double t = 1.0 - u * u;
          kern *= (15.0 / 16.0) * t * t;
        }
        SXX[hh] += kern * dxdx;
        SXy[hh] += kern * dxdy;
      }
    }
  }

  NumericVector theta(H);
  for (int hh = 0; hh < H; hh++)
    theta[hh] = (SXX[hh] < 1e-300) ? R_NaN : SXy[hh] / SXX[hh];

  return theta;
}
