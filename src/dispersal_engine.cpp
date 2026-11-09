#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix cpp_accumulate_dispersal(NumericMatrix A,
                                       IntegerMatrix cells,
                                       IntegerVector idx,
                                       NumericVector theta,
                                       NumericVector rad,
                                       bool keep_origin) {
  const R_xlen_t n = idx.size();
  const int nr = A.nrow();
  const int nc = A.ncol();
  const int ncell = cells.nrow();

  double* pa = A.begin();
  const int* pc = cells.begin();
  const int* pi = idx.begin();
  const double* pt = theta.begin();
  const double* pr = rad.begin();

  for (R_xlen_t i = 0; i < n; ++i) {
    const int k = pi[i] - 1;
    const int r0 = pc[k] - 1;
    const int c0 = pc[k + ncell] - 1;   // segunda coluna da matriz cells

    // R::fround(x, 0) = round() do R (half to even), nao std::round
    const int d_lat = (int) R::fround(pr[i] * std::sin(pt[i]), 0.0);
    const int d_lon = (int) R::fround(pr[i] * std::cos(pt[i]), 0.0);

    const int r = r0 + d_lat;
    const int c = c0 + d_lon;

    if (r >= 0 && r < nr && c >= 0 && c < nc) {
      pa[r + nr * c] += 1.0;
    } else if (keep_origin) {
      pa[r0 + nr * c0] += 1.0;
    }
  }
  return A;
}
