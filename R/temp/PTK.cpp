#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector cpp_PTK(NumericVector Y, NumericVector Pk, int Nmin) {
  // int n = x.size();
  // IntegerVector ind =  seq(1, n);
  // IntegerVector ind =  seq_along(Y);
  // NumericVector x2 = Y[ind <= k];
  double PTx = -9999.9;
  double KPx = 0;
  double Tx;

  int N = Y.size();
  
  double EY1; 
  double EY2;
  double ALL = sum(Y);

  double mean, var, std;

  for (int k = Nmin; k <= N-Nmin; k++) {
    EY1 = 0;
    for(int i = 0; i < k; i++) EY1 += Y[i];
    EY1 = EY1 / k;
    // EY2 = 0;
    // double EY1 = mean(Y[1:k]);
    EY2 = (ALL - EY1 * k) / (N - k);
    
    var = 0;
    for(int i = 0; i < N; i++) {
      mean = i < k ? EY1 : EY2;
      var += pow(Y[i] - mean, 2);
    }

    std = sqrt(var / (N - 2));
    double Tk = sqrt(k * (N - k) * 1.0 / N) * abs(EY1 - EY2) / std;
    double PTk = Tk * Pk[k - 1];
    
    // Rprintf("var = %f, sd = %f, Tk = %f, PTk = %f \n", var, std, Tk, PTk);
    if (PTk > PTx) {
      PTx = PTk;
      KPx = k;
      Tx = Tk;
    }
  }

  return NumericVector::create(
    Named("PTx") = PTx, 
    Named("KPx") = KPx, 
    Named("Tx")  = Tx);
}

// PTK <- function(Y, Pk, Nmin) {
//   #  search input vector, return PTx.max and corresponding KPx
//   PTx <- (-9999.9)
//   KPx <- 0
//   N <- length(Y)
//   ALL <- sum(Y)

//   for (k in Nmin:(N - Nmin)) {
//     EY1 <- mean(Y[1:k])
//     EY2 <- (ALL - EY1*k) / (N - k) # mean(Y[(k + 1):N])

//     var <- sum(c((Y[1:k] - EY1)^2, (Y[(k + 1):N] - EY2)^2))
//     std <- sqrt(var / (N - 2))

//     Tk <- sqrt(k * (N - k) / N) * abs(EY1 - EY2) / std
//     PTk <- Tk * Pk[k]
//     if (PTk > PTx) {
//       PTx <- PTk
//       KPx <- k
//       Tx <- Tk
//     }
//   }
//   list(PTx = PTx, KPx = KPx, Tx = Tx)
// }
