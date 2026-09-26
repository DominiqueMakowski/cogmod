# Nanoseconds per call of the transcendentals the LNR density needs, with the
# RTools toolchain CmdStan uses. Arguments vary so nothing is hoisted; the sum
# is returned so nothing is optimised away.
Rcpp::cppFunction(includes = "#include <cmath>\n#include <chrono>", code = '
Rcpp::NumericVector bench(int n) {
  std::vector<double> x(1024);
  for (int i = 0; i < 1024; i++) x[i] = -3.0 + 6.0 * i / 1024.0;
  const char* names[] = {"add", "log", "exp", "erfc", "log1p", "sqrt"};
  Rcpp::NumericVector out(6);
  volatile double sink = 0;
  for (int f = 0; f < 6; f++) {
    double s = 0;
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int r = 0; r < n; r++) {
      for (int i = 0; i < 1024; i++) {
        double v = x[i] + 1e-9 * r;
        switch (f) {
          case 0: s += v + 1.0; break;
          case 1: s += std::log(v + 4.0); break;
          case 2: s += std::exp(v); break;
          case 3: s += std::erfc(v); break;
          case 4: s += std::log1p(v + 3.5); break;
          case 5: s += std::sqrt(v + 4.0); break;
        }
      }
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    sink = sink + s;
    out[f] = std::chrono::duration<double, std::nano>(t1 - t0).count() / (1024.0 * n);
  }
  out.attr("names") = Rcpp::CharacterVector(names, names + 6);
  return out;
}')
print(round(bench(20000L), 1))
