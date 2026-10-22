// Standard Jacobi polynomials P_n^(alpha, beta) on [-1, 1].

#pragma once

namespace wasl {

// P_n^(alpha, beta)(x), orthogonal with respect to the weight
// (1 - x)^alpha (1 + x)^beta. Base cases P_0 = 1 and
// P_1 = ((alpha - beta) + (alpha + beta + 2) x) / 2, then the three-term
// recurrence (Abramowitz-Stegun 22.7.1). The recurrence coefficients stay
// nonzero for the nonnegative alpha, beta used by the Dubiner basis.
inline double jacobi_p(double x, int n, double alpha, double beta) noexcept {
  if (n <= 0)
    return 1.0;
  double p_prev = 1.0;                                             // P_0
  double p_curr = 0.5 * (alpha - beta + (alpha + beta + 2.0) * x); // P_1
  for (int k = 1; k < n; ++k) {
    const double ab = 2.0 * k + alpha + beta;
    const double a = 2.0 * (k + 1) * (k + alpha + beta + 1.0) * ab;
    const double b = (ab + 1.0) * (ab * (ab + 2.0) * x + (alpha * alpha - beta * beta));
    const double c = 2.0 * (k + alpha) * (k + beta) * (ab + 2.0);
    const double p_next = (b * p_curr - c * p_prev) / a;
    p_prev = p_curr;
    p_curr = p_next;
  }
  return p_curr;
}

// d/dx P_n^(alpha,beta)(x) = (n + alpha + beta + 1) / 2 * P_{n-1}^(alpha+1,beta+1)(x),
// and zero for the degree-zero mode.
inline double jacobi_p_derivative(double x, int n, double alpha, double beta) noexcept {
  if (n <= 0)
    return 0.0;
  return 0.5 * (n + alpha + beta + 1.0) * jacobi_p(x, n - 1, alpha + 1.0, beta + 1.0);
}

} // namespace wasl
