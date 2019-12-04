#include "Sparse.h"

Sparse::Sparse (double N, double alpha, double beta, double gamma) :
N(N), alpha(alpha), beta(beta), gamma(gamma)
{}

Sparse::Sparse (Sparse const& A) : N(N), alpha(alpha), beta(beta), gamma(gamma)
{}

Sparse::~Sparse() {}

void gemv(double a, Sparse const& A, Vector& x) {
  
}