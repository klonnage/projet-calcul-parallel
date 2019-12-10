#include "Sparse.h"
#include "Vector.h"
#include <cblas.h>

Sparse::Sparse (int Nx, int Ny, double alpha, double beta, double gamma) :
Nx(Nx), Ny(Ny), alpha(alpha), beta(beta), gamma(gamma)
{}

Sparse::Sparse (Sparse const& A) : Nx(Nx), Ny(Ny), alpha(alpha), beta(beta), gamma(gamma)
{}

Sparse::~Sparse() {}

Vector Sparse::operator*(Vector const& v) const{
  Vector res(v.size);
  /* The following does not create problem as x won't be modified */
  spmv(1., *this, const_cast<Vector&>(v), 0., res);
}

void spmv(double a, Sparse const& A, Vector& x, double b, Vector& y) {
  y.scale(b);

  /* Scale by the main diag of A */
  cblas_daxpy(x.size(), a * A.alpha, x.data(), 1, y.data(), 1);

  /* Scale by sub diagonal of A */
  cblas_daxpy(x.size() - 1, a * A.beta, x.data() + 1, 1, y.data(), 1);
  cblas_daxpy(x.size() - 1, a * A.beta, x.data(), 1, y.data() + 1, 1);

  /* Scale by gamma */
  cblas_daxpy(x.size() - A.Ny, a * A.gamma, x.data() + A.Ny, 1, y.data(), 1);
  cblas_daxpy(x.size() - A.Ny, a * A.gamma, x.data(), 1, y.data() + A.Ny, 1);
}