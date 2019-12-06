#pragma once
#include "Vector.h"

class Sparse {
    double alpha, beta, gamma;
    int Nx, Ny;
public:
    Sparse () = delete;
    Sparse (int Nx, int Ny, double alpha, double beta, double gamma);
    Sparse (Sparse const& A);
    ~Sparse();

    friend void spmv(double a, Sparse const& A, Vector& x, double b, Vector& y);

    Vector operator*(Vector const& v) const;
};