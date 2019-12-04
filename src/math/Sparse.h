#pragma once
#include "Vector.h"

class Sparse {
    double alpha, beta, gamma;
    double N;
public:
    Sparse () = delete;
    Sparse (double N, double alpha, double beta, double gamma);
    Sparse (Sparse const& A);
    ~Sparse();

    friend void gemv(double a, Sparse const& A, Vector& x);

    Vector operator*(Vector const& v) const;
};