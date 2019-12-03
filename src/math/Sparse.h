#pragma once
#include "Vector.h"

class Sparse {
    double alpha, beta, gamma;
    double Nx, Lx;
public:
    Sparse () = delete;
    Sparse (double Nx, double Lx);
    ~Sparse();

    Vector operator*(Vector const& v) const;
};