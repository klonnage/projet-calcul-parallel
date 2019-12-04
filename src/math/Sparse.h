#pragma once
#include "Vector.h"

class Sparse {
    double alpha, beta, gamma;
    double N;
public:
    Sparse () = delete;
    Sparse (double N, double alpha, double beta, double gamma);
    ~Sparse();

    Vector operator*(Vector const& v) const;
};