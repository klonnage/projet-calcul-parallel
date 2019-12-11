#pragma once
#include "math/Sparse.h"
#include "math/Vector.h"
#include <cstdlib>
#include <cmath>

void GC_sparse_parallel(Sparse const& A,
                        Vector const& b,
                        Vector& x,
                        Vector const& x0,
                        int kmax,
                        double eps);