#include "math/Sparse.h"
#include "math/Vector.h"

void GC_sparse_parallel(Sparse const* A,
                        Vector const* b,
                        Vector* x,
                        Vector const* x0,
                        int kmax,
                        double eps) {
    Vector r = *b;
    Vector prod;
    if (x0)
        prod = (*A) * (*x0);
    r -= prod;
}