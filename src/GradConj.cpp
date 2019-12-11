#include "math/Sparse.h"
#include "math/Vector.h"
#include <cstdlib>
#include <cmath>

void GC_sparse_parallel(Sparse const& A,
                        Vector const& b,
                        Vector& x,
                        Vector const& x0,
                        int kmax,
                        double eps) {
    double alpha, beta, gamma;

    Vector r = b;
    Vector p = r;
    Vector prod;
    if (x0.size() != 0)
        prod = A * x0;
    r -= prod;

    /* Computed here and updated at the end of the loop (avoid twice computing the same thing) */
    double rdotr = r.dot(r);
    beta = rdotr;
    for (size_t k = 0; k < kmax && beta > eps; ++k) {
        spmv(1., A, prod, 0., prod);
        alpha = rdotr / p.dot(prod);

        /* x <- x + alpha * p */
        axpy(alpha, p, x);
        /* r <- r - alpha * A * p */
        r -= prod.scale(alpha);

        double rndotrn = r.dot(r);
        /* gamma <- (rn * rn)/(r * r) */
        gamma = rndotrn / rdotr;

        /* r <- r + beta * p */
        p.scale(gamma);
        p += r;

        beta = sqrt(rdotr);
        rdotr = rndotrn;
    }
}