#include "math/Sparse.h"
#include "math/Vector.h"
#include <cstdlib>
#include <cmath>
#include <iostream>

void GC_sparse_parallel(Sparse const& A,
                        Vector const& b,
                        Vector& x,
                        Vector const& x0,
                        int kmax,
                        double eps) {
    double alpha, beta, gamma;

    Vector r1, r2(b.size());
    r1 = b;
    Vector p;
    Vector prod;
    if (x0.size() != 0)
        prod = A * x0;
    r1 -= prod;
    
    p = r1;
    /* Computed here and updated at the end of the loop (avoid twice computing the same thing) */
    beta = r1.dot(r1);
    for (int k = 0; k < kmax && beta > eps; ++k) {
        double rdotr = r1.dot(r1);
        spmv(1., A, p, 0., prod);
        std::cout << "spmv : "  << prod << std::endl;
        alpha = rdotr / p.dot(prod);

        /* x <- x + alpha * p */
        axpy(alpha, p, x);
        /* r <- r - alpha * A * p */
        r2 = r1;
        axpy(-alpha, prod, r2);

        double rndotrn = r2.dot(r2);
        /* gamma <- (rn * rn)/(r * r) */
        gamma = rndotrn / rdotr;

        /* r <- r + beta * p */
        p.scale(gamma);
        p += r2;

        beta = sqrt(rndotrn);
        r1 = r2;
    }
}
