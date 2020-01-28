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
    if (x0.size() != 0) {
        prod = A * x0;
        r1 -= prod;
        x = (x0);
    } else {
        x.set_value(1.);
    }
    
    p = r1;
    beta = r1.dot(r1);
    int k;
    for (k = 0; k < kmax; ++k) {
        double rdotr = r1.dot(r1);
        spmv(1., A, p, 0., prod);
        alpha = rdotr / p.dot(prod);

        /* x <- x + alpha * p */
        axpy(alpha, p, x);
        r2 = r1;
        /* r <- r - alpha * A * p */
        axpy(-alpha, prod, r2);

        double rndotrn = r2.dot(r2);
        beta = sqrt(rndotrn);
        //std::cout << "beta = " << beta << std::endl;
        if (beta <= eps) {
            break;
        }
        /* gamma <- (rn * rn)/(r * r) */
        gamma = rndotrn / rdotr;

        /* r <- r + beta * p */
        p.scale(gamma);
        p += r2;

        r1 = r2;
    }
    std::cout << k << std::endl;
    std::cout << beta << std::endl;
}
