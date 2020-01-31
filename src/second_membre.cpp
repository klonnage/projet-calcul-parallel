#include "second_membre.h"
#include "Vector.h"
#include <iostream>

void calcul_second_membre(Vector &F, int Nlime, int Ncol, double dx, double dy, double D, double dt, Vector const& g, Vector const& h, Vector const& termeSource) {
  for (int i = 0; i < Nlime; ++i) {
    for (int j = 0; j < Ncol; ++j) {
      F[i*Ncol + j] = termeSource[i*Ncol + j];
      /* First row */
      if(i == 0) {
        F[j] += g[j]*D/(dy*dy);
        if (j == 0) {
          F[0] += h[0]*D/(dx*dx);
        } else if (j == Ncol-1) {
          F[j] += h[Nlime]*D/(dx*dx);
        }
      }
      /* Last row */
      else if (i == Nlime-1) {
        F[i*Ncol + j] += g[j + Ncol]*D/(dy*dy);
        if (j == 0) {
          F[i*Ncol] += h[i]*D/(dx*dx);
        } else if (j == Ncol-1) {
          F[i*Ncol + j] += h[i + Nlime]*D/(dx*dx);
        }
      }
      /* First and last col (without first and last elt) */
      else {
        /* First col */
        if (j == 0) {
          F[i*Ncol] += h[i] * D / (dx*dx);
        }
        /* Last col */
        else if(j == Ncol-1) {
          F[i*Ncol + j] += h[i + Nlime]*D/(dx*dx);
        }
      }
    }
  }
}
