#include "functions.h"

#include <cmath>
#include <cstddef>
#include <mpi.h>
#include <cstring>

void g(int me, int Ncol, double dx, double Ly, Vector& gme, int mode) {
    int i;
    double x;

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    switch (mode)
    {
    case 1:
        gme.set_value(0.);
        break;
    case 2:
        gme.set_value(0.);
        if (me == 0 || me == size) {
            /* Compute shift in gme plus constant value */
            int bool_first = (me == 0);
            int shift = (!bool_first) * Ncol;
            double cos_border = (bool_first) ? 1 : cos(Ly);
            for (i = 0; i < Ncol; ++i) {
                gme[shift + i] = sin((double)(i) * dx) + cos_border;
            }

        }
        break;
    case 3:
        gme.set_value(0.);
        break;
    default:
        break;
    }
}

void h(int me, int Nlime, int i1, double dy, double Lx, Vector& hme, int mode) {
    int i;
    double y;

    switch (mode)
    {
    case 1:
        hme.set_value(0.);
        break;
    case 2:
        y = (double)(i - i1 + 1)*dy;
        for (i = 0; i < Nlime; ++i) {
            hme[Nlime + i] = cos(y) + sin(Lx);
        }
        for (i = 0; i < Nlime; ++i) {
            hme[i]         = cos(y);
        }
        break;
    case 3:
        hme.set_value(1.);
        break;
    default:
        break;
    }
}



void fsource(int me, int Ncol, int Nlime, int i1, double dx, double dy, double Lx, double Ly, double t, Vector& fsourceme, int mode) {
    switch (mode)
    {
    case 1:
        for (int i = 0; i < Nlime; ++i) {
            for (int j = 0; j < Ncol; ++j) {
                double x = (double)(i + i1 - 1) * dx;
                double y = (double)(j) * dy;
                fsourceme[i * Nlime + j] = 2 * (y - y**2 + x - x**2);
            }
        }
        break;
    case 2:
        for (int i = 0; i < Nlime; ++i) {
            for (int j = 0; j < Ncol; ++j) {
                double x = (double)(i + i1 - 1) * dx;
                double y = (double)(j) * dy;
                fsourceme[i * Nlime + j] = sin(x) + cos(y);
            }
        }
        
        break;

    case 3:
        for (int i = 0; i < Nlime; ++i) {
            for (int j = 0; j < Ncol; ++j) {
                double x = (double)(i + i1 - 1) * dx;
                double y = (double)(j) * dy;
                double tmpx = (x - Lx*0.5);
                double tmpy = (y - Ly*0.5);
                fsourceme[i * Nlime + j] = exp( - tmpx * tmpx  - tmpy * tmpy) * cos (t * M_PI * 0.5);
            }
        }

        break;

    default:
        break;
    }
}


void charge( int rowCount, int np, int rank /* aka "me" */, int &iBegin, int &iEnd )
{
    int loadSize = rowCount / np;
    int leftover = rowCount % np;

    if ( leftover == 0 ) {
        iBegin = rank * loadSize;
        iEnd   = iBegin + loadSize - 1;
    }
    else {
        if ( rank < leftover ) {
            iBegin = rank * ( loadSize + 1 );
            iEnd   = iBegin + loadSize;
        }
        else {
            iBegin = leftover * ( loadSize + 1 ) + ( rank - leftover ) * loadSize;
            iEnd   = iBegin + loadSize - 1;
        }
    }
}

