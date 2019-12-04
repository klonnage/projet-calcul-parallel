#include "functions.h"

#include <cmath>
#include <cstddef>

Mat fsta( double Lx, double Ly, int Nx, int Ny )
{
    double dx = Lx / static_cast<double>( Nx );
    double dy = Ly / static_cast<double>( Ny );
    Mat    result( Nx, Ny, 0. );

    for ( size_t i = 0; i < Nx; ++i ) {
        for ( size_t j = 0; j < Ny; ++j ) {
            double x       = static_cast<double>( i ) * dx;
            double y       = static_cast<double>( j ) * dy;
            result( i, j ) = 2 * ( x - x * x + y - y * y );
        }
    }

    return result;
}

Mat gsta( double Lx, double Ly, int Nx, int Ny )
{
    return Mat( Nx, Ny, 0. );
}

Mat hsta( double Lx, double Ly, int Nx, int Ny )
{
    return Mat( Nx, Ny, 0. );
}

Mat fper( double Lx, double Ly, int Nx, int Ny )
{
    Mat    ret( Nx, Ny, 0.0 );
    double dx = Lx / Nx, dy = Ly / Ny;

    for ( int i = 0; i < Nx; ++i ) {
        for ( int j = 0; j < Ny; ++j ) {
            double x    = static_cast<double>( i ) * dx;
            double y    = static_cast<double>( j ) * dy;
            ret( i, j ) = std::sin( x ) + std::cos( y );
        }
    }
}

Mat gper( double Lx, double Ly, int Nx, int Ny ) {}

Mat hper( double Lx, double Ly, int Nx, int Ny ) {}

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


