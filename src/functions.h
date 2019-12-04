#pragma once
#include "math/Mat.h"

Mat fsta(double Lx, double Ly, int Nx, int Ny);

Mat gsta(double Lx, double Ly, int Nx, int Ny);

Mat hsta(double Lx, double Ly, int Nx, int Ny);

Mat fper(double Lx, double Ly, int Nx, int Ny);

Mat gper(double Lx, double Ly, int Nx, int Ny);

Mat hper(double Lx, double Ly, int Nx, int Ny);

void charge(int rowCount, int np, int rank /* aka "me" */, int& iBegin, int& iEnd);
