#pragma once
#include "math/Vector.h"

void g(int me, int Ncol, double dx, double Ly, Vector& gme, int mode);

void h(int me, int Nlime, int i1, double dy, double Lx, Vector& hme, int mode);

void fsource(int me, int Ncol, int Nlime, int i1, double dx, double dy, double Lx, double Ly, double t, Vector& fsourceme, int mode);

void solution(int me, int Ncol, int Nlime, int i1, double dx, double dy, double Lx, double Ly, double t, Vector& solutionme, int mode);

void charge(int rowCount, int np, int rank /* aka "me" */, int& iBegin, int& iEnd);