#pragma once
#include "Vector.h"

void calcul_second_membre(Vector &F, int Nlime, int Ncol, double dx, double dy, double D, Vector const& g, Vector const& h, Vector const& termeSource);