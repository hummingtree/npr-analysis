#ifndef G1CHECK_H
#define G1CHECK_H

#include "DoubleWilsonMatrix.h"
#include "util.h"
#include <array>
#include <complex>

DoubleWilsonMatrix BuildLowestOrderG1TensorStructure(const double q[4], double qsq, Parity parity);

class NPRSettings;

std::array<std::complex<double>, 8> ComputeTreeLevelG1ProjectedGreensFuncs(const NPRSettings &sett);

void CheckG1TreeLevelProjectedGreensFunctions(const NPRSettings &sett);

#endif
