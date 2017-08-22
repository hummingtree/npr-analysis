#ifndef A0_H
#define A0_H

#include "JackknifeDatabase.h"
#include "Matrix.h"
#include <complex>

void Test_DoA0();

void DoA0_2p28GeV(const JackknifeDatabase<Matrix<std::complex<double>, 7>> &jack_R_lat_to_RI);

#endif
