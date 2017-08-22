#ifndef NEW_NPR_SIMPLEG1_WITHSUBTRACTIONS_H
#define NEW_NPR_SIMPLEG1_WITHSUBTRACTIONS_H

#include "JackknifeDatabase.h"
#include "Matrix.h"
#include "NPRSettings.h"

struct NPRResults
{
    JackknifeDatabase<Matrix<std::complex<double>, 7>> jack_Z_7x7;
    JackknifeDatabase<Matrix<std::complex<double>, 8>> jack_Z_8x8;
    JackknifeDatabase<Matrix<std::complex<double>, 7>> jack_R;
};

JackknifeDatabase<Matrix<std::complex<double>, 7>> ComputeRMatrixJackknife(
        const JackknifeDatabase<Matrix<std::complex<double>, 8>> &Z_8x8,
        const JackknifeDatabase<Matrix<std::complex<double>, 7>> &Z_7x7);

NPRResults NPR_SimpleG1_WithSubtractions(NPRSettings &sett);

NPRResults NPR_SimpleG1_WithSubtractions_TwoMoms(NPRSettings &sett1, NPRSettings &sett2);

NPRResults NPR_SimpleG1_WithSubtractions_TwoMomsForG1(NPRSettings &sett1, NPRSettings &sett2);

#endif
