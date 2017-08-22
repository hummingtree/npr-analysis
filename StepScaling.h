#ifndef STEP_SCALING_H
#define STEP_SCALING_H

#include "JackknifeDatabase.h"
#include "Matrix.h"

// Can step scale 7x7 or 8x8 Z-factors
template<int N>
JackknifeDatabase<Matrix<std::complex<double>, N>> StepScaleJackknifeZfactors(
    const JackknifeDatabase<Matrix<std::complex<double>, N>> &jack_lat_coarse_to_RI_low,
    const JackknifeDatabase<Matrix<std::complex<double>, N>> &jack_lat_fine_to_RI_low,
    const JackknifeDatabase<Matrix<std::complex<double>, N>> &jack_lat_fine_to_RI_high)
{
    printf("==== STEP SCALING of %dx%d Z-FACTORS ====\n", N, N);

    assert(jack_lat_fine_to_RI_low.Size() == jack_lat_fine_to_RI_high.Size());

    Matrix<std::complex<double>, N> central_lat_coarse_to_RI_low = jack_lat_coarse_to_RI_low.CentralValue();
    Matrix<std::complex<double>, N> central_lat_fine_to_RI_low = jack_lat_fine_to_RI_low.CentralValue();
    Matrix<std::complex<double>, N> central_lat_fine_to_RI_high = jack_lat_fine_to_RI_high.CentralValue();

    // build jackknife database for RI_low -> RI_high using fine data
    JackknifeDatabase<Matrix<std::complex<double>, N>> jack_low_to_high;
    Matrix<std::complex<double>, N> central_low_to_high = central_lat_fine_to_RI_high * central_lat_fine_to_RI_low.Inverse();
    jack_low_to_high.Add(central_low_to_high);
    for (int i = 1; i < jack_lat_fine_to_RI_low.Size(); ++i) {
	Matrix<std::complex<double>, N> sample = jack_lat_fine_to_RI_high[i] * jack_lat_fine_to_RI_low[i].Inverse();
	jack_low_to_high.Add(sample);
    }
    assert(jack_low_to_high.Size() == jack_lat_fine_to_RI_low.Size());

    // a superjackknife database of the full coarse lattice operator -> high scale RI operator
    // Z-factor matrix, which is the result of the step scaling
    JackknifeDatabase<Matrix<std::complex<double>, N>> jack_lat_coarse_to_RI_high;
    Matrix<std::complex<double>, N> central_lat_coarse_to_RI_high = central_low_to_high * central_lat_coarse_to_RI_low;
    jack_lat_coarse_to_RI_high.Add(central_lat_coarse_to_RI_high);
    // portion of the superjackknife where the coarse lattice varies
    for (int i = 1; i < jack_lat_coarse_to_RI_low.Size(); ++i) { // start at i=1 to skip central
	Matrix<std::complex<double>, N> sample = central_low_to_high * jack_lat_coarse_to_RI_low[i];
	jack_lat_coarse_to_RI_high.Add(sample);
    }
    // portion of the superjackknife where the fine lattice varies
    for (int i = 1; i < jack_low_to_high.Size(); ++i) { // start at i=1 to skip central
	Matrix<std::complex<double>, N> sample = jack_low_to_high[i] * central_lat_coarse_to_RI_low;
	jack_lat_coarse_to_RI_high.Add(sample);
    }
    
    printf("\n");
    printf("Nbins_coarse = %d\n", jack_lat_coarse_to_RI_low.Size() - 1);
    printf("Nbins_fine = %d\n", jack_low_to_high.Size() - 1);
    printf("\n");

    PrintRealPartMatrixWithError("real low_to_high", jack_low_to_high);
    LatexPrintRealPartMatrixWithErrorEnforcingRepresentations("real low_to_high", jack_low_to_high);
    PrintImaginaryPartMatrixWithError("imaginary low_to_high", jack_low_to_high);

    printf("\n");

    PrintRealPartMatrixWithError("real lat_coarse_to_RI_high", jack_lat_coarse_to_RI_high);
    LatexPrintRealPartMatrixWithErrorEnforcingRepresentations("real lat_coarse_to_RI_high", jack_lat_coarse_to_RI_high);
    PrintImaginaryPartMatrixWithError("imaginary lat_coarse_to_RI_high", jack_lat_coarse_to_RI_high);
    
    printf("==== DONE STEP SCALING of %dx%d Z-FACTORS ====\n", N, N);

    return jack_lat_coarse_to_RI_high;
}

#endif
