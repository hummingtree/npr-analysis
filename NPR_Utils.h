#ifndef NPR_UTILS_H
#define NPR_UTILS_H

#include "util.h"
#include "DoubleWilsonMatrix.h"
#include "JackknifeDatabase.h"
#include "WilsonMatrix.h"

#include <array>

std::array<DoubleWilsonMatrix, 7> BuildProjectorSpinColorStructures(Parity parity);

std::array<DoubleWilsonMatrix, 7> BuildQslashProjectorSpinColorStructures(const double q[4], double qsq, Parity parity);
std::array<DoubleWilsonMatrix, 7> BuildDaiqianQslashProjectorSpinColorStructures(const double q[4], double qsq, Parity parity);

void TestQslashProjectors();

DoubleWilsonMatrix BuildHFProjectorSpinColorStructure(Parity parity);

class BK_pscs: public std::array<DoubleWilsonMatrix, 5> {};

BK_pscs build_BK_gammaMu_pscs(Parity parity);
BK_pscs build_BK_Qslash_pscs(const double q[4], double qsq, Parity parity);

template<class T>
std::array<ConfSampleDatabase<T>, 4> MakeLRDiagrams(
    const std::array<ConfSampleDatabase<T>, 4> &lr_diagrams, Parity parity)
{
    std::array<ConfSampleDatabase<T>, 4> ret;

    const ConfSampleDatabase<T> &VV = lr_diagrams[GammaVV];
    const ConfSampleDatabase<T> &AA = lr_diagrams[GammaAA];
    const ConfSampleDatabase<T> &VA = lr_diagrams[GammaVA];
    const ConfSampleDatabase<T> &AV = lr_diagrams[GammaAV];

    switch (parity) {
	case BOTH_PARITIES:
	    ret[GammaLL] = VV + AA - VA - AV;
	    ret[GammaRR] = VV + AA + VA + AV;
	    ret[GammaLR] = VV - AA + VA - AV;
	    ret[GammaRL] = VV - AA - VA + AV;
	    break;

	case POSITIVE_PARITY:
	    ret[GammaLL] = VV + AA;
	    ret[GammaRR] = VV + AA;
	    ret[GammaLR] = VV - AA;
	    ret[GammaRL] = VV - AA;
	    break;

	case NEGATIVE_PARITY:
	    ret[GammaLL] = -VA - AV;
	    ret[GammaRR] = VA + AV;
	    ret[GammaLR] = VA - AV;
	    ret[GammaRL] = -VA + AV;
	    break;
    }

    return ret;
}

template<class T>
ConfSampleDatabase<T> MakeLeftHandedDiagrams(
    const std::array<ConfSampleDatabase<T>, 2> &lr_diagrams, Parity parity)
{
    const ConfSampleDatabase<T> &V = lr_diagrams[0];
    const ConfSampleDatabase<T> &A = lr_diagrams[1];

    switch (parity) {
	case BOTH_PARITIES: return V - A;
	case POSITIVE_PARITY: return V;
	case NEGATIVE_PARITY: default: return (-1.0) * A;
    }
}



// Computes subtraction coefficients for Nsub two-quark operators.
// The inputs are 
// - the amputated two-quark diagram of the unsubtracted operator
// - the amputated two-quark diagrams of the subtraction operators
// - the spin matrices onto which to project.
// The goal is to construct a subtracted operator whose amputated two-quark
// diagram has zero projection onto all of the given projectors
template<std::size_t Nsub>
std::array<std::complex<double>, Nsub> ComputeSubtractionCoefficients(
    const WilsonMatrix &amputated_twoq_diagram,
    const std::array<WilsonMatrix, Nsub> &amputated_twoq_subtraction_diagrams,
    const std::array<SpinMatrix, Nsub> &projectors)
{
    // goal: Tr ( projectors[i] * (amputated_twoq_diagram - sum_j (B_j * amputated_twoq_subtraction_diagrams[j])) ) = 0
    // for all i. This can be written
    // A_ij B_j = c_i
    // where
    // A_ij = Tr (projectors[i] * amputated_twoq_subtraction_diagrams[j])
    // c_i = Tr (projectors[i] * amputated_twoq_diagram)

    Matrix<std::complex<double>, Nsub> A;
    std::array<std::complex<double>, Nsub> c;

    for (unsigned i = 0; i < Nsub; ++i) {
	c[i] = (projectors[i] * amputated_twoq_diagram).Trace();
	for (unsigned j = 0; j < Nsub; ++j) {
	    A[i][j] = (projectors[i] * amputated_twoq_subtraction_diagrams[j]).Trace();
	}
    }

    std::array<std::complex<double>, Nsub> B;
    B = A.Inverse() * c;
    return B;
}

template<std::size_t Nsub>
WilsonMatrix ComputeSubtractedTwoQuarkGreensFunction(const WilsonMatrix &amputated_twoq_diagram,
    const std::array<WilsonMatrix, Nsub> &amputated_twoq_subtraction_diagrams,
    const std::array<std::complex<double>, Nsub> &subtraction_coefficients)
{
    WilsonMatrix ret = amputated_twoq_diagram;
    for (unsigned i = 0; i < Nsub; ++i) {
	ret -= subtraction_coefficients[i] * amputated_twoq_subtraction_diagrams[i];
    }
    return ret;
}

template<int N>
void PrintRealMatrixWithError(const char *name,
    const Matrix<double, N> &vals,
    const Matrix<double, N> &errs)
{
    printf("%s:\n", name);
    for (int row = 0; row < N; ++row) {
	for (int col = 0; col < N; ++col) {
	    double val = vals[row][col];
	    double err = errs[row][col];
	    if (val == 0 && err == 0) {
		printf("%15s", "0 ");
	    } else if (std::abs(val) < 0.0000001) {
		printf("%15s", "~0 ");
	    } else {
		char val_with_err[512];
		SprintfValueAndError(val_with_err, val, err);
		printf("%15s", val_with_err);
	    }
	}
	printf("\n");
    }
}

template<int N>
void PrintRealPartMatrixWithError(const char *name, const JackknifeDatabase<Matrix<std::complex<double>, N> > &jack_db)
{
    Matrix<double, N> vals = RealMatrix(jack_db.CentralValue());
    Matrix<double, N> errs = RealMatrix(jack_db.Error());
    PrintRealMatrixWithError(name, vals, errs);
}

template<int N>
void PrintImaginaryPartMatrixWithError(const char *name, const JackknifeDatabase<Matrix<std::complex<double>, N> > &jack_db)
{
    Matrix<double, N> vals = ImaginaryMatrix(jack_db.CentralValue());
    Matrix<double, N> errs = ImaginaryMatrix(jack_db.Error());
    PrintRealMatrixWithError(name, vals, errs);
}


template<int N>
void LatexPrintRealMatrixWithErrorEnforcingRepresentations(const char* name,
    const Matrix<double, N> &vals, const Matrix<double, N> &errs)
{
    printf("LATEX VERSION!!!!! of %s:\n", name);
    assert(N == 7 || N == 8);
    std::array<int, 8> reps{{ 1, 2, 2, 2, 2, 3, 3, 2 }}; // which SU(3)_L(x)SU(3)_R representation each row/column corresponds to

    for (int row = 0; row < N; ++row) {
	for (int col = 0; col < N; ++col) {
	    double val = vals[row][col];
	    double err = errs[row][col];
	    if (reps[row] == reps[col]) {
		if (std::abs(val) < 0.0000001) {
		    printf("%15s", "0");
		} else {
		    char val_with_err[512];
		    SprintfValueAndError(val_with_err, val, err);
		    printf("%15s", val_with_err);
		}
	    } else {
		if (std::abs(val) >= 0.0000001) {
		    printf("%15s", "ASDF");
		} else {
		    printf("%15s", "");
		}
	    }
	    if (col < N - 1) {
		printf(" &");
	    }
	}
	printf(" \\\\\n");
    }
}

template<int N>
void LatexPrintRealPartMatrixWithErrorEnforcingRepresentations(const char *name, 
    const JackknifeDatabase<Matrix<std::complex<double>, N> > &jack_db)
{
    Matrix<double, N> vals = RealMatrix(jack_db.CentralValue());
    Matrix<double, N> errs = RealMatrix(jack_db.Error());
    LatexPrintRealMatrixWithErrorEnforcingRepresentations(name, vals, errs);
}

template<class T>
void Zero(ConfSampleDatabase<T> &db)
{
    for (T& t : db.samples) t.Zero();
}

#endif
