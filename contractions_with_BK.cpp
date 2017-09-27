#include "DoubleWilsonMatrix.h"
#include "WilsonMatrix.h"
#include "Matrix.h"
#include "util.h"
#include <array>
#include <complex>
#include "NPR_Utils.h"

namespace contrations_with_BK {
// 
std::array<std::complex<double>, 5> do_contractions_VVpAA(
	const DoubleWilsonMatrix &amputated_vertex,
    const BK_pscs& pscs)
{
	std::array<std::complex<double>, 5> rtn;	
// Now only the GammaMu scheme.
	DoubleWilsonMatrix Pi_trtr = amputated_vertex; // The trace-trace structure. (4.104) in C. Kelly's thesis.
	DoubleWilsonMatrix Pi_tr = amputated_vertex; // The trace structure.
	// Pi_tr.Swap01();
	Pi_tr.Swap13();
	// P_gammaMu.Swap01();
	// P_gammaMu.Swap23();
	for(int i = 0; i < 5; i++){
		DoubleWilsonMatrix P_gammaMu = pscs[i]; // using the POSITIVE parity
		rtn[i] = 2.*(Pi_trtr - Pi_tr).Project(P_gammaMu);
	}
	return rtn;
}

Matrix<std::complex<double>, 5> do_contractions_DSeq2(
	const std::array<DoubleWilsonMatrix, 5> &amputated_vertex,
    const BK_pscs& pscs)
{
	Matrix<std::complex<double>, 5> rtn;	
// Now only the GammaMu scheme.
	for(int i = 0; i < 5; i++){
		DoubleWilsonMatrix Pi_trtr = amputated_vertex[i]; // The trace-trace structure. (4.104) in C. Kelly's thesis.
		DoubleWilsonMatrix Pi_tr = amputated_vertex[i]; // The trace structure.
		// Pi_tr.Swap01();
		Pi_tr.Swap13();
		// P_gammaMu.Swap01();
		// P_gammaMu.Swap23();
		for(int j = 0; j < 5; j++){
			DoubleWilsonMatrix P_gammaMu = pscs[j]; // using the POSITIVE parity
			rtn[i][j] = 2.*(Pi_trtr - Pi_tr).Project(P_gammaMu);
		}
	
	}
	return rtn;
}

} // end namespace

