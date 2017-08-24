#include "DoubleWilsonMatrix.h"
#include "WilsonMatrix.h"
#include "Matrix.h"
#include "util.h"
#include <array>
#include <complex>



namespace contrations_with_BK {
// 
std::complex<double> do_contractions_VVpAA(
	const DoubleWilsonMatrix &amputated_vertex,
    const std::array<DoubleWilsonMatrix, 7> &projector_spin_color_structures)
{
// Now only the GammaMu scheme.
	DoubleWilsonMatrix P_gammaMu = projector_spin_color_structures[0];// using the POSITIVE parity
	DoubleWilsonMatrix Pi_trtr = amputated_vertex; // The trace-trace structure. (4.104) in C. Kelly's thesis.
	DoubleWilsonMatrix Pi_tr = amputated_vertex; // The trace structure.
	Pi_tr.Swap23();
	Pi_tr.Swap13();
	// P_gammaMu.Swap01();
	// P_gammaMu.Swap23();
	return (Pi_trtr - Pi_tr).Project(P_gammaMu);
}


} // end namespace

