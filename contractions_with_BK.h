#include "Matrix.h"
#include "NPR_Utils.h"
#include <array>
#include <complex>

class DoubleWilsonMatrix;
class WilsonMatrix;

namespace contrations_with_BK {
	std::array<std::complex<double>, 5> do_contractions_VVpAA(
	const DoubleWilsonMatrix &amputated_vertex,
    const BK_pscs& pscs);
}
