#include "Matrix.h"
#include <array>
#include <complex>

class DoubleWilsonMatrix;
class WilsonMatrix;

namespace contrations_with_BK {
	std::complex<double> do_contractions_VVpAA(
	const DoubleWilsonMatrix &amputated_vertex,
    const std::array<DoubleWilsonMatrix, 7> &projector_spin_color_structures);
}
