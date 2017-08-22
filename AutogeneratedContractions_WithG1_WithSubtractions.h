#include "Matrix.h"
#include <array>
#include <complex>

class DoubleWilsonMatrix;
class WilsonMatrix;

namespace ContractionsWithG1 {
    Matrix<std::complex<double>, 7> DoContractionsFourQuarkOpFourQuarkExt(
	const std::array<std::array<DoubleWilsonMatrix, 4>, 6> &amputated_diagrams,
	const std::array<DoubleWilsonMatrix, 7> &projector_spin_color_structures);

    std::array<std::complex<double>, 7> DoContractionsFourQuarkOpHFExt(
	const std::array<std::array<DoubleWilsonMatrix, 4>, 6> &amputated_diagrams,
	const DoubleWilsonMatrix &HF_projector_spin_color_structure);

    std::array<std::complex<double>, 7> DoContractionsTwoQuarkOpFourQuarkExt(
	const DoubleWilsonMatrix &amputated_sd_diagram,
	const std::array<DoubleWilsonMatrix, 7> &projector_spin_color_structures);

    std::array<WilsonMatrix, 7> DoContractionsFourQuarkOpTwoQuarkExt(
	const std::array<std::array<WilsonMatrix, 4>, 6> &amputated_twoq_diagrams);

    std::complex<double> DoContractionTwoQuarkOpHFExt(
	const DoubleWilsonMatrix &amputated_sd_diagram,
	const DoubleWilsonMatrix &HF_projector_spin_color_structure);

    Matrix<std::complex<double>, 7> DoQslashContractionsFourQuarkOpFourQuarkExt(
	const std::array<std::array<DoubleWilsonMatrix, 4>, 6> &amputated_diagrams,
	const std::array<DoubleWilsonMatrix, 7> &projector_spin_color_structures);

    std::array<std::complex<double>, 7> DoQslashContractionsTwoQuarkOpFourQuarkExt(
	const DoubleWilsonMatrix &amputated_sd_diagram,
	const std::array<DoubleWilsonMatrix, 7> &projector_spin_color_structures);
}
