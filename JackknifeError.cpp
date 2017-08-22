#include "JackknifeError.h"

#include "SpinMatrix.h"

double JackknifeError(const std::vector<double> &jackknife_values)
{
    if (jackknife_values.size() == 0) return 0.0;

    int N = jackknife_values.size() - 1;
    double central_value = jackknife_values[0];

    double sum_sq = 0.0;
    for (unsigned i = 1; i < jackknife_values.size(); ++i) {
	double diff = jackknife_values[i] - central_value;
	sum_sq += diff * diff;
    }
    return std::sqrt(N * sum_sq / (N - 1));
}

std::complex<double> JackknifeError(const std::vector<std::complex<double> > &jackknife_values)
{
    int Njack = jackknife_values.size();
    std::vector<double> real_values(Njack);
    std::vector<double> imag_values(Njack);
    for (int i = 0; i < Njack; ++i) {
	real_values[i] = jackknife_values[i].real();
	imag_values[i] = jackknife_values[i].imag();
    }
    std::complex<double> err = std::complex<double>(JackknifeError(real_values), JackknifeError(imag_values));
    return err;
}

SpinMatrix JackknifeError(const std::vector<SpinMatrix> &jackknife_values)
{
    std::vector<Matrix<std::complex<double>, 4>> mat_values;
    for (unsigned i = 0; i < jackknife_values.size(); ++i) {
	mat_values.push_back(jackknife_values[i].data);
    }
    SpinMatrix ret;
    ret.data = JackknifeError(mat_values);
    return ret;
}
