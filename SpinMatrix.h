#ifndef SPIN_MATRIX_H
#define SPIN_MATRIX_H

#include "Matrix.h"
#include <complex>

class SpinMatrix
{
public:
    SpinMatrix();

    Matrix<std::complex<double>, 4> data;

    std::complex<double>& operator()(int s1, int s2);
    const std::complex<double>& operator()(int s1, int s2) const;

    SpinMatrix operator*(const SpinMatrix &rhs) const;
    SpinMatrix& operator*=(std::complex<double> coeff);
    SpinMatrix operator+(const SpinMatrix &rhs) const;
    SpinMatrix operator-(const SpinMatrix &rhs) const;

    static SpinMatrix One();
    static SpinMatrix Gamma5();
    static SpinMatrix Gamma(int mu);
    static SpinMatrix GammaMuGamma5(int mu);
    static SpinMatrix Sigma(int mu, int nu);
    static SpinMatrix SigmaGamma5(int mu, int nu);

    static SpinMatrix Slash(const double p[4]);
};

SpinMatrix operator*(std::complex<double> coeff, const SpinMatrix &sm);



#endif