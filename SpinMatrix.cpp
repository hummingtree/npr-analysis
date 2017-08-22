#include "SpinMatrix.h"

SpinMatrix::SpinMatrix()
{
    for (int s1 = 0; s1 < 4; ++s1) {
	for (int s2 = 0; s2 < 4; ++s2) {
	    data[s1][s2] = 0.0;
	}
    }
}

std::complex<double>& SpinMatrix::operator()(int s1, int s2)
{
    return this->data[s1][s2];
}

const std::complex<double>& SpinMatrix::operator()(int s1, int s2) const
{
    return this->data[s1][s2];
}

SpinMatrix SpinMatrix::operator*(const SpinMatrix &rhs) const
{
    SpinMatrix ret;
    for (int s1 = 0; s1 < 4; ++s1) {
	for (int s2 = 0; s2 < 4; ++s2) {
	    for (int x = 0; x < 4; ++x) {
		ret.data[s1][s2] += this->data[s1][x] * rhs.data[x][s2];
	    }
	}
    }
    return ret;
}

SpinMatrix& SpinMatrix::operator*=(std::complex<double> coeff)
{
    for (int s1 = 0; s1 < 4; ++s1) {
	for (int s2 = 0; s2 < 4; ++s2) {
	    this->data[s1][s2] *= coeff;
	}
    }
    return *this;
}

SpinMatrix operator*(std::complex<double> coeff, const SpinMatrix &rhs)
{
    SpinMatrix ret = rhs;
    ret *= coeff;
    return ret;
}

SpinMatrix SpinMatrix::operator+(const SpinMatrix &rhs) const
{
    SpinMatrix ret;
    for (int s1 = 0; s1 < 4; ++s1) {
	for (int s2 = 0; s2 < 4; ++s2) {
	    ret.data[s1][s2] = this->data[s1][s2] + rhs.data[s1][s2];
	}
    }
    return ret;
}

SpinMatrix SpinMatrix::operator-(const SpinMatrix &rhs) const
{
    SpinMatrix ret;
    for (int s1 = 0; s1 < 4; ++s1) {
	for (int s2 = 0; s2 < 4; ++s2) {
	    ret.data[s1][s2] = this->data[s1][s2] - rhs.data[s1][s2];
	}
    }
    return ret;
}

SpinMatrix SpinMatrix::One()
{
    SpinMatrix ret;
    for (int i = 0; i < 4; ++i) ret.data[i][i] = 1.0;
    return ret;
}

SpinMatrix SpinMatrix::Gamma5()
{
    SpinMatrix ret;
    ret.data[0][0] = ret.data[1][1] = 1.0;
    ret.data[2][2] = ret.data[3][3] = -1.0;
    return ret;
}

SpinMatrix SpinMatrix::Gamma(int mu)
{
    SpinMatrix ret;

    switch (mu) {
	case 0:
	    ret.data[0][3] = ret.data[1][2] = std::complex<double>(0, 1.0);
	    ret.data[2][1] = ret.data[3][0] = std::complex<double>(0, -1.0);
	    break;

	case 1:
	    ret.data[0][3] = ret.data[3][0] = -1.0;
	    ret.data[1][2] = ret.data[2][1] = 1.0;
	    break;

	case 2:
	    ret.data[0][2] = ret.data[3][1] = std::complex<double>(0, 1.0);
	    ret.data[1][3] = ret.data[2][0] = std::complex<double>(0, -1.0);
	    break;

	case 3:
	    ret.data[0][2] = ret.data[1][3] = ret.data[2][0] = ret.data[3][1] = 1.0;
	    break;

	default:
	    printf("unknown mu = %d\n", mu);
	    exit(-1);
    }

    return ret;
}

SpinMatrix SpinMatrix::GammaMuGamma5(int mu)
{
    return Gamma(mu) * Gamma5();
}

SpinMatrix SpinMatrix::Sigma(int mu, int nu)
{
    return 0.5 * (Gamma(mu) * Gamma(nu) - Gamma(nu) * Gamma(mu));
}

SpinMatrix SpinMatrix::SigmaGamma5(int mu, int nu)
{
    return Sigma(mu, nu) * Gamma5();
}

SpinMatrix SpinMatrix::Slash(const double p[4])
{
    SpinMatrix ret;
    for (int mu = 0; mu < 4; ++mu) ret = ret + p[mu] * Gamma(mu);
    return ret;
}