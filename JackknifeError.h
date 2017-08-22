#ifndef JACKKNIFE_ERROR_H
#define JACKKNIFE_ERROR_H

#include <vector>
#include <complex>
#include <array>
#include <cassert>

#include "Matrix.h"

double JackknifeError(const std::vector<double> &jackknife_values);

std::complex<double> JackknifeError(const std::vector<std::complex<double> > &jackknife_values);

class SpinMatrix;
SpinMatrix JackknifeError(const std::vector<SpinMatrix> &jackknife_values);


template<class T>
std::vector<T> JackknifeError(const std::vector<std::vector<T> > &jackknife_values)
{
    int Njack = jackknife_values.size();
    int Nvals = jackknife_values[0].size();
    std::vector<T> errs(Nvals);

    for (int val = 0; val < Nvals; ++val) {
	std::vector<T> jack_vals(Njack);
	for (int jack = 0; jack < Njack; ++jack) {
	    jack_vals[jack] = jackknife_values[jack][val];
	}
	errs[val] = JackknifeError(jack_vals);
    }

    return errs;
}




template<class T, int N>
std::vector<T> ListFromMatrix(const Matrix<T, N> &mat)
{
    std::vector<T> ret(N*N);
    for (int row = 0; row < N; ++row) {
	for (int col = 0; col < N; ++col) {
	    ret[row*N + col] = mat[row][col];
	}
    }
    return ret;
}

template<class T, int N>
Matrix<T, N> MatrixFromList(const std::vector<T> &vec)
{
    Matrix<T, N> ret;
    for (int row = 0; row < N; ++row) {
	for (int col = 0; col < N; ++col) {
	    ret[row][col] = vec[row*N + col];
	}
    }
    return ret;
}

template<class T, int N>
Matrix<T, N> JackknifeError(const std::vector<Matrix<T, N> > &jackknife_mats)
{
    int Njack = jackknife_mats.size();
    std::vector<std::vector<T> > flat_jackknife_mats(Njack);

    for (int jack = 0; jack < Njack; ++jack) {
	flat_jackknife_mats[jack] = ListFromMatrix<T, N>(jackknife_mats[jack]);
    }

    std::vector<T> flat_errs = JackknifeError(flat_jackknife_mats);

    Matrix<T, N> errs = MatrixFromList<T, N>(flat_errs);
    return errs;
}

template<class T, std::size_t N>
std::vector<T> ListFromArray(const std::array<T, N> &arr)
{
    return std::vector<T>(arr.begin(), arr.end());
}

template<class T, std::size_t N>
std::array<T, N> ArrayFromList(const std::vector<T> &vec)
{
    assert(vec.size() == N);
    std::array<T, N> ret;
    std::copy(vec.begin(), vec.end(), ret.begin());
    return ret;
}

template<class T, std::size_t N>
std::array<T, N> JackknifeError(const std::vector<std::array<T, N> > &jackknife_arrs)
{
    int Njack = jackknife_arrs.size();
    std::vector<std::vector<T> > jackknife_vecs(Njack);

    for (int jack = 0; jack < Njack; ++jack) {
	jackknife_vecs[jack] = ListFromArray<T, N>(jackknife_arrs[jack]);
    }

    std::vector<T> vec_errs = JackknifeError(jackknife_vecs);

    std::array<T, N> errs = ArrayFromList<T, N>(vec_errs);
    return errs;
}

#endif
