#ifndef MATRIX_H
#define MATRIX_H

#include <array>
#include <iostream>
#include <complex>
#include <iomanip>
#include <cstdio>
#include <cstdlib>

template<class T, int N>
class Row
{
public:
    std::array<T, N> data;

    T& operator[](int i) { return data[i]; }
    const T& operator[](int i) const { return data[i]; }

    Row<T,N> operator*(T coeff) const
    {
	Row<T, N> ret;
	for (int i = 0; i < N; ++i) ret[i] = (*this)[i] * coeff;
	return ret;
    }

    Row<T, N> operator+(Row<T, N> rhs) const
    {
	Row<T, N> ret;
	for (int i = 0; i < N; ++i) ret[i] = (*this)[i] + rhs[i];
	return ret;
    }

    Row<T, N>& operator=(const std::array<T, N> &rhs)
    {
	this->data = rhs;
	return *this;
    }

    void Print() const;
};

template<class T> 
void PrintMatrixEntry(T x)
{
    std::cout << std::setw(30);
    if (std::abs(x) < 1e-10) std::cout << 0;
    else std::cout << x;
}

template<> 
inline void PrintMatrixEntry<double>(double x)
{
    std::cout << std::setw(15);
    if (std::abs(x) < 1e-10) std::cout << 0;
    else std::cout << x;
}

template<class T, int N>
void Row<T,N>::Print() const
{
    for (int i = 0; i < N; ++i) {
	PrintMatrixEntry(data[i]);
    }
    std::cout << "\n";
}

template<class T, int N>
class Matrix
{
public:
    std::array<Row<T, N>, N> data;

    Row<T, N>& operator[](int i) { return data[i]; }
    const Row<T,N>& operator[](int i) const { return data[i]; }

    void Zero();
    void Identity();

    Matrix<T, N> operator*(const Matrix<T, N> &rhs) const;
    std::array<T, N> operator*(const std::array<T, N> &rhs) const;
    Matrix<T, N> operator+(const Matrix<T, N> &rhs) const;
    Matrix<T, N> operator-(const Matrix<T, N> &rhs) const;

    Matrix<T, N>& operator+=(const Matrix<T, N> &rhs);

    Matrix<T, N>& operator*=(std::complex<double> rhs);

    Matrix<T, N> Inverse() const;
    Matrix<T, N> Dagger() const;
    T Trace() const;

    void Print(const char *name) const;
};

template<class T, int N>
void Matrix<T, N>::Zero()
{
    for (int i = 0; i < N; ++i) {
	for (int j = 0; j < N; ++j) {
	    this->data[i][j] = 0.0;
	}
    }
}

template<class T, int N>
void Matrix<T, N>::Identity()
{
    for (int i = 0; i < N; ++i) {
	for (int j = 0; j < N; ++j) {
	    this->data[i][j] = (i == j ? 1.0 : 0.0);
	}
    }
}

template<class T, int N>
Matrix<T, N> Matrix<T, N>::operator*(const Matrix<T, N> &rhs) const
{
    Matrix<T, N> ret;
    ret.Zero();
    for (int i = 0; i < N; ++i) {
	for (int j = 0; j < N; ++j) {
	    for (int k = 0; k < N; ++k) {
		ret.data[i][j] += this->data[i][k] * rhs.data[k][j];
	    }
	}
    }
    return ret;
}

template<class T, int N>
std::array<T, N> Matrix<T, N>::operator*(const std::array<T, N> &rhs) const
{
    std::array<T, N> ret;
    for (int i = 0; i < N; ++i) {
	ret[i] = 0;
	for (int j = 0; j < N; ++j) {
	    ret[i] += this->data[i][j] * rhs[j];
	}
    }
    return ret;
}

template<class T, int N>
Matrix<T, N>& Matrix<T, N>::operator+=(const Matrix<T, N> &rhs)
{
    for (int row = 0; row < N; ++row) {
	for (int col = 0; col < N; ++col) {
	    this->data[row][col] += rhs.data[row][col];
	}
    }
    return *this;
}

template<class T, int N>
Matrix<T, N>& Matrix<T, N>::operator*=(std::complex<double> rhs)
{
    for (int row = 0; row < N; ++row) {
	for (int col = 0; col < N; ++col) {
	    this->data[row][col] *= rhs;
	}
    }
    return *this;
}

template<class T, int N>
Matrix<T, N> Matrix<T, N>::operator+(const Matrix<T, N> &rhs) const
{
    Matrix<T, N> ret;
    for (int i = 0; i < N; ++i) {
	for (int j = 0; j < N; ++j) {
	    ret.data[i][j] = this->data[i][j] + rhs.data[i][j];
	}
    }
    return ret;
}

template<class T, int N>
Matrix<T, N> Matrix<T, N>::operator-(const Matrix<T, N> &rhs) const
{
    Matrix<T, N> ret;
    for (int i = 0; i < N; ++i) {
	for (int j = 0; j < N; ++j) {
	    ret.data[i][j] = this->data[i][j] - rhs.data[i][j];
	}
    }
    return ret;
}

template<class T, int N>
Matrix<T, N> Matrix<T, N>::Inverse() const
{
    Matrix<T, N> tmp = *this;
    Matrix<T, N> ret;
    ret.Identity();

    const double eps = 1e-10;

    // eliminate entries below diagonal
    for (int column = 0; column < N; ++column) {
	// first set (column,column) = 1
	for (int row = column; row < N; ++row) {
	    if (std::abs(tmp[row][column]) > eps) {
		T scale = (1.0 - tmp[column][column]) / tmp[row][column];
		tmp[column] = tmp[column] + tmp[row] * scale;
		ret[column] = ret[column] + ret[row] * scale;
		break;
	    }
	    if (row == N - 1) {
		printf("failure A\n");
		exit(-1);
	    }
	}
	// Now tmp[column][column] should be 1

	// next eliminate all zeros below the diagonal in column
	for (int row = column + 1; row < N; ++row) {
	    T scale = -tmp[row][column] / tmp[column][column];
	    tmp[row] = tmp[row] + tmp[column] * scale;
	    ret[row] = ret[row] + ret[column] * scale;
	}
	// Now column should be all zeros below the diagonal
    }
    // Now tmp should be all zeros below the diagonal

    // eliminate entries above the diagonal
    for (int column = N - 1; column > 0; --column) {
	for (int row = column - 1; row >= 0; --row) {
	    T scale = -tmp[row][column] / tmp[column][column];
	    tmp[row] = tmp[row] + tmp[column] * scale;
	    ret[row] = ret[row] + ret[column] * scale;
	}
    }

    // Now tmp should be the identity, and ret should be the inverse matrix
    return ret;
}

template<class T, int N>
void Matrix<T, N>::Print(const char *name) const
{
    printf("%s:\n", name);
    for (int i = 0; i < N; ++i) data[i].Print();
}

template<class T, int N>
Matrix<T, N> RealMatrix(const Matrix<std::complex<T>, N> &complex_mat)
{
    Matrix<T, N> ret;
    for (int row = 0; row < N; ++row) {
	for (int col = 0; col < N; ++col) {
	    ret[row][col] = std::real(complex_mat[row][col]);
	}
    }
    return ret;
}

template<class T, int N>
Matrix<T, N> ImaginaryMatrix(const Matrix<std::complex<T>, N> &complex_mat)
{
    Matrix<T, N> ret;
    for (int row = 0; row < N; ++row) {
	for (int col = 0; col < N; ++col) {
	    ret[row][col] = std::imag(complex_mat[row][col]);
	}
    }
    return ret;
}

template<class T, int N>
void ZeroOutImaginaryPart(Matrix<std::complex<T>, N> &complex_mat)
{
    for (int row = 0; row < N; ++row) {
	for (int col = 0; col < N; ++col) {
	    complex_mat[row][col] = std::real(complex_mat[row][col]);
	}
    }
}

template<class T, int N>
Matrix<T, N> operator*(std::complex<double> coeff, const Matrix<T, N> &mat)
{
    Matrix<T, N> ret;
    for (int row = 0; row < N; ++row) {
	for (int col = 0; col < N; ++col) {
	    ret[row][col] = coeff * mat[row][col];
	}
    }
    return ret;
}

template<class T, int N>
Matrix<T, N> Matrix<T, N>::Dagger() const
{
    Matrix<T, N> ret;
    for (int row = 0; row < N; ++row) {
	for (int col = 0; col < N; ++col) {
	    ret[row][col] = std::conj((*this)[col][row]);
	}
    }
    return ret;
}
template<class T, int N>
T Matrix<T, N>::Trace() const
{
    T ret = 0;
    for (int i = 0; i < N; ++i) {
	ret += (*this)[i][i];
    }
    return ret;
}

#endif
