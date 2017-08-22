#ifndef JACKKNIFE_DATABASE_H
#define JACKKNIFE_DATABASE_H

#include "SampleLoader.h"
#include "JackknifeError.h"

#include <vector>
#include <complex>
#include <cassert>

///////////////// SAMPLEDATABASE ///////////////////

// A useful base class for storing a set of samples
// and allowing some arithmetic to be done on all
// of them in parallel (for example we form a sample database
// which is the sum of two others).
template<class T>
class SampleDatabase
{
public:
    std::vector<T> samples;

    SampleDatabase() : samples() {} // I shouldn't have to have this
    T& operator[](int i) { return samples[i]; }
    const T& operator[](int i) const { return samples[i]; }
    int Size() const { return samples.size(); }
    void Add(const T &sample) { samples.push_back(sample); }
    void Resize(int size) { samples.resize(size); }
};



///////////////// CONFSAMPLEDATABASE ///////////////////

template<class T> class JackknifeDatabase;

// Samples are loaded directly into a ConfSampleDatabase
template<class T>
class ConfSampleDatabase : public SampleDatabase<T>
{
public:
    ConfSampleDatabase() {}

    static ConfSampleDatabase<T> Load(const SampleLoader<T> &loader, const std::vector<int> &confs);

    static ConfSampleDatabase<T> LoadBinned(const BinningSampleLoader<T> &loader, const std::vector<int> &confs, int bin_size);

    T MeanOnSample(const std::vector<int> &sample) const;

    JackknifeDatabase<T> ComputeJackknifeMeans(const std::vector<std::vector<int>> &jackknife_samples) const;
};


std::vector<std::vector<int>> MakeJackknifeSamples(int N, int block_size);


template<class T>
ConfSampleDatabase<T> ConfSampleDatabase<T>::Load(const SampleLoader<T> &loader, const std::vector<int> &confs)
{
    ConfSampleDatabase<T> ret;
    printf("ConfSampleDatabase: going to load %d samples... ", (int)confs.size());
    for (int conf : confs) {
	ret.samples.push_back(loader.LoadSample(conf));
    }
    printf("done.\n");
    return ret;
}

template<class T>
ConfSampleDatabase<T> ConfSampleDatabase<T>::LoadBinned(const BinningSampleLoader<T> &loader, const std::vector<int> &confs, int bin_size)
{
    printf("ConfSampleDatabase: going to load %d samples with bin_size=%d... ", (int)confs.size(), bin_size);
    assert((confs.size() % bin_size) == 0);
    int Nbins = confs.size() / bin_size;
    ConfSampleDatabase<T> ret;
    for (int bin = 0; bin < Nbins; ++bin) {
	std::vector<int> bin_confs;
	for (int elem = 0; elem < bin_size; ++elem) {
	    bin_confs.push_back(confs[bin * bin_size + elem]);
	}
	ret.samples.push_back(loader.LoadAverageSample(bin_confs));
    }
    printf("done.\n");
    return ret;
}


template<class T>
T ConfSampleDatabase<T>::MeanOnSample(const std::vector<int> &sample_indices) const
{
    assert(sample_indices.size() >= 1);
    T ret = this->samples[sample_indices[0]];
    for (unsigned i = 1; i < sample_indices.size(); ++i) {
	ret += this->samples[sample_indices[i]];
    }
    ret *= 1.0 / sample_indices.size();
    return ret;
}

template<class T>
JackknifeDatabase<T> ConfSampleDatabase<T>::ComputeJackknifeMeans(const std::vector<std::vector<int>> &jackknife_samples) const
{
    JackknifeDatabase<T> ret;
    for (const std::vector<int> &sample_indices : jackknife_samples) {
	ret.Add(MeanOnSample(sample_indices));
    }
    return ret;
}


///////////////// ARITHMETIC FOR CONFSAMPLEDATABASE //////////////////




template<class T> ConfSampleDatabase<T> operator+(const ConfSampleDatabase<T> &a, const ConfSampleDatabase<T> &b)
{
    assert(a.Size() == b.Size());
    ConfSampleDatabase<T> ret;
    for (int i = 0; i < a.Size(); ++i) ret.Add(a[i] + b[i]);
    return ret;
}

template<class T> ConfSampleDatabase<T> operator-(const ConfSampleDatabase<T> &a, const ConfSampleDatabase<T> &b)
{
    assert(a.Size() == b.Size());
    ConfSampleDatabase<T> ret;
    for (int i = 0; i < a.Size(); ++i) ret.Add(a[i] - b[i]);
    return ret;
}

template<class T> void operator+=(ConfSampleDatabase<T> &a, const ConfSampleDatabase<T> &b)
{
    assert(a.Size() == b.Size());
    for (int i = 0; i < a.Size(); ++i) a[i] += b[i];
}

template<class T> void operator-=(ConfSampleDatabase<T> &a, const ConfSampleDatabase<T> &b)
{
    assert(a.Size() == b.Size());
    for (int i = 0; i < a.Size(); ++i) a[i] -= b[i];
}

template<class T> ConfSampleDatabase<T> operator*(std::complex<double> coeff, const ConfSampleDatabase<T> &a)
{
    ConfSampleDatabase<T> ret;
    for (int i = 0; i < a.Size(); ++i) ret.Add(coeff * a[i]);
    return ret;
}

template<class T> void operator*=(ConfSampleDatabase<T> &a, std::complex<double> coeff)
{
    for (int i = 0; i < a.Size(); ++i) a[i] *= coeff;
}

// unary minus
template<class T> ConfSampleDatabase<T> operator-(const ConfSampleDatabase<T> &a)
{
    ConfSampleDatabase<T> ret;
    for (int i = 0; i < a.Size(); ++i) ret.Add(-a[i]);
    return ret;
}


///////////////// JACKKNIFESAMPLEDATABASE ///////////////////

// Stores a set of jackknife sample means
template<class T>
class JackknifeDatabase : public SampleDatabase<T>
{
public:
    T CentralValue() const;
    T Error() const;
};

template<class T>
T JackknifeDatabase<T>::CentralValue() const
{
    return this->samples[0];
}


template<class T>
T JackknifeDatabase<T>::Error() const
{
    return JackknifeError(this->samples);
}

void PrintDoubleWithError(const char *name, const JackknifeDatabase<double> &jack_db, bool newline = true);

void PrintComplexWithError(const char *name, const JackknifeDatabase<std::complex<double> > &jack_db, bool newline = true);

template<std::size_t N>
void PrintComplexArrayWithError(const char* name, const std::array<JackknifeDatabase<std::complex<double>>, N> &jack_dbs)
{
    printf("%s:\n", name);
    for (unsigned i = 0; i < N; ++i) PrintComplexWithError("", jack_dbs[i]);
}

template<std::size_t N>
void PrintComplexArrayWithError(const char* name, const JackknifeDatabase<std::array<std::complex<double>, N>> &jack_db)
{
    printf("%s:\n", name);
    std::array<std::complex<double>, N> centrals = jack_db.CentralValue();
    std::array<std::complex<double>, N> errs = jack_db.Error();
    for (unsigned i = 0; i < N; ++i) {
        printf("[%d]    %e + i %e   +/-   %e + i %e\n", i, centrals[i].real(), centrals[i].imag(), errs[i].real(), errs[i].imag());
    }
}

class SpinMatrix;
void PrintSpinMatrixWithError(const char *name, const JackknifeDatabase<SpinMatrix> &jack_db);




///////////////// ARITHMETIC FOR JACKKNIFEDATABASE //////////////////




template<class T> JackknifeDatabase<T> operator+(const JackknifeDatabase<T> &a, const JackknifeDatabase<T> &b)
{
    assert(a.Size() == b.Size());
    JackknifeDatabase<T> ret;
    for (int i = 0; i < a.Size(); ++i) ret.Add(a[i] + b[i]);
    return ret;
}

template<class T, std::size_t N> JackknifeDatabase<std::array<T, N>> operator+(
        const JackknifeDatabase<std::array<T, N>> &a, const JackknifeDatabase<std::array<T, N>> &b)
{
    assert(a.Size() == b.Size());
    JackknifeDatabase<std::array<T, N>> ret;
    for (int i = 0; i < a.Size(); ++i) {
        std::array<T, N> arr;
        for (unsigned j = 0; j < N; ++j) {
            arr[j] = a[i][j] + b[i][j];
        }
        ret.Add(arr);
    }
    return ret;
}


template<class T> JackknifeDatabase<T> operator-(const JackknifeDatabase<T> &a, const JackknifeDatabase<T> &b)
{
    assert(a.Size() == b.Size());
    JackknifeDatabase<T> ret;
    for (int i = 0; i < a.Size(); ++i) ret.Add(a[i] - b[i]);
    return ret;
}

template<class T> void operator+=(JackknifeDatabase<T> &a, const JackknifeDatabase<T> &b)
{
    assert(a.Size() == b.Size());
    for (int i = 0; i < a.Size(); ++i) a[i] += b[i];
}

template<class T> void operator-=(JackknifeDatabase<T> &a, const JackknifeDatabase<T> &b)
{
    assert(a.Size() == b.Size());
    for (int i = 0; i < a.Size(); ++i) a[i] -= b[i];
}

template<class T> JackknifeDatabase<T> operator*(std::complex<double> coeff, const JackknifeDatabase<T> &a)
{
    JackknifeDatabase<T> ret;
    for (int i = 0; i < a.Size(); ++i) ret.Add(coeff * a[i]);
    return ret;
}

template<class T, std::size_t N> JackknifeDatabase<std::array<T, N>> operator*(std::complex<double> coeff, const JackknifeDatabase<std::array<T, N>> &a)
{
    JackknifeDatabase<std::array<T, N>> ret;
    for (int i = 0; i < a.Size(); ++i) {
        std::array<T, N> arr;
        for (unsigned j = 0; j < N; ++j) {
            arr[j] = coeff * a[i][j];
        }
        ret.Add(arr);
    }
    return ret;
}


template<class T> void operator*=(JackknifeDatabase<T> &a, std::complex<double> coeff)
{
    for (int i = 0; i < a.Size(); ++i) a[i] *= coeff;
}

// unary minus
template<class T> JackknifeDatabase<T> operator-(const JackknifeDatabase<T> &a)
{
    JackknifeDatabase<T> ret;
    for (int i = 0; i < a.Size(); ++i) ret.Add(-a[i]);
    return ret;
}

#endif
