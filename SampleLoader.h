#ifndef SAMPLE_LOADER_H
#define SAMPLE_LOADER_H

#include <vector>
#include <string>

// For individual applications, write a subclass of SampleLoader
// to parse whatever files we have.
template<class T>
class SampleLoader
{
public:
    virtual T LoadSample(int conf) const = 0;
};


class WilsonMatrix;

class WilsonMatrixLoader : public SampleLoader<WilsonMatrix>
{
private:
    const char* format;
public:
    WilsonMatrixLoader(const char* format) : format(format) {}
    WilsonMatrix LoadSample(int conf) const;
};


class DoubleWilsonMatrix;

class DoubleWilsonMatrixLoader : public SampleLoader<DoubleWilsonMatrix>
{
private:
    const char* format;
public:
    DoubleWilsonMatrixLoader(const char* format) : format(format) {}
    DoubleWilsonMatrix LoadSample(int conf) const;
};


template<class T> 
class BinningSampleLoader
{
public: 
    virtual T LoadAverageSample(const std::vector<int> &confs) const = 0;
};

class BinningWilsonMatrixLoader : public BinningSampleLoader<WilsonMatrix>
{
private:
    WilsonMatrixLoader inner_loader;

public:
    BinningWilsonMatrixLoader(const char* format) : inner_loader(format) {}
    WilsonMatrix LoadAverageSample(const std::vector<int> &confs) const;
};

class BinningDoubleWilsonMatrixLoader : public BinningSampleLoader<DoubleWilsonMatrix>
{
private:
    DoubleWilsonMatrixLoader inner_loader;

public:
    BinningDoubleWilsonMatrixLoader(const char* format) : inner_loader(format) {}
    DoubleWilsonMatrix LoadAverageSample(const std::vector<int> &confs) const;
};



#endif
