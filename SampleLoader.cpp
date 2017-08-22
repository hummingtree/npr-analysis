#include "SampleLoader.h"
#include "WilsonMatrix.h"
#include "DoubleWilsonMatrix.h"
#include "util.h"

#include <fstream>
#include <sstream>
#include <string>
#include <cstring>




WilsonMatrix WilsonMatrixLoader::LoadSample(int conf) const
{
    char filename[512];
    sprintf(filename, format, conf);
    //printf("loading %s\n", filename);
    WilsonMatrix ret(filename); 
    return ret;
}


DoubleWilsonMatrix DoubleWilsonMatrixLoader::LoadSample(int conf) const
{
    char filename[512];
    sprintf(filename, format, conf);
    //printf("loading %s\n", filename);
    if (EndsWithBin(filename)) {
	DoubleWilsonMatrix ret;
	ret.LoadFromBinaryFileReversingEndianness(filename);
	return ret;
    } else {
	return DoubleWilsonMatrix(filename);
    }
}



WilsonMatrix BinningWilsonMatrixLoader::LoadAverageSample(const std::vector<int> &confs) const
{
    //printf("BinningWilsonMatrixLoader going to load %lu confs: ", confs.size());
    //for (int conf : confs) printf("%d ", conf);
    //printf("\n");

    WilsonMatrix ret;
    for (int conf : confs) {
	ret += inner_loader.LoadSample(conf);
    }
    ret *= 1.0 / confs.size();
    return ret;
}

DoubleWilsonMatrix BinningDoubleWilsonMatrixLoader::LoadAverageSample(const std::vector<int> &confs) const
{
    //printf("BinningDoubleWilsonMatrixLoader going to load %lu confs: ", confs.size());
    //for (int conf : confs) printf("%d ", conf);
    //printf("\n");

    DoubleWilsonMatrix ret;
    for (int conf : confs) {
	ret += inner_loader.LoadSample(conf);
    }
    ret *= 1.0 / confs.size();
    return ret;
}








