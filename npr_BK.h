#ifndef NPR_BK_H
#define NPR_BK_H

#include "JackknifeDatabase.h"
#include "Matrix.h"
#include "NPRSettings.h"

struct NPRResults
{
    JackknifeDatabase<std::complex<double>> jack_Z_BK;
};

void npr_BK(NPRSettings &sett, char* BK_lat_src);

void npr_DSeq2(NPRSettings &sett, char* BK_lat_src);

#endif
