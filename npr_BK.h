#ifndef NPR_BK_H
#define NPR_BK_H

#include "JackknifeDatabase.h"
#include "Matrix.h"
#include "NPRSettings.h"

struct NPRResults
{
    JackknifeDatabase<std::complex<double>> jack_Z_BK;
};

NPRResults npr_BK(NPRSettings &sett);

#endif
