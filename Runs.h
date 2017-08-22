#ifndef RUNS_H
#define RUNS_H

#include "NPR_SimpleG1_WithSubtractions.h"

NPRResults NPR_SimpleG1_WithSubtractions_16cubed();

NPRResults NPR_SimpleG1_WithSubtractions_32cubed_Test(Scheme scheme, Parity parity);

NPRResults NPR_SimpleG1_WithSubtractions_24cubed_low(Scheme scheme, Parity parity);
NPRResults NPR_SimpleG1_WithSubtractions_24cubed_low_TwoMomsForG1(Scheme scheme, Parity parity);

NPRResults NPR_SimpleG1_WithSubtractions_24cubed_high(Scheme scheme, Parity parity);

NPRResults NPR_SimpleG1_WithSubtractions_32cubed_low(Scheme scheme, Parity parity);
NPRResults NPR_SimpleG1_WithSubtractions_32cubed_low_TwoMoms10hits(Scheme scheme, Parity parity);

NPRResults DoStepScaling_24_32(Scheme scheme_fourq, Scheme scheme_Zq, Parity parity);

#endif
