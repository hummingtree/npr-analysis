#include <cstdio>
#include <fstream>

#include "tests.h"
// #include "Runs.h"

#include "npr_BK.h"

void run_BK()
{
    NPRSettings sett;
    //const char* dir = "/home/gregm/fitting/NPR/G1_NPR/data/16x32x16_b2.13_ms0.032_ml0.01/justG1_NPR_mass0.0100";
    sett.dir = "../data/24x64x16I/BK_mass0.0050";
    //sett.sub_dir = "/home/gregm/fitting/NPR/G1_NPR/data/16x32x16_b2.13_ms0.032_ml0.01/c1_0.0_justG1_NPR_mass0.0100";
    sett.sub_dir = sett.dir;
    sett.c1_str = "";
    for (int conf = 1000; conf <= 1101; conf += 10) sett.confs.push_back(conf);
    sett.mom1 = {{ 1, 1, 2, 4 }}; 
    sett.mom2 = {{ 2, 1, 2, -2 }};
    sett.cont_mom1 = {{ -4.5, 0., 4.5, 0. }}; 
    sett.cont_mom2 = {{ 0., 4.5, 4.5, 0. }}; 
    sett.Ls = {{ 24, 24, 24, 64 }};
    sett.parity = POSITIVE_PARITY;
    sett.scheme = SchemeGammaMu;
    sett.do_disconnected = true;
    sett.do_subtractions = true;
    sett.enforce_reality = true;
    sett.bin_size = 1;
	
    npr_BK(sett);
}

int main()
{
    //RunAllTests();
    
    //NPR_SimpleG1_WithSubtractions_16cubed();
    
    //NPR_SimpleG1_WithSubtractions_32cubed_Test(SchemeGammaMu, POSITIVE_PARITY);
    //NPR_SimpleG1_WithSubtractions_32cubed_Test(SchemeGammaMu, NEGATIVE_PARITY);
    //NPR_SimpleG1_WithSubtractions_32cubed_Test(SchemeQslash, POSITIVE_PARITY);
    //NPR_SimpleG1_WithSubtractions_32cubed_Test(SchemeQslash, NEGATIVE_PARITY);

    //NPR_SimpleG1_WithSubtractions_24cubed_high(SchemeQslash, POSITIVE_PARITY);
    //NPR_SimpleG1_WithSubtractions_24cubed_high(SchemeQslash, NEGATIVE_PARITY);
    //NPR_SimpleG1_WithSubtractions_32cubed_low_TwoMoms10hits(SchemeGammaMu);

    //printf("ONE MOM FOR G1:\n");
    //NPR_SimpleG1_WithSubtractions_24cubed_low(SchemeGammaMu, NEGATIVE_PARITY);
    //printf("TWO MOMS FOR G1:\n");
    //NPR_SimpleG1_WithSubtractions_24cubed_low_TwoMomsForG1(SchemeGammaMu, NEGATIVE_PARITY);

    //DoStepScaling_24_32(SchemeGammaMu, SchemeGammaMu, NEGATIVE_PARITY);
    //DoStepScaling_24_32(SchemeQslash, SchemeQslash, NEGATIVE_PARITY);

	run_BK();

    return 0;
}
