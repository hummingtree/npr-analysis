#include <cstdio>
#include <fstream>

#include "tests.h"
// #include "Runs.h"

#include "npr_BK.h"

void run_BK_24I()
{
    NPRSettings sett;
    //const char* dir = "/home/gregm/fitting/NPR/G1_NPR/data/16x32x16_b2.13_ms0.032_ml0.01/justG1_NPR_mass0.0100";
    sett.dir = "../data/24x64x16I/BK_mass0.0200";
//	sett.dir = "../data/greg_test_24I";
    //sett.sub_dir = "/home/gregm/fitting/NPR/G1_NPR/data/16x32x16_b2.13_ms0.032_ml0.01/c1_0.0_justG1_NPR_mass0.0100";
    sett.sub_dir = sett.dir;
    sett.c1_str = "";
    for (int conf = 1000; conf <= 1199; conf += 20) sett.confs.push_back(conf);
    sett.mom1 = {{ 0, 2, 2, 0 }}; 
    sett.mom2 = {{ 2, 2, 0, 0}};
    sett.cont_mom1 = {{ -1.695, 0., 1.695, 0. }}; 
    sett.cont_mom2 = {{ 0., 1.695, 1.695, 0. }}; 
    sett.Ls = {{ 24, 24, 24, 64 }};
    sett.parity = POSITIVE_PARITY;
    sett.scheme = SchemeGammaMu;
    sett.do_disconnected = true;
    sett.do_subtractions = true;
    sett.enforce_reality = true;
    sett.bin_size = 1;
	
    npr_BK(sett, NULL);
    
	sett.scheme = SchemeQslash;

    npr_BK(sett, NULL);
}

void run_BK_32I()
{
    NPRSettings sett;
    //const char* dir = "/home/gregm/fitting/NPR/G1_NPR/data/16x32x16_b2.13_ms0.032_ml0.01/justG1_NPR_mass0.0100";
    sett.dir = "../data/32x64x16I/BK_mass0.0080";
//	sett.dir = "../data/greg_test_24I";
    //sett.sub_dir = "/home/gregm/fitting/NPR/G1_NPR/data/16x32x16_b2.13_ms0.032_ml0.01/c1_0.0_justG1_NPR_mass0.0100";
    sett.sub_dir = sett.dir;
    sett.c1_str = "";
    for (int conf = 1000; conf <= 1199; conf += 20) sett.confs.push_back(conf);
    sett.mom1 = {{ 0, 2, 2, 0 }}; 
    sett.mom2 = {{ 2, 2, 0, 0}};
    sett.cont_mom1 = {{ -1.692, 0., 1.692, 0. }}; 
    sett.cont_mom2 = {{ 0., 1.692, 1.692, 0. }}; 
    sett.Ls = {{ 32, 32, 32, 64 }};
    sett.parity = POSITIVE_PARITY;
    sett.scheme = SchemeGammaMu;
    sett.do_disconnected = true;
    sett.do_subtractions = true;
    sett.enforce_reality = true;
    sett.bin_size = 1;
	
    npr_BK(sett, NULL);
    
	sett.scheme = SchemeQslash;

    npr_BK(sett, NULL);
}

void run_BK_24ID()
{
    NPRSettings sett;
    //const char* dir = "/home/gregm/fitting/NPR/G1_NPR/data/16x32x16_b2.13_ms0.032_ml0.01/justG1_NPR_mass0.0100";
    sett.dir = "../data/24x64x24ID/BK_mass0.0011";
//	sett.dir = "../data/greg_test_24I";
    //sett.sub_dir = "/home/gregm/fitting/NPR/G1_NPR/data/16x32x16_b2.13_ms0.032_ml0.01/c1_0.0_justG1_NPR_mass0.0100";
    sett.sub_dir = sett.dir;
    sett.c1_str = "";
    for (int conf = 300; conf <= 440; conf += 4) sett.confs.push_back(conf);
    sett.mom1 = {{ 0, 2, 2, 0 }}; 
    sett.mom2 = {{ 2, 2, 0, 0}};
    sett.cont_mom1 = {{ -3., 0., 3., 0. }}; 
    sett.cont_mom2 = {{ 0., 3., 3., 0. }}; 
    sett.Ls = {{ 24, 24, 24, 64 }};
    sett.parity = POSITIVE_PARITY;
    sett.scheme = SchemeGammaMu;
    sett.do_disconnected = true;
    sett.do_subtractions = true;
    sett.enforce_reality = true;
    sett.bin_size = 3;
	
    npr_BK(sett, "../data/24x64x24ID/BK_24ID");
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

	run_BK_32I();

    return 0;
}
