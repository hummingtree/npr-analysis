#include "Runs.h"

#include "NPR_SimpleG1_WithSubtractions.h"
#include "StepScaling.h"

using namespace std;

NPRResults NPR_SimpleG1_WithSubtractions_16cubed()
{
    NPRSettings sett;
    //const char* dir = "/home/gregm/fitting/NPR/G1_NPR/data/16x32x16_b2.13_ms0.032_ml0.01/justG1_NPR_mass0.0100";
    sett.dir = "/home/gregm/fitting/NPR/G1_NPR/data/16x32x16_b2.13_ms0.032_ml0.01/NPR_full_good_mass0.01";
    //sett.sub_dir = "/home/gregm/fitting/NPR/G1_NPR/data/16x32x16_b2.13_ms0.032_ml0.01/c1_0.0_justG1_NPR_mass0.0100";
    sett.sub_dir = sett.dir;
    sett.c1_str = "";
    for (int conf = 1000; conf <= 8000; conf += 40) sett.confs.push_back(conf);
    sett.mom1 = {{ 1, 1, 2, 4 }}; 
    sett.mom2 = {{ 2, 1, 2, -2 }};
    sett.Ls = {{ 16, 16, 16, 32 }};
    sett.parity = NEGATIVE_PARITY;
    sett.scheme = SchemeGammaMu;
    sett.do_disconnected = true;
    sett.do_subtractions = true;
    sett.enforce_reality = true;
    sett.bin_size = 8;

    return NPR_SimpleG1_WithSubtractions(sett);
}

NPRResults NPR_SimpleG1_WithSubtractions_32cubed_Test(Scheme scheme, Parity parity)
{
    NPRSettings sett;
    sett.dir = "/home/gregm/fitting/NPR/G1_NPR/data/DSDR_b1.75_ms0.045_mu0.001_32x64x32/zbai_NPR_data_5hits_Mobius_converted";
    sett.sub_dir = "/home/gregm/fitting/NPR/G1_NPR/data/DSDR_b1.75_ms0.045_mu0.001_32x64x32/justG1_NPR_mass0.0100";
    sett.c1_str = "";
    for (int conf = 1066; conf <= 1066; conf += 4) sett.confs.push_back(conf);
    sett.mom1 = {{ -4, -2, 2, 0 }}; 
    sett.mom2 = {{ 0, -2, 4, -4 }};
    sett.Ls = {{ 32, 32, 32, 64 }};
    sett.parity = parity;
    sett.scheme = scheme;
    sett.do_disconnected = true;
    sett.do_subtractions = true;
    sett.enforce_reality = false;
    sett.bin_size = 1;

    return NPR_SimpleG1_WithSubtractions(sett);
}

NPRResults NPR_SimpleG1_WithSubtractions_24cubed_low(Scheme scheme, Parity parity)
{
    NPRSettings sett;
    sett.dir = sett.sub_dir = "/home/gregm/fitting/NPR/G1_NPR/data/24nt64_ms0.04_mu0.005/combined_ANL_BNL_fullG1_NPR_mass0.0100";
    sett.c1_str = "";
    for (int conf = 1000; conf <= 8910; conf += 10) sett.confs.push_back(conf);
    sett.mom1 = {{ 0, 2, 2, 0 }};
    sett.mom2 = {{ 2, 2, 0, 0 }};
    sett.Ls = {{ 24, 24, 24, 64 }};
    sett.parity = parity;
    sett.scheme = scheme;
    sett.do_disconnected = true;
    sett.do_subtractions = true;
    sett.enforce_reality = true;
    sett.bin_size = 8;

    return NPR_SimpleG1_WithSubtractions(sett);
}

NPRResults NPR_SimpleG1_WithSubtractions_24cubed_low_TwoMomsForG1(Scheme scheme, Parity parity)
{
    NPRSettings sett;
    sett.dir = sett.sub_dir = "/home/gregm/fitting/NPR/G1_NPR/data/24nt64_ms0.04_mu0.005/combined_ANL_BNL_fullG1_NPR_mass0.0100";
    sett.c1_str = "";
    for (int conf = 1000; conf <= 8900/*8910*/; conf += 10) {
        if (!(2010 <= conf && conf <= 2070)) {
            sett.confs.push_back(conf);
        }
    }
    sett.mom1 = {{ 0, 2, 2, 0 }};
    sett.mom2 = {{ 2, 2, 0, 0 }};
    sett.Ls = {{ 24, 24, 24, 64 }};
    sett.parity = parity;
    sett.scheme = scheme;
    sett.do_disconnected = true;
    sett.do_subtractions = true;
    sett.enforce_reality = true;
    sett.bin_size = 8;

    NPRSettings sett2 = sett;
    sett2.dir = sett2.sub_dir = "/home/gregm/fitting/NPR/G1_NPR/data/24nt64_ms0.04_mu0.005/low_justG1_newmom_mass0.0100";
    sett2.c1_str = "_c1_-0.331";
    sett2.mom1 = {{ 2, 0, 0, -4 }};
    sett2.mom2 = {{ 0, 0, -2, -4 }};

    return NPR_SimpleG1_WithSubtractions_TwoMomsForG1(sett, sett2);
}

NPRResults NPR_SimpleG1_WithSubtractions_24cubed_high(Scheme scheme, Parity parity)
{
    NPRSettings sett;
    sett.dir = sett.sub_dir = "/home/gregm/fitting/NPR/G1_NPR/data/24nt64_ms0.04_mu0.005/combined_ANL_BNL_fullG1_NPR_mass0.0100";
    sett.c1_str = "";
    for (int conf = 1000; conf <= 8910; conf += 10) {
        sett.confs.push_back(conf);
    }
    sett.mom1 = {{ 2, 4, -2, 0 }};
    sett.mom2 = {{ 4, 2, 2, 0 }};
    sett.Ls = {{ 24, 24, 24, 64 }};
    sett.parity = parity;
    sett.scheme = scheme;
    sett.do_disconnected = true;
    sett.do_subtractions = true;
    sett.enforce_reality = true;
    sett.bin_size = 8;

    return NPR_SimpleG1_WithSubtractions(sett);
}

NPRResults NPR_SimpleG1_WithSubtractions_32cubed_low_Ziyuandata(Scheme scheme, Parity parity)
{
    NPRSettings sett;
    sett.dir = "/home/gregm/fitting/NPR/G1_NPR/data/DSDR_b1.75_ms0.045_mu0.001_32x64x32/zbai_NPR_data_5hits_Mobius_converted";
    sett.sub_dir = "/home/gregm/fitting/NPR/G1_NPR/data/DSDR_b1.75_ms0.045_mu0.001_32x64x32/justG1_NPR_mass0.0100";
    sett.c1_str = "";
    for (int conf = 1066; conf <= 2166; conf += 4) sett.confs.push_back(conf);
    sett.mom1 = {{ -4, -2, 2, 0 }};
    sett.mom2 = {{ 0, -2, 4, -4 }};
    sett.Ls = {{ 32, 32, 32, 64 }};
    sett.parity = parity;
    sett.scheme = scheme;
    sett.do_disconnected = true;
    sett.do_subtractions = true;
    sett.enforce_reality = true;
    sett.bin_size = 2;

    return NPR_SimpleG1_WithSubtractions(sett);
}

NPRResults NPR_SimpleG1_WithSubtractions_32cubed_low(Scheme scheme, Parity parity)
{
    NPRSettings sett;
    sett.dir = "/home/gregm/fitting/NPR/G1_NPR/data/DSDR_b1.75_ms0.045_mu0.001_32x64x32/zbai_NPR_data_5hits_Mobius_converted";
    sett.sub_dir = "/home/gregm/fitting/NPR/G1_NPR/data/DSDR_b1.75_ms0.045_mu0.001_32x64x32/justG1_NPR_mass0.0100";
    sett.c1_str = "";
    for (int conf = 1066; conf <= 2166; conf += 4) sett.confs.push_back(conf);
    sett.mom1 = {{ -4, -2, 2, 0 }};
    sett.mom2 = {{ 0, -2, 4, -4 }};
    sett.Ls = {{ 32, 32, 32, 64 }};
    sett.parity = parity;
    sett.scheme = scheme;
    sett.do_disconnected = true;
    sett.do_subtractions = true;
    sett.enforce_reality = true;
    sett.bin_size = 2;

    return NPR_SimpleG1_WithSubtractions(sett);
}
NPRResults NPR_SimpleG1_WithSubtractions_32cubed_low_TwoMoms10hits(Scheme scheme, Parity parity)
{
    NPRSettings sett;
    sett.dir = sett.sub_dir = "/home/gregm/fitting/NPR/G1_NPR/data/DSDR_b1.75_ms0.045_mu0.001_32x64x32/twomoms_10hits_NPR_fullG1_mass0.0100";
    sett.c1_str = "_c1_-0.331";
    for (int conf = 1000; conf <= 2404; conf += 4) {
        if (conf != 1268 && conf != 1844) {
            sett.confs.push_back(conf);
        }
    }
    sett.mom1 = {{ -4, -2, 2, 0 }};
    sett.mom2 = {{ 0, -2, 4, -4 }};
    sett.Ls = {{ 32, 32, 32, 64 }};
    sett.parity = parity;
    sett.scheme = scheme;
    sett.do_disconnected = true;
    sett.do_subtractions = true;
    sett.enforce_reality = true;
    sett.bin_size = 2;

    NPRSettings sett2 = sett;
    sett2.mom1 = {{ 0, 4, -2, -4 }};
    sett2.mom2 = {{ -2, 0, -2, -8 }};
    printf("sett.confs.size() = %d, sett2.confs.size() = %d\n", (int)sett.confs.size(), (int)sett2.confs.size());

    return NPR_SimpleG1_WithSubtractions_TwoMoms(sett, sett2);
}


NPRResults DoStepScaling_24_32(Scheme scheme_fourq, Scheme scheme_Zq, Parity parity)
{
    printf("Start step scaling calculation.\n");

    printf("\n\nGetting 32^3 Z-factors...\n\n");
    //NPRResults results_32_low = NPR_SimpleG1_WithSubtractions_32cubed_low(scheme_fourq, parity);
    NPRResults results_32_low = NPR_SimpleG1_WithSubtractions_32cubed_low_TwoMoms10hits(scheme_fourq, parity);

    printf("\n\nGetting 24^3 low Z-factors...\n\n");
    NPRResults results_24_low = NPR_SimpleG1_WithSubtractions_24cubed_low(scheme_fourq, parity);

    printf("\n\nGetting 24^3 high Z-factors...\n\n");
    NPRResults results_24_high = NPR_SimpleG1_WithSubtractions_24cubed_high(scheme_fourq, parity);

    // properly we should include the Zq's in the jackknifes, but their variance is quite
    // small so we will just treat them as constants
    double Zq_high_over_low_24;
    double Zq_low_32;
    if (scheme_Zq == SchemeGammaMu) {
        printf("USING GammaMu WAVE FUNCTION RENORMALIZATION SCHEME FOR ZQ\n");
        Zq_high_over_low_24 = 0.9918522151257575; // by me, from ZV, GammaMu scheme
        Zq_low_32 = 0.719583; // by me, from ZV, GammaMu scheme, assuming ZV=0.6728 from 1208.4412
    } else { // scheme == SchemeQslash
        printf("USING Qslash WAVE FUNCTION RENORMALIZATION SCHEME FOR ZQ\n");
        Zq_high_over_low_24 = 0.96016057970715996; // by me, from ZV, Qslash scheme
        Zq_low_32 = 0.795763; // by me, from ZV, Qslash scheme, assuming ZV=0.6728 from 1208.4412
    }

    printf("\n\nDone getting all the Z-factors. Now combining to get step-scaled Z-factors...\n\n");

    // 7x7 step scaling
    JackknifeDatabase<Matrix<complex<double>, 7>> jack_Z_lat_32_to_RI_high_7x7 = 
        StepScaleJackknifeZfactors(Zq_low_32*Zq_low_32 * results_32_low.jack_Z_7x7, 
                                   results_24_low.jack_Z_7x7, 
                                   Zq_high_over_low_24*Zq_high_over_low_24 * results_24_high.jack_Z_7x7);

    // 8x8 step scaling
    JackknifeDatabase<Matrix<complex<double>, 8>> jack_Z_lat_32_to_RI_high_8x8 = 
        StepScaleJackknifeZfactors(Zq_low_32*Zq_low_32 * results_32_low.jack_Z_8x8, 
                                   results_24_low.jack_Z_8x8, 
                                   Zq_high_over_low_24*Zq_high_over_low_24 * results_24_high.jack_Z_8x8);

    // Make a jackknife of the Rmatrix
    printf("Doing full step-scaling R-matrix\n");
    JackknifeDatabase<Matrix<complex<double>, 7>> jack_step_scaled_Rmatrix = ComputeRMatrixJackknife(jack_Z_lat_32_to_RI_high_8x8, jack_Z_lat_32_to_RI_high_7x7);

    return NPRResults { jack_Z_lat_32_to_RI_high_7x7, jack_Z_lat_32_to_RI_high_8x8, jack_step_scaled_Rmatrix };

    //printf("A0 with 7x7 Z-factors:\n");
    //DoA0_2p28GeV(jack_Z_lat_32_to_RI_high_7x7);
    //printf("A0 with true R matrix:\n");
    //DoA0_2p28GeV(jack_step_scaled_Rmatrix);
    //printf("A0 with difference matrix:\n");
    //DoA0_2p28GeV(jack_step_scaled_Rmatrix_minus_Zfactors7x7);
}
