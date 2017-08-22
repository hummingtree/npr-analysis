#include "NPR_Utils.h"

using namespace std;

array<DoubleWilsonMatrix, 7> BuildProjectorSpinColorStructures(Parity parity)
{
    DoubleWilsonMatrix VV_diag;
    DoubleWilsonMatrix AA_diag;
    DoubleWilsonMatrix VA_diag;
    DoubleWilsonMatrix AV_diag;
    DoubleWilsonMatrix VV_mixed;
    DoubleWilsonMatrix AA_mixed;
    DoubleWilsonMatrix VA_mixed;
    DoubleWilsonMatrix AV_mixed;
    for (int mu = 0; mu < 4; ++mu) {
	SpinMatrix V = SpinMatrix::Gamma(mu);
	SpinMatrix A = SpinMatrix::GammaMuGamma5(mu);
	VV_diag += DoubleWilsonMatrix::Construct(V, V, COLOR_STRUCTURE_DIAGONAL);
	AA_diag += DoubleWilsonMatrix::Construct(A, A, COLOR_STRUCTURE_DIAGONAL);
	VA_diag += DoubleWilsonMatrix::Construct(V, A, COLOR_STRUCTURE_DIAGONAL);
	AV_diag += DoubleWilsonMatrix::Construct(A, V, COLOR_STRUCTURE_DIAGONAL);
	VV_mixed += DoubleWilsonMatrix::Construct(V, V, COLOR_STRUCTURE_MIXED);
	AA_mixed += DoubleWilsonMatrix::Construct(A, A, COLOR_STRUCTURE_MIXED);
	VA_mixed += DoubleWilsonMatrix::Construct(V, A, COLOR_STRUCTURE_MIXED);
	AV_mixed += DoubleWilsonMatrix::Construct(A, V, COLOR_STRUCTURE_MIXED);
    }

    DoubleWilsonMatrix LL_diag;
    DoubleWilsonMatrix LR_diag;
    DoubleWilsonMatrix LL_mixed;
    DoubleWilsonMatrix LR_mixed;

    switch (parity) {
	case BOTH_PARITIES:
	    LL_diag = VV_diag + AA_diag - VA_diag - AV_diag;
	    LR_diag = VV_diag - AA_diag + VA_diag - AV_diag;
	    LL_mixed = VV_mixed + AA_mixed - VA_mixed - AV_mixed;
	    LR_mixed = VV_mixed - AA_mixed + VA_mixed - AV_mixed;
	    break;

	case POSITIVE_PARITY:
	    LL_diag = VV_diag + AA_diag;
	    LR_diag = VV_diag - AA_diag;
	    LL_mixed = VV_mixed + AA_mixed;
	    LR_mixed = VV_mixed - AA_mixed;
	    break;

	case NEGATIVE_PARITY:
	    // seems we need an overall minus sign here
	    LL_diag = -(-VA_diag - AV_diag);
	    LR_diag = -(VA_diag - AV_diag);
	    LL_mixed = -(-VA_mixed - AV_mixed);
	    LR_mixed = -(VA_mixed - AV_mixed);
	    break;

	default:
	    printf("unknown parity!\n");
	    exit(-1);
    }

    // The spin - color structure for the seven projectors used in the DeltaS = 1 NPR
    // These are defined in https://rbc.phys.columbia.edu/rbc_ukqcd/individual_postings/qliu/2011/k2pipiNPR/k2pipiNPR_2.pdf
    array<DoubleWilsonMatrix, 7> projectors;
    projectors[0] = LL_diag;
    projectors[1] = LL_mixed;
    projectors[2] = LL_diag;
    projectors[3] = LR_diag;
    projectors[4] = LR_mixed;
    projectors[5] = LR_diag;
    projectors[6] = LR_mixed;
    return projectors;
}

// The G1-like projector
DoubleWilsonMatrix BuildHFProjectorSpinColorStructure(Parity parity)
{
    DoubleWilsonMatrix VV_mixed;
    DoubleWilsonMatrix AV_mixed;
    for (int mu = 0; mu < 4; ++mu) {
	SpinMatrix V = SpinMatrix::Gamma(mu);
	SpinMatrix A = SpinMatrix::GammaMuGamma5(mu);
	VV_mixed += DoubleWilsonMatrix::Construct(V, V, COLOR_STRUCTURE_MIXED);
	AV_mixed += DoubleWilsonMatrix::Construct(A, V, COLOR_STRUCTURE_MIXED);
    }

    DoubleWilsonMatrix LV_mixed;

    switch (parity) {
	case BOTH_PARITIES:
	    LV_mixed = VV_mixed - AV_mixed;
	    break;

	case POSITIVE_PARITY:
	    LV_mixed = VV_mixed;
	    break;

	case NEGATIVE_PARITY:
	    // seems we need an overall minus sign here
	    LV_mixed = -(-AV_mixed);
	    break;

	default:
	    printf("unknown parity!\n");
	    exit(-1);
    }

    return LV_mixed;
}

static DoubleWilsonMatrix BuildQslashProjector(const double q[4], double qsq, double sign, 
    ColorStructure color_struct, Parity parity)
{
    assert(sign == +1.0 || sign == -1.0);
    SpinMatrix qslash = SpinMatrix::Slash(q);
    SpinMatrix qslash_g5 = qslash * SpinMatrix::Gamma5();
    if (parity == POSITIVE_PARITY) {
        return (1 / qsq) * 
	    (DoubleWilsonMatrix::Construct(qslash, qslash, color_struct)
	    + sign * DoubleWilsonMatrix::Construct(qslash_g5, qslash_g5, color_struct));
    } else {
        printf("\n\n\n\nWARNING: UNTESTED NEGATIVE PARITY QSLASH PROJECTORS!!!!!!!!!!\n");
        return (1 / qsq) * 
	    (DoubleWilsonMatrix::Construct(qslash_g5, qslash, color_struct)
	    + sign * DoubleWilsonMatrix::Construct(qslash, qslash_g5, color_struct));
    }
}

static DoubleWilsonMatrix BuildQslashProjectorOnePlus(const double q[4], double qsq, Parity parity)
{
    return BuildQslashProjector(q, qsq, +1.0, COLOR_STRUCTURE_DIAGONAL, parity);
}

static DoubleWilsonMatrix BuildQslashProjectorOneMinus(const double q[4], double qsq, Parity parity)
{
    return BuildQslashProjector(q, qsq, -1.0, COLOR_STRUCTURE_DIAGONAL, parity);
}

static DoubleWilsonMatrix BuildQslashProjectorTwoPlus(const double q[4], double qsq, Parity parity)
{
    return BuildQslashProjector(q, qsq, +1.0, COLOR_STRUCTURE_MIXED, parity);
}

static DoubleWilsonMatrix BuildQslashProjectorTwoMinus(const double q[4], double qsq, Parity parity)
{
    return BuildQslashProjector(q, qsq, -1.0, COLOR_STRUCTURE_MIXED, parity);
}


array<DoubleWilsonMatrix, 7> BuildQslashProjectorSpinColorStructures(const double q[4], double qsq, Parity parity)
{
    array<DoubleWilsonMatrix, 7> ret;

    const int Nc = 3;

    // (27,1) projector, Eq. (99) in Lehner & Sturm 1104.4948
    ret[0] = (1.0 / (64 * Nc * (Nc + 1))) * BuildQslashProjectorOnePlus(q, qsq, parity);

    const double coeff = 1.0 / (32 * Nc * (Nc*Nc - 1));
    
    // (8, 8) projectors, Eq. 100
    ret[5] = coeff * (Nc * BuildQslashProjectorOneMinus(q, qsq, parity) - BuildQslashProjectorTwoMinus(q, qsq, parity));
    ret[6] = coeff * (-1.0 * BuildQslashProjectorOneMinus(q, qsq, parity) + Nc * BuildQslashProjectorTwoMinus(q, qsq, parity));

    // (8,1) projectors, Eq. 102
    ret[1] = coeff * ((3 * Nc - 2) * BuildQslashProjectorOnePlus(q, qsq, parity) + (2 * Nc - 3) * BuildQslashProjectorTwoPlus(q, qsq, parity));
    ret[2] = coeff * ((2 * Nc - 3) * BuildQslashProjectorOnePlus(q, qsq, parity) + (3 * Nc - 2) * BuildQslashProjectorTwoPlus(q, qsq, parity));
    ret[3] = ret[5];
    ret[4] = ret[6];

    return ret;
}

array<DoubleWilsonMatrix, 7> BuildDaiqianQslashProjectorSpinColorStructures(const double q[4], double qsq, Parity parity)
{
    array<DoubleWilsonMatrix, 7> ret;

    // (27, 1)
    ret[0] = BuildQslashProjectorOnePlus(q, qsq, parity);

    // (8,1) 
    ret[1] = BuildQslashProjectorOnePlus(q, qsq, parity);
    ret[2] = BuildQslashProjectorTwoPlus(q, qsq, parity);
    ret[3] = BuildQslashProjectorOneMinus(q, qsq, parity);
    ret[4] = BuildQslashProjectorTwoMinus(q, qsq, parity);

    // (8, 8)
    ret[5] = BuildQslashProjectorOneMinus(q, qsq, parity);
    ret[6] = BuildQslashProjectorTwoMinus(q, qsq, parity);


    return ret;
}


void TestQslashProjectors()
{
    const double q[4] = { 1.0, 2.0, 3.0, 4.0 };
    const double qsq = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    Parity parity = POSITIVE_PARITY;

    DoubleWilsonMatrix P1plus = BuildQslashProjectorOnePlus(q, qsq, parity);
    DoubleWilsonMatrix P1minus = BuildQslashProjectorOneMinus(q, qsq, parity);
    DoubleWilsonMatrix P2plus = BuildQslashProjectorTwoPlus(q, qsq, parity);
    DoubleWilsonMatrix P2minus = BuildQslashProjectorTwoMinus(q, qsq, parity);

    const vector<DoubleWilsonMatrix> basic_Ps { P1plus, P1minus, P2plus, P2minus };
    const vector<const char*> basic_P_names { "P1plus", "P1minus", "P2plus", "P2minus" };

    SpinMatrix qslash = SpinMatrix::Slash(q);
    SpinMatrix qslashg5 = qslash * SpinMatrix::Gamma5();

    DoubleWilsonMatrix qqdiag = (1.0 / qsq) * DoubleWilsonMatrix::Construct(qslash, qslash, COLOR_STRUCTURE_DIAGONAL);
    DoubleWilsonMatrix qqmixed = (1.0 / qsq) * DoubleWilsonMatrix::Construct(qslash, qslash, COLOR_STRUCTURE_MIXED);
    DoubleWilsonMatrix q5q5diag = (1.0 / qsq) * DoubleWilsonMatrix::Construct(qslashg5, qslashg5, COLOR_STRUCTURE_DIAGONAL);
    DoubleWilsonMatrix q5q5mixed = (1.0 / qsq) * DoubleWilsonMatrix::Construct(qslashg5, qslashg5, COLOR_STRUCTURE_MIXED);

    DoubleWilsonMatrix ggdiag;
    DoubleWilsonMatrix ggmixed;
    DoubleWilsonMatrix g5g5diag;
    DoubleWilsonMatrix g5g5mixed;
    for (int mu = 0; mu < 4; ++mu) {
        SpinMatrix gmu = SpinMatrix::Gamma(mu);
        SpinMatrix gmug5 = gmu * SpinMatrix::Gamma5();
        ggdiag += DoubleWilsonMatrix::Construct(gmu, gmu, COLOR_STRUCTURE_DIAGONAL);
        ggmixed += DoubleWilsonMatrix::Construct(gmu, gmu, COLOR_STRUCTURE_MIXED);
        g5g5diag += DoubleWilsonMatrix::Construct(gmug5, gmug5, COLOR_STRUCTURE_DIAGONAL);
        g5g5mixed += DoubleWilsonMatrix::Construct(gmug5, gmug5, COLOR_STRUCTURE_MIXED);
    }

    const vector<DoubleWilsonMatrix> tests { ggdiag, ggmixed, g5g5diag, g5g5mixed, qqdiag, qqmixed, q5q5diag, q5q5mixed };
    const vector<const char*> test_names { "ggdiag", "ggmixed", "g5g5diag", "g5g5mixed", "qqdiag", "qqmixed", "q5q5diag", "q5q5mixed" };

    for (int p = 0; p < 4; ++p) {
        printf("basic projector %s:\n", basic_P_names[p]);
        DoubleWilsonMatrix proj = basic_Ps[p];
        for (int t = 0; t < 8; ++t) {
            complex<double> result =basic_Ps[p].Project(tests[t]);
            printf("  on test state %s = %f + i %f\n", test_names[t], result.real(), result.imag());
        }
        printf("\n");
    }
}


