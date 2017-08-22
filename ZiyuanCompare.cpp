#if 0
#include "Matrix.h"
#include "DoubleWilsonMatrix.h"
#include "util.h"





#include "ZiyuanCompareAutogeneratedContractions.h"



static std::array<DoubleWilsonMatrix, 7> BuildProjectorSpinColorStructures(Parity parity)
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
    DoubleWilsonMatrix LV_mixed;

    switch (parity) {
	case BOTH_PARITIES:
	    LL_diag = VV_diag + AA_diag - VA_diag - AV_diag;
	    LR_diag = VV_diag - AA_diag + VA_diag - AV_diag;
	    LL_mixed = VV_mixed + AA_mixed - VA_mixed - AV_mixed;
	    LR_mixed = VV_mixed - AA_mixed + VA_mixed - AV_mixed;
	    LV_mixed = VV_mixed - AV_mixed;
	    break;

	case POSITIVE_PARITY:
	    LL_diag = VV_diag + AA_diag;
	    LR_diag = VV_diag - AA_diag;
	    LL_mixed = VV_mixed + AA_mixed;
	    LR_mixed = VV_mixed - AA_mixed;
	    LV_mixed = VV_mixed;
	    break;

	case NEGATIVE_PARITY:
	    // seems we need an overall minus sign here
	    LL_diag = -(-VA_diag - AV_diag);
	    LR_diag = -(VA_diag - AV_diag);
	    LL_mixed = -(-VA_mixed - AV_mixed);
	    LR_mixed = -(VA_mixed - AV_mixed);
	    LV_mixed = -AV_mixed;
	    break;

	default:
	    printf("unknown parity!\n");
	    exit(-1);
    }

    // The spin - color structure for the seven projectors used in the DeltaS = 1 NPR
    // These are defined in https://rbc.phys.columbia.edu/rbc_ukqcd/individual_postings/qliu/2011/k2pipiNPR/k2pipiNPR_2.pdf
    std::array<DoubleWilsonMatrix, 7> projectors;
    projectors[0] = LL_diag;
    projectors[1] = LL_mixed;
    projectors[2] = LL_diag;
    projectors[3] = LR_diag;
    projectors[4] = LR_mixed;
    projectors[5] = LR_diag;
    projectors[6] = LR_mixed;
    return projectors;
}


static WilsonMatrix LoadPropagator(const char *dir, int traj, const int mom[4])
{
    char filename[1024];
    sprintf(filename, "%s/external_leg_p%d_%d_%d_%d_traj%d.dat",
	dir, mom[0], mom[1], mom[2], mom[3], traj);

    WilsonMatrix ret;
    ret.LoadFromFile(filename);
    return ret;
}


static DoubleWilsonMatrix LoadFourQuarkDiagram(const char *dir, const char *diagram_name,
    int traj, const int mom1[4], const int mom2[4])
{
    char filename[1024];
    sprintf(filename, "%s/fourq_%s_pa%d_%d_%d_%d_pb%d_%d_%d_%d_traj%d.dat",
	dir, diagram_name, mom1[0], mom1[1], mom1[2], mom1[3],
	mom2[0], mom2[1], mom2[2], mom2[3], traj);

    return DoubleWilsonMatrix(filename);
}

static DoubleWilsonMatrix LoadVAFourQuarkDiagram(const char *dir, const char *diagram_name,
    VAStructure va, int traj, const int mom1[4], const int mom2[4])
{
    char full_diagram_name[128];
    sprintf(full_diagram_name, "%s_%s", diagram_name, VA_names[va]);
    return LoadFourQuarkDiagram(dir, full_diagram_name, traj, mom1, mom2);
}

static std::array<DoubleWilsonMatrix, 4> LoadLRFourQuarkDiagrams(const char *dir, const char *diagram_name,
    Parity parity, int traj, const int mom1[4], const int mom2[4])
{
    DoubleWilsonMatrix VV = LoadVAFourQuarkDiagram(dir, diagram_name, GammaVV, traj, mom1, mom2);
    DoubleWilsonMatrix AA = LoadVAFourQuarkDiagram(dir, diagram_name, GammaAA, traj, mom1, mom2);
    DoubleWilsonMatrix VA = LoadVAFourQuarkDiagram(dir, diagram_name, GammaVA, traj, mom1, mom2);
    DoubleWilsonMatrix AV = LoadVAFourQuarkDiagram(dir, diagram_name, GammaAV, traj, mom1, mom2);

    std::array<DoubleWilsonMatrix, 4> ret;
    switch (parity) {
	case BOTH_PARITIES:
	    ret[GammaLL] = VV + AA - VA - AV;
	    ret[GammaRR] = VV + AA + VA + AV;
	    ret[GammaLR] = VV - AA + VA - AV;
	    ret[GammaRL] = VV - AA - VA + AV;
	    break;

	case POSITIVE_PARITY:
	    ret[GammaLL] = VV + AA;
	    ret[GammaRR] = VV + AA;
	    ret[GammaLR] = VV - AA;
	    ret[GammaRL] = VV - AA;
	    break;

	case NEGATIVE_PARITY:
	    ret[GammaLL] = -VA - AV;
	    ret[GammaRR] = VA + AV;
	    ret[GammaLR] = VA - AV;
	    ret[GammaRL] = -VA + AV;
	    break;
    }

    return ret;
}

void ZiyuanCompare()
{
    printf("\n------------- ZiyuanCompare start ---------------\n\n");

    const int mom1[] = { 1, 1, 2, 4 };
    const int mom2[] = { 2, 1, 2, -2 };

    const Parity parity = POSITIVE_PARITY;

    //const char* dir = "../../explore/data/Ziyuan-compare";
    const char* dir = "../data/16x32x16_b2.13_ms0.032_ml0.01/for_Ziyuan_compare";
    const int conf = 3620;

    WilsonMatrix prop1 = LoadPropagator(dir, conf, mom1);
    WilsonMatrix prop2 = LoadPropagator(dir, conf, mom2);

    std::array<std::array<DoubleWilsonMatrix, 4>, 6> unamputated_diagrams;
    for (int diagram = 0; diagram < 6; ++diagram) {
	unamputated_diagrams[diagram] = LoadLRFourQuarkDiagrams(dir, diagram_names[diagram], parity, conf, mom1, mom2);
    }


    std::array<std::array<DoubleWilsonMatrix, 4>, 6> amputated_diagrams;
    for (int diagram = 0; diagram < 6; ++diagram) {
	for (int lr = 0; lr < 4; ++lr) {
	    amputated_diagrams[diagram][lr] = Amputate(unamputated_diagrams[diagram][lr], prop1, prop2, prop1, prop2);
	}
    }

    // Taken from https://rbc.phys.columbia.edu/rbc_ukqcd/individual_postings/qliu/2011/k2pipiNPR/k2pipiNPR_2.pdf
    // I have verified this tree level mixing matrix with a free lattice calculation.
    // The Z factor matrix is given by
    // Z / Zq ^ 2 = F inv(M)
    // where M is the matrix of mixings(projections of amputated Green's functions)
    // and F is the tree level mixing matrix below.
    const Matrix<std::complex<double>, 7> tree_level_mixings = {{{
        {{{3072, 3072, 0, 0, 0, 0, 0}}},
        {{{537.6, -230.4, 1152, 0, 0, 0, 0}}},
        {{{-230.4, 537.6, 384, 0, 0, 0, 0}}},
        {{{0, 0, 0, 1152, 384, 3456, 1152}}},
        {{{0, 0, 0, 384, 1152, 1152, 3456}}},
        {{{0, 0, 0, 1152, 384, 0, 0}}},
        {{{0, 0, 0, 384, 1152, 0, 0}}}
    }}};

    std::array<DoubleWilsonMatrix, 7> projector_spin_color_structures = BuildProjectorSpinColorStructures(parity);

    Matrix<std::complex<double>, 7> mixings = 
	ContractionsZiyuanCompare::DoAllDeltaSEquals1ContractionsZiyuanCompare(amputated_diagrams, projector_spin_color_structures);

    Matrix<std::complex<double>, 7> inv_Zfactors = mixings * tree_level_mixings.Inverse();

    inv_Zfactors.Print("inv_Zfactors");

    printf("\n------------- ZiyuanCompare end ---------------\n\n");
}
#endif