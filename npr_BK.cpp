#include "npr_BK.h"

#include "DoubleWilsonMatrix.h"
#include "util.h"
#include "NPR_Utils.h"
#include "A0.h"

#include "contractions_with_BK.h"

#include <cassert>

using namespace std;

static const vector<int> non_8_1_operators { 0,   // Q'_1  (27,1)
                                             5,   // Q'_7  (8,8)
                                             6 }; // Q'_8  (8,8)

static void compute_external_legs_jackknife(
        JackknifeDatabase<WilsonMatrix> &jack_prop1,
        JackknifeDatabase<WilsonMatrix> &jack_prop2,
        const NPRSettings &sett)
{
    // Load external legs
    printf("loading external legs...\n");
    char prop1_format[512];
    char prop2_format[512];
	sprintf(prop1_format, "%s/rev_external_leg_p%.3f_%.3f_%.3f_%.3f_traj%%d.dat", sett.dir, sett.cont_mom1[0], sett.cont_mom1[1], sett.cont_mom1[2], sett.cont_mom1[3]);
	sprintf(prop2_format, "%s/external_leg_p%.3f_%.3f_%.3f_%.3f_traj%%d.dat", sett.dir, sett.cont_mom2[0], sett.cont_mom2[1], sett.cont_mom2[2], sett.cont_mom2[3]);
//	sprintf(prop1_format, "%s/external_leg_p%d_%d_%d_%d_traj%%d.dat", sett.dir, sett.mom1[0], sett.mom1[1], sett.mom1[2], sett.mom1[3]);
//	sprintf(prop2_format, "%s/external_leg_p%d_%d_%d_%d_traj%%d.dat", sett.dir, sett.mom2[0], sett.mom2[1], sett.mom2[2], sett.mom2[3]);

    ConfSampleDatabase<WilsonMatrix> prop1_db = ConfSampleDatabase<WilsonMatrix>::LoadBinned(BinningWilsonMatrixLoader(prop1_format), sett.confs, sett.bin_size);
    ConfSampleDatabase<WilsonMatrix> prop2_db = ConfSampleDatabase<WilsonMatrix>::LoadBinned(BinningWilsonMatrixLoader(prop2_format), sett.confs, sett.bin_size);

    printf("doing external leg jackknife...\n");
    jack_prop1.Resize(sett.Njack);
    jack_prop2.Resize(sett.Njack);
#pragma omp parallel for
    for (int jack = 0; jack < sett.Njack; ++jack) {
	const vector<int> &sample = sett.jackknife_samples[jack];
	jack_prop1[jack] = prop1_db.MeanOnSample(sample);
	jack_prop2[jack] = prop2_db.MeanOnSample(sample);
    }
    printf("done external legs\n");
}

static JackknifeDatabase<array<DoubleWilsonMatrix, 5>> compute_amputated_DSeq2_vertex(
        const JackknifeDatabase<WilsonMatrix> &jack_prop1,
        const JackknifeDatabase<WilsonMatrix> &jack_prop2,
        const NPRSettings &sett)
{

// VV +- AA
    printf("loading unamputated VVpmAA vertex...\n");
    array<ConfSampleDatabase<DoubleWilsonMatrix>, 2> unamputated_VV_AND_AA_dbs;
    
	for (int vvaa = 0; vvaa < 2; ++vvaa) {
		char fourq_diagram_format[512];
		if (vvaa == 0) {
		    sprintf(fourq_diagram_format, "%s/fourq_VV%s_pa%.3f_%.3f_%.3f_%.3f_pb%.3f_%.3f_%.3f_%.3f_traj%%d.bin",
		        sett.sub_dir, sett.c1_str, sett.cont_mom1[0], sett.cont_mom1[1], sett.cont_mom1[2], sett.cont_mom1[3], 
				sett.cont_mom2[0], sett.cont_mom2[1], sett.cont_mom2[2], sett.cont_mom2[3]);
		} else {
		    sprintf(fourq_diagram_format, "%s/fourq_AA%s_pa%.3f_%.3f_%.3f_%.3f_pb%.3f_%.3f_%.3f_%.3f_traj%%d.bin",
		        sett.sub_dir, sett.c1_str, sett.cont_mom1[0], sett.cont_mom1[1], sett.cont_mom1[2], sett.cont_mom1[3], 
				sett.cont_mom2[0], sett.cont_mom2[1], sett.cont_mom2[2], sett.cont_mom2[3]);
		}
		unamputated_VV_AND_AA_dbs[vvaa] =
		    ConfSampleDatabase<DoubleWilsonMatrix>::LoadBinned(BinningDoubleWilsonMatrixLoader(fourq_diagram_format), sett.confs, sett.bin_size);
    }
	ConfSampleDatabase<DoubleWilsonMatrix> unamputated_VVpAA_db = unamputated_VV_AND_AA_dbs[0] + unamputated_VV_AND_AA_dbs[1];
	ConfSampleDatabase<DoubleWilsonMatrix> unamputated_VVmAA_db = unamputated_VV_AND_AA_dbs[0] - unamputated_VV_AND_AA_dbs[1];
	
// END OF VV +- AA

// SS -+ PP

	printf("loading unamputated SSmpPP vertex...\n");
    array<ConfSampleDatabase<DoubleWilsonMatrix>, 2> unamputated_SS_AND_PP_dbs;
    
	for (int sspp = 0; sspp < 2; ++sspp) {
		char fourq_diagram_format[512];
		if (sspp == 0) {
		    sprintf(fourq_diagram_format, "%s/fourq_SS%s_pa%.3f_%.3f_%.3f_%.3f_pb%.3f_%.3f_%.3f_%.3f_traj%%d.bin",
		        sett.sub_dir, sett.c1_str, sett.cont_mom1[0], sett.cont_mom1[1], sett.cont_mom1[2], sett.cont_mom1[3], 
				sett.cont_mom2[0], sett.cont_mom2[1], sett.cont_mom2[2], sett.cont_mom2[3]);
		} else {
		    sprintf(fourq_diagram_format, "%s/fourq_PP%s_pa%.3f_%.3f_%.3f_%.3f_pb%.3f_%.3f_%.3f_%.3f_traj%%d.bin",
		        sett.sub_dir, sett.c1_str, sett.cont_mom1[0], sett.cont_mom1[1], sett.cont_mom1[2], sett.cont_mom1[3], 
				sett.cont_mom2[0], sett.cont_mom2[1], sett.cont_mom2[2], sett.cont_mom2[3]);
		}
		unamputated_SS_AND_PP_dbs[sspp] =
		    ConfSampleDatabase<DoubleWilsonMatrix>::LoadBinned(BinningDoubleWilsonMatrixLoader(fourq_diagram_format), sett.confs, sett.bin_size);
    }
	ConfSampleDatabase<DoubleWilsonMatrix> unamputated_SSpPP_db = unamputated_SS_AND_PP_dbs[0] + unamputated_SS_AND_PP_dbs[1];
	ConfSampleDatabase<DoubleWilsonMatrix> unamputated_SSmPP_db = unamputated_SS_AND_PP_dbs[0] - unamputated_SS_AND_PP_dbs[1];

// END OF SS -+ PP

// TT

	printf("loading unamputated TT vertex...\n");
    //const char* gammaAs[] = { "V", "A" };
    ConfSampleDatabase<DoubleWilsonMatrix> unamputated_TT_db;
	
	char fourq_diagram_format[512];
	sprintf(fourq_diagram_format, "%s/fourq_TT%s_pa%.3f_%.3f_%.3f_%.3f_pb%.3f_%.3f_%.3f_%.3f_traj%%d.bin",
	    sett.sub_dir, sett.c1_str, sett.cont_mom1[0], sett.cont_mom1[1], sett.cont_mom1[2], sett.cont_mom1[3], 
		sett.cont_mom2[0], sett.cont_mom2[1], sett.cont_mom2[2], sett.cont_mom2[3]);
	unamputated_TT_db =
	    ConfSampleDatabase<DoubleWilsonMatrix>::LoadBinned(BinningDoubleWilsonMatrixLoader(fourq_diagram_format), sett.confs, sett.bin_size);

// END OF TT
    
    printf("amputating DSeq2 vertices\n");
    JackknifeDatabase<array<DoubleWilsonMatrix, 5>> jack_amputated_DSeq2_vertex;
	jack_amputated_DSeq2_vertex.Resize(sett.Njack);	

#pragma omp parallel for
    for (int jack = 0; jack < sett.Njack; ++jack) {
		const vector<int> &sample = sett.jackknife_samples[jack];
		const WilsonMatrix &prop1 = jack_prop1[jack];
		const WilsonMatrix &prop2 = jack_prop2[jack];
	
		// Amputate the  diagram of G1
		jack_amputated_DSeq2_vertex[jack][0] = Amputate(unamputated_VVpAA_db.MeanOnSample(sample), prop1, prop2, prop1, prop2);
		jack_amputated_DSeq2_vertex[jack][1] = Amputate(unamputated_VVmAA_db.MeanOnSample(sample), prop1, prop2, prop1, prop2);
		jack_amputated_DSeq2_vertex[jack][2] = Amputate(unamputated_SSmPP_db.MeanOnSample(sample), prop1, prop2, prop1, prop2);
		jack_amputated_DSeq2_vertex[jack][3] = Amputate(unamputated_SSpPP_db.MeanOnSample(sample), prop1, prop2, prop1, prop2);
		jack_amputated_DSeq2_vertex[jack][4] = Amputate(unamputated_TT_db.MeanOnSample(sample), prop1, prop2, prop1, prop2);
    
//		jack_amputated_DSeq2_vertex[jack][0] = build_BK_gammaMu_pscs(POSITIVE_PARITY)[0];
//		jack_amputated_DSeq2_vertex[jack][1] = build_BK_gammaMu_pscs(POSITIVE_PARITY)[1];
//		jack_amputated_DSeq2_vertex[jack][2] = build_BK_gammaMu_pscs(POSITIVE_PARITY)[2];
//		jack_amputated_DSeq2_vertex[jack][3] = build_BK_gammaMu_pscs(POSITIVE_PARITY)[3];
//		jack_amputated_DSeq2_vertex[jack][4] = build_BK_gammaMu_pscs(POSITIVE_PARITY)[4];
}
    printf("done amputating DSeq2 vertex...\n");

    return jack_amputated_DSeq2_vertex;
}

static JackknifeDatabase<DoubleWilsonMatrix> compute_amputated_VVpAA_vertex(
        const JackknifeDatabase<WilsonMatrix> &jack_prop1,
        const JackknifeDatabase<WilsonMatrix> &jack_prop2,
        const NPRSettings &sett)
{
    printf("loading unamputated VVpAA vertex...\n");
    //const char* gammaAs[] = { "V", "A" };
    array<ConfSampleDatabase<DoubleWilsonMatrix>, 2> unamputated_VV_AND_AA_dbs;
//	array<ConfSampleDatabase<WilsonMatrix>, 2> unamputated_twoq_VA_G1_dbs;
    for (int vvaa = 0; vvaa < 2; ++vvaa) {
	char fourq_diagram_format[512];
	//sprintf(fourq_diagram_format, "%s/fourq_G1%s_%s_pa%d_%d_%d_%d_pb%d_%d_%d_%d_traj%%d.bin",
	//    dir, pc_str, gammaAs[va], mom1[0], mom1[1], mom1[2], mom1[3], mom2[0], mom2[1], mom2[2], mom2[3]);
		if (vvaa == 0) {
		    sprintf(fourq_diagram_format, "%s/fourq_VV%s_pa%.3f_%.3f_%.3f_%.3f_pb%.3f_%.3f_%.3f_%.3f_traj%%d.bin",
		        sett.sub_dir, sett.c1_str, sett.cont_mom1[0], sett.cont_mom1[1], sett.cont_mom1[2], sett.cont_mom1[3], 
				sett.cont_mom2[0], sett.cont_mom2[1], sett.cont_mom2[2], sett.cont_mom2[3]);
//			sprintf(fourq_diagram_format, "%s/fourq_fully_connected_color_diag_V_V%s_pa%d_%d_%d_%d_pb%d_%d_%d_%d_traj%%d.bin",
//		        sett.sub_dir, sett.c1_str, sett.mom1[0], sett.mom1[1], sett.mom1[2], sett.mom1[3], 
//				sett.mom2[0], sett.mom2[1], sett.mom2[2], sett.mom2[3]);
		} else {
//			sprintf(fourq_diagram_format, "%s/fourq_fully_connected_color_diag_A_A%s_pa%d_%d_%d_%d_pb%d_%d_%d_%d_traj%%d.bin",
//		        sett.sub_dir, sett.c1_str, sett.mom1[0], sett.mom1[1], sett.mom1[2], sett.mom1[3], 
//				sett.mom2[0], sett.mom2[1], sett.mom2[2], sett.mom2[3]);

		    sprintf(fourq_diagram_format, "%s/fourq_AA%s_pa%.3f_%.3f_%.3f_%.3f_pb%.3f_%.3f_%.3f_%.3f_traj%%d.bin",
		        sett.sub_dir, sett.c1_str, sett.cont_mom1[0], sett.cont_mom1[1], sett.cont_mom1[2], sett.cont_mom1[3], 
				sett.cont_mom2[0], sett.cont_mom2[1], sett.cont_mom2[2], sett.cont_mom2[3]);
		}
		unamputated_VV_AND_AA_dbs[vvaa] =
		    ConfSampleDatabase<DoubleWilsonMatrix>::LoadBinned(BinningDoubleWilsonMatrixLoader(fourq_diagram_format), sett.confs, sett.bin_size);
    }
    // VV + AA
	ConfSampleDatabase<DoubleWilsonMatrix> unamputated_VVpAA_db = unamputated_VV_AND_AA_dbs[0] - unamputated_VV_AND_AA_dbs[1];
	
//	MakeLeftHandedDiagrams<DoubleWilsonMatrix>(unamputated_fourq_VA_G1_dbs, sett.parity);
    
    printf("amputating VVpAA vertex\n");
    JackknifeDatabase<DoubleWilsonMatrix> jack_amputated_VVpAA_vertex;
    jack_amputated_VVpAA_vertex.Resize(sett.Njack);
#pragma omp parallel for
    for (int jack = 0; jack < sett.Njack; ++jack) {
		const vector<int> &sample = sett.jackknife_samples[jack];
		const WilsonMatrix &prop1 = jack_prop1[jack];
		const WilsonMatrix &prop2 = jack_prop2[jack];
	
		// Amputate the four-quark diagram of G1
		jack_amputated_VVpAA_vertex[jack] = Amputate(unamputated_VVpAA_db.MeanOnSample(sample), prop1, prop2, prop1, prop2);
    }
    printf("done amputating VVpAA vertex...\n");

    return jack_amputated_VVpAA_vertex;
}

static JackknifeDatabase<array<complex<double>, 5>> compute_projected_VVpAA_vertex( // projected indicates the vertex is already amputated
        const JackknifeDatabase<DoubleWilsonMatrix> &jack_amputated_VVpAA_vertex,
        const NPRSettings &sett)
{
    printf("computing projected VVpAA vertex...\n");
    JackknifeDatabase<array<complex<double>, 5>> jack_projected_VVpAA_vertex;
    jack_projected_VVpAA_vertex.Resize(sett.Njack);
#pragma omp parallel for
    for (int jack = 0; jack < sett.Njack; ++jack) {
		const DoubleWilsonMatrix &amputated_VVpAA_vertex = jack_amputated_VVpAA_vertex[jack];

		array<complex<double>, 5> projected_VVpAA_vertex; 
		// TODO:currently only the GammaMu scheme
		if(sett.scheme == SchemeGammaMu){
			BK_pscs gammaMu_pscs = build_BK_gammaMu_pscs(sett.parity); // for qslash scheme. seems redundant but have to b/c of cont_q		
			projected_VVpAA_vertex = 
				contrations_with_BK::do_contractions_VVpAA(amputated_VVpAA_vertex, gammaMu_pscs);
		}else{
			BK_pscs qslash_pscs = build_BK_Qslash_pscs(sett.cont_q, sett.cont_qsq, sett.parity); // for qslash scheme. seems redundant but have to b/c of cont_q		
			projected_VVpAA_vertex = 
				contrations_with_BK::do_contractions_VVpAA(amputated_VVpAA_vertex, qslash_pscs);
		}

		jack_projected_VVpAA_vertex[jack] = projected_VVpAA_vertex;
    }
    printf("done computing projected VVpAA vertex.\n");

    return jack_projected_VVpAA_vertex;
}

static JackknifeDatabase<Matrix<complex<double>, 5>> compute_projected_DSeq2_vertex( // projected indicates the vertex is already amputated
        const JackknifeDatabase<array<DoubleWilsonMatrix, 5>> &jack_amputated_DSeq2_vertex,
        const NPRSettings &sett)
{
    printf("computing projected DSeq2 vertex...\n");

	JackknifeDatabase<Matrix<complex<double>, 5>> jack_projected_DSeq2_vertex;
    jack_projected_DSeq2_vertex.Resize(sett.Njack);
#pragma omp parallel for
    for (int jack = 0; jack < sett.Njack; ++jack) {
		const array<DoubleWilsonMatrix, 5> &amputated_DSeq2_vertex = jack_amputated_DSeq2_vertex[jack];

		Matrix<complex<double>, 5> projected_DSeq2_vertex; 
		if(sett.scheme == SchemeGammaMu){
			BK_pscs gammaMu_pscs = build_BK_gammaMu_pscs(sett.parity); 		
			projected_DSeq2_vertex = 
				contrations_with_BK::do_contractions_DSeq2(amputated_DSeq2_vertex, gammaMu_pscs);
		}else{
			BK_pscs qslash_pscs = build_BK_Qslash_pscs(sett.cont_q, sett.cont_qsq, sett.parity); 		
			projected_DSeq2_vertex = 
				contrations_with_BK::do_contractions_DSeq2(amputated_DSeq2_vertex, qslash_pscs);
		}

		jack_projected_DSeq2_vertex[jack] = projected_DSeq2_vertex;
    }
    printf("done computing projected DSeq2 vertex.\n");

    return jack_projected_DSeq2_vertex;
}

static void build_VVpAA_vertex( // projected and amputated VVpAA vertex
    JackknifeDatabase<array<complex<double>, 5>> &jack_projected_VVpAA_vertex,
    NPRSettings &sett)
{
    // Compute external legs
    JackknifeDatabase<WilsonMatrix> jack_prop1, jack_prop2;
    compute_external_legs_jackknife(jack_prop1, jack_prop2, sett);

    // Load and amputate
    JackknifeDatabase<DoubleWilsonMatrix> jack_amputated_VVpAA_vertex = compute_amputated_VVpAA_vertex(jack_prop1, jack_prop2, sett);

    // Project
	jack_projected_VVpAA_vertex = compute_projected_VVpAA_vertex(jack_amputated_VVpAA_vertex, sett);
}

static void build_DSeq2_vertex( // projected and amputated VVpAA vertex
    JackknifeDatabase<Matrix<complex<double>, 5>> &jack_projected_DSeq2_vertex,
    NPRSettings &sett)
{
    // Compute external legs
    JackknifeDatabase<WilsonMatrix> jack_prop1, jack_prop2;
    compute_external_legs_jackknife(jack_prop1, jack_prop2, sett);

    // Load and amputate
    JackknifeDatabase<array<DoubleWilsonMatrix, 5>> jack_amputated_DSeq2_vertex = compute_amputated_DSeq2_vertex(jack_prop1, jack_prop2, sett);

    // Project
	jack_projected_DSeq2_vertex = compute_projected_DSeq2_vertex(jack_amputated_DSeq2_vertex, sett);
}

static void build_VA_vertex(
	JackknifeDatabase<complex<double>> &jack_projected_V_vertex,
	JackknifeDatabase<complex<double>> &jack_projected_A_vertex,
	NPRSettings &sett)
{
	// Compute external legs
	JackknifeDatabase<WilsonMatrix> jack_prop1, jack_prop2;
	compute_external_legs_jackknife(jack_prop1, jack_prop2, sett);

	// Load and amputate
	printf("Loading ZV/ZA diagrams...\n");
    array<ConfSampleDatabase<WilsonMatrix>, 4> unamputated_V_db;
    array<ConfSampleDatabase<WilsonMatrix>, 4> unamputated_A_db;
	for (int mu = 0; mu < 4; ++mu) {
		
		char V_format[512];
		char A_format[512];

		sprintf(V_format, "%s/twoq_gamma%d_pa%.3f_%.3f_%.3f_%.3f_pb%.3f_%.3f_%.3f_%.3f_traj%%d.dat",
				sett.sub_dir, mu, sett.cont_mom1[0], sett.cont_mom1[1], sett.cont_mom1[2], sett.cont_mom1[3], 
				sett.cont_mom2[0], sett.cont_mom2[1], sett.cont_mom2[2], sett.cont_mom2[3]);
		sprintf(A_format, "%s/twoq_gamma%dgamma5_pa%.3f_%.3f_%.3f_%.3f_pb%.3f_%.3f_%.3f_%.3f_traj%%d.dat",
				sett.sub_dir, mu, sett.cont_mom1[0], sett.cont_mom1[1], sett.cont_mom1[2], sett.cont_mom1[3], 
				sett.cont_mom2[0], sett.cont_mom2[1], sett.cont_mom2[2], sett.cont_mom2[3]);
	
		unamputated_V_db[mu] = ConfSampleDatabase<WilsonMatrix>::LoadBinned(BinningWilsonMatrixLoader(V_format), sett.confs, sett.bin_size);
		unamputated_A_db[mu] = ConfSampleDatabase<WilsonMatrix>::LoadBinned(BinningWilsonMatrixLoader(A_format), sett.confs, sett.bin_size);
    }
	
	printf("amputating VA vertex\n");
	
	jack_projected_V_vertex.Resize(sett.Njack);
	jack_projected_A_vertex.Resize(sett.Njack);

	// TODO: qslash!!!
	SpinMatrix cont_qslash = SpinMatrix::Slash(sett.cont_q);

#pragma omp parallel for
    for (int jack = 0; jack < sett.Njack; ++jack) {
		const vector<int> &sample = sett.jackknife_samples[jack];
		const WilsonMatrix &prop1 = jack_prop1[jack];
		const WilsonMatrix &prop2 = jack_prop2[jack];
	
		// Amputate the two quark vertex.
		jack_projected_V_vertex[jack] = 0.;
		jack_projected_A_vertex[jack] = 0.;
		for(int mu = 0; mu < 4; mu++){
			if(sett.scheme == SchemeGammaMu){
				jack_projected_V_vertex[jack] += (Amputate(unamputated_V_db[mu].MeanOnSample(sample), prop1, prop2) * SpinMatrix::Gamma(mu)).Trace() / 48.;
				jack_projected_A_vertex[jack] += (Amputate(unamputated_A_db[mu].MeanOnSample(sample), prop1, prop2) * SpinMatrix::Gamma5() * SpinMatrix::Gamma(mu)).Trace() / 48.;
			}else{
				// TODO: Implement Qslash scheme
				jack_projected_V_vertex[jack] += (Amputate(unamputated_V_db[mu].MeanOnSample(sample), prop1, prop2) 
													* cont_qslash).Trace() * sett.cont_q[mu]/ (12.*sett.cont_qsq);
				jack_projected_A_vertex[jack] += (Amputate(unamputated_A_db[mu].MeanOnSample(sample), prop1, prop2) 
													* SpinMatrix::Gamma5() * cont_qslash).Trace() * sett.cont_q[mu] / (12.*sett.cont_qsq);
			}
		}
    }
    printf("done amputating VA vertex...\n");

	// Project vertex
	// TODO: Currently only GammaMu scheme.
	

}

void npr_DSeq2(NPRSettings &sett, char* BK_lat_src)
{
	sett.Init();
	
	JackknifeDatabase<Matrix<complex<double>, 5>> jack_projected_DSeq2_vertex;
	jack_projected_DSeq2_vertex.Resize(sett.Njack);
	JackknifeDatabase<complex<double>> jack_projected_V_vertex;
	JackknifeDatabase<complex<double>> jack_projected_A_vertex;

	build_DSeq2_vertex(jack_projected_DSeq2_vertex, sett);
	build_VA_vertex(jack_projected_V_vertex, jack_projected_A_vertex, sett);
	
	PrintComplexWithError("V     ", jack_projected_V_vertex, true);
	PrintComplexWithError("A     ", jack_projected_A_vertex, true);

	PrintRealPartMatrixWithError("real Pi_5x5", jack_projected_DSeq2_vertex);

	JackknifeDatabase<Matrix<complex<double>, 5>> jack_DSeq2_Z;
	jack_DSeq2_Z.Resize(sett.Njack);

	JackknifeDatabase<Matrix<complex<double>, 5>> jack_DSeq2_M;
	jack_DSeq2_M.Resize(sett.Njack);
//
//
	double Z_A = 0.73457; // https://rbc.phys.columbia.edu/rbc_ukqcd/individual_postings/dmurphy/coarse_ensembles/spectrum.pdf page 2

	for (int jack = 0; jack < sett.Njack; ++jack) {
		double VpAd2 = (jack_projected_V_vertex[jack].real() + jack_projected_A_vertex[jack].real()) / 2.;
		jack_DSeq2_M[jack] = (1. / (VpAd2 * VpAd2)) * jack_projected_DSeq2_vertex[jack];
		if(sett.scheme == SchemeGammaMu){
			jack_DSeq2_Z[jack] = sett.tree_level_greens_funcs_DSeq2_5x5_GammaMu * jack_DSeq2_M[jack].Inverse();
		}else{
			jack_DSeq2_Z[jack] = sett.tree_level_greens_funcs_DSeq2_5x5_Qslash * jack_DSeq2_M[jack].Inverse();
		}
	}
	
	PrintRealPartMatrixWithError("real Z_5x5", jack_DSeq2_Z);
//
//	PrintDoubleWithError("ZBK   ", jack_projected_ZBK_vertex, true);

//	printf("VVpAA: %.8f+I%.8f\n", jack_projected_VVpAA_vertex.CentralValue().real(), jack_projected_VVpAA_vertex.CentralValue().imag());
//	printf("V:     %.8f+I%.8f\n", jack_projected_V_vertex.CentralValue().real(), jack_projected_V_vertex.CentralValue().imag());
//	printf("A:     %.8f+I%.8f\n", jack_projected_A_vertex.CentralValue().real(), jack_projected_A_vertex..imag());
	
//	if(BK_lat_src != NULL){	
//		JackknifeDatabase<double> jack_projected_BK_lat_vertex;
//		JackknifeDatabase<double> jack_projected_BK_SMOM_vertex;
//		jack_projected_BK_lat_vertex.Resize(sett.Njack);	
//		jack_projected_BK_SMOM_vertex.Resize(sett.Njack);	
//		ReadJackknifeDatabaseFromFile(jack_projected_BK_lat_vertex, BK_lat_src, sett.Njack);
//		for (int jack = 0; jack < sett.Njack; ++jack) {
//			jack_projected_BK_SMOM_vertex[jack] = jack_projected_BK_lat_vertex[jack] * jack_projected_ZBK_vertex[jack];
//		}
//		PrintDoubleWithError("BK_lat ", jack_projected_BK_lat_vertex, true);
//		PrintDoubleWithError("BK_SMOM", jack_projected_BK_SMOM_vertex, true);
//	}


}

void npr_BK(NPRSettings &sett, char* BK_lat_src)
{
	sett.Init();
	
	JackknifeDatabase<array<complex<double>, 5>> jack_projected_VVpAA_vertex;
	JackknifeDatabase<complex<double>> jack_projected_VVpAA0_vertex;
	JackknifeDatabase<complex<double>> jack_projected_VVpAA1_vertex;
	JackknifeDatabase<complex<double>> jack_projected_VVpAA2_vertex;
	JackknifeDatabase<complex<double>> jack_projected_VVpAA3_vertex;
	JackknifeDatabase<complex<double>> jack_projected_VVpAA4_vertex;
	jack_projected_VVpAA0_vertex.Resize(sett.Njack);
	jack_projected_VVpAA1_vertex.Resize(sett.Njack);
	jack_projected_VVpAA2_vertex.Resize(sett.Njack);
	jack_projected_VVpAA3_vertex.Resize(sett.Njack);
	jack_projected_VVpAA4_vertex.Resize(sett.Njack);
//	JackknifeDatabase<complex<double>> jack_projected_BK_vertex;
	JackknifeDatabase<complex<double>> jack_projected_V_vertex;
	JackknifeDatabase<complex<double>> jack_projected_A_vertex;
	
//	jack_projected_BK_vertex.Resize(Njack);

	build_VVpAA_vertex(jack_projected_VVpAA_vertex, sett);
	build_VA_vertex(jack_projected_V_vertex, jack_projected_A_vertex, sett);

	for (int jack = 0; jack < sett.Njack; ++jack) {
		jack_projected_VVpAA0_vertex[jack] = jack_projected_VVpAA_vertex[jack][0];
		jack_projected_VVpAA1_vertex[jack] = jack_projected_VVpAA_vertex[jack][1];
		jack_projected_VVpAA2_vertex[jack] = jack_projected_VVpAA_vertex[jack][2];
		jack_projected_VVpAA3_vertex[jack] = jack_projected_VVpAA_vertex[jack][3];
		jack_projected_VVpAA4_vertex[jack] = jack_projected_VVpAA_vertex[jack][4];
	}
	
	PrintComplexWithError("VVpAA ", jack_projected_VVpAA0_vertex, true);
	PrintComplexWithError("VVpAA1", jack_projected_VVpAA1_vertex, true);
	PrintComplexWithError("VVpAA2", jack_projected_VVpAA2_vertex, true);
	PrintComplexWithError("VVpAA3", jack_projected_VVpAA3_vertex, true);
	PrintComplexWithError("VVpAA4", jack_projected_VVpAA4_vertex, true);
	PrintComplexWithError("V     ", jack_projected_V_vertex, true);
	PrintComplexWithError("A     ", jack_projected_A_vertex, true);

	JackknifeDatabase<double> jack_projected_ZBK_vertex;
	jack_projected_ZBK_vertex.Resize(sett.Njack);

	for (int jack = 0; jack < sett.Njack; ++jack) {
		double VpAd2 = (jack_projected_V_vertex[jack].real() + jack_projected_A_vertex[jack].real()) / 2.;
		jack_projected_ZBK_vertex[jack] = VpAd2 * VpAd2 / jack_projected_VVpAA0_vertex[jack].real();
	}

	PrintDoubleWithError("ZBK   ", jack_projected_ZBK_vertex, true);

//	printf("VVpAA: %.8f+I%.8f\n", jack_projected_VVpAA_vertex.CentralValue().real(), jack_projected_VVpAA_vertex.CentralValue().imag());
//	printf("V:     %.8f+I%.8f\n", jack_projected_V_vertex.CentralValue().real(), jack_projected_V_vertex.CentralValue().imag());
//	printf("A:     %.8f+I%.8f\n", jack_projected_A_vertex.CentralValue().real(), jack_projected_A_vertex..imag());
	
	if(BK_lat_src != NULL){	
		JackknifeDatabase<double> jack_projected_BK_lat_vertex;
		JackknifeDatabase<double> jack_projected_BK_SMOM_vertex;
		jack_projected_BK_lat_vertex.Resize(sett.Njack);	
		jack_projected_BK_SMOM_vertex.Resize(sett.Njack);	
		ReadJackknifeDatabaseFromFile(jack_projected_BK_lat_vertex, BK_lat_src, sett.Njack);
		for (int jack = 0; jack < sett.Njack; ++jack) {
			jack_projected_BK_SMOM_vertex[jack] = jack_projected_BK_lat_vertex[jack] * jack_projected_ZBK_vertex[jack];
		}
		PrintDoubleWithError("BK_lat ", jack_projected_BK_lat_vertex, true);
		PrintDoubleWithError("BK_SMOM", jack_projected_BK_SMOM_vertex, true);
	}


}

