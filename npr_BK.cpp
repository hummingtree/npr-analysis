#include "npr_BK.h"

#include "DoubleWilsonMatrix.h"
#include "util.h"
#include "NPR_Utils.h"
#include "A0.h"

#include "contrations_with_BK.h"

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
		} else {
		    sprintf(fourq_diagram_format, "%s/fourq_AA%s_pa%.3f_%.3f_%.3f_%.3f_pb%.3f_%.3f_%.3f_%.3f_traj%%d.bin",
		        sett.sub_dir, sett.c1_str, sett.cont_mom1[0], sett.cont_mom1[1], sett.cont_mom1[2], sett.cont_mom1[3], 
				sett.cont_mom2[0], sett.cont_mom2[1], sett.cont_mom2[2], sett.cont_mom2[3]);
		}
		unamputated_VV_AND_AA_dbs[vvaa] =
		    ConfSampleDatabase<DoubleWilsonMatrix>::LoadBinned(BinningDoubleWilsonMatrixLoader(fourq_diagram_format), sett.confs, sett.bin_size);
    }
    // VV + AA
	ConfSampleDatabase<DoubleWilsonMatrix> unamputated_VVpAA_db = unamputated_VV_AND_AA_dbs[0] + unamputated_VV_AND_AA_dbs[1];
	
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

static JackknifeDatabase<complex<double>> compute_projected_VVpAA_vertex( // projected indicates the vertex is already amputated
        const JackknifeDatabase<DoubleWilsonMatrix> &jack_amputated_VVpAA_vertex,
        const NPRSettings &sett)
{
    printf("computing projected VVpAA vertex...\n");
    JackknifeDatabase<complex<double>> jack_projected_VVpAA_vertex;
    jack_projected_VVpAA_vertex.Resize(sett.Njack);
#pragma omp parallel for
    for (int jack = 0; jack < sett.Njack; ++jack) {
		const DoubleWilsonMatrix &amputated_VVpAA_vertex = jack_amputated_VVpAA_vertex[jack];

		complex<double> projected_VVpAA_vertex; 
		
		// TODO:currently only the GammaMu scheme
		projected_VVpAA_vertex = 
			contrations_with_BK::do_contractions_VVpAA(amputated_VVpAA_vertex, sett.projector_spin_color_structures);

		jack_projected_VVpAA_vertex[jack] = projected_VVpAA_vertex;
    }
    printf("done computing projected VVpAA vertex.\n");

    return jack_projected_VVpAA_vertex;
}

static void build_VVpAA_vertex( // projected and amputated VVpAA vertex
    JackknifeDatabase<complex<double>> &jack_projected_VVpAA_vertex,
    NPRSettings &sett)
{
    sett.Init();

    // Compute external legs
    JackknifeDatabase<WilsonMatrix> jack_prop1, jack_prop2;
    compute_external_legs_jackknife(jack_prop1, jack_prop2, sett);

    // Load and amputate
    JackknifeDatabase<DoubleWilsonMatrix> jack_amputated_VVpAA_vertex = compute_amputated_VVpAA_vertex(jack_prop1, jack_prop2, sett);

    // Project
	jack_projected_VVpAA_vertex = compute_projected_VVpAA_vertex(jack_amputated_VVpAA_vertex, sett);
}

static void build_VA_vertex(
	JackknifeDatabase<complex<double>> &jack_projected_VA_vertex,
	NPRSettings &sett)
{
	sett.Init();

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
				sett.sub_dir, sett.cont_mom1[0], sett.cont_mom1[1], sett.cont_mom1[2], sett.cont_mom1[3], 
				sett.cont_mom2[0], sett.cont_mom2[1], sett.cont_mom2[2], sett.cont_mom2[3]);
		sprintf(A_format, "%s/twoq_gamma%dgamma5_pa%.3f_%.3f_%.3f_%.3f_pb%.3f_%.3f_%.3f_%.3f_traj%%d.dat",
				sett.sub_dir, sett.cont_mom1[0], sett.cont_mom1[1], sett.cont_mom1[2], sett.cont_mom1[3], 
				sett.cont_mom2[0], sett.cont_mom2[1], sett.cont_mom2[2], sett.cont_mom2[3]);
	
		unamputated_V_db[mu] = ConfSampleDatabase<WilsonMatrix>::LoadBinned(BinningWilsonMatrixLoader(V_format), sett.confs, sett.bin_size);
		unamputated_A_db[mu] = ConfSampleDatabase<WilsonMatrix>::LoadBinned(BinningWilsonMatrixLoader(A_format), sett.confs, sett.bin_size);
    }
	
	printf("amputating VA vertex\n");
    JackknifeDatabase<complex<double>> jack_projected_V_vertex;
	jack_projected_V_vertex.Resize(sett.Njack);
    JackknifeDatabase<complex<double>> jack_projected_A_vertex;
	jack_projected_A_vertex.Resize(sett.Njack);

#pragma omp parallel for
    for (int jack = 0; jack < sett.Njack; ++jack) {
		const vector<int> &sample = sett.jackknife_samples[jack];
		const WilsonMatrix &prop1 = jack_prop1[jack];
		const WilsonMatrix &prop2 = jack_prop2[jack];
	
		// Amputate the two quark vertex.
		jack_projected_V_vertex[jack].Zero();
		jack_projected_A_vertex[jack].Zero();
		for(int mu = 0; mu < 4; mu++){
			jack_projected_V_vertex[jack] += (Amputate(unamputated_V_db[mu].MeanOnSample(sample), prop1, prop2) * SpinMatrix::Gamma(mu)).Trace();
			jack_projected_A_vertex[jack] += (Amputate(unamputated_A_db[mu].MeanOnSample(sample), prop1, prop2) * SpinMatrix::Gamma(mu) * SpinMatrix::Gamma5()).Trace();
		}
    }
    printf("done amputating VA vertex...\n");

	// Project vertex
	// TODO: Currently only GammaMu scheme.
	

}

NPRResults npr_BK(NPRSettings &sett)
{
	JackknifeDatabase<complex<double>> &jack_projected_VVpAA_vertex;
	JackknifeDatabase<complex<double>> &jack_projected_VA_vertex;
	
	build_VVpAA_vertex(jack_projected_VVpAA_vertex, sett);
	build_VA_vertex()(jack_projected_VA_vertex, sett);
}


#endif
