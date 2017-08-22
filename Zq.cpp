#include "Zq.h"

#include "util.h"
#include "JackknifeDatabase.h"
#include "WilsonMatrix.h"
#include <complex>
#include <array>
#include <unordered_set>

using namespace std;


static void DoZq(
    const char* ZVZA_dir,
    const char* leg_dir,
    const vector<int> &confs,
    const int mom1[4],
    const int mom2[4],
    const int Ls[4],
    const int bin_size,
    const double known_ZV,
    const double known_ZA)
{
    const int Nconf = confs.size();
    printf("Nconf = %d\n", Nconf);
    printf("bin_size = %d\n", bin_size);
    assert((Nconf % bin_size) == 0);
    const int Nbins = Nconf / bin_size;
    printf("Nbins = %d\n", Nbins);

    double p1[4], p2[4], q[4];
    double qsq = 0;
    for (int mu = 0; mu < 4; ++mu) {
	p1[mu] = 2 * PI * mom1[mu] / Ls[mu];
	p2[mu] = 2 * PI * mom2[mu] / Ls[mu];
	q[mu] = p1[mu] - p2[mu];
	qsq += q[mu] * q[mu];
    }
    SpinMatrix qslash = SpinMatrix::Slash(q);

    // Load external legs
    printf("loading external legs...\n");
    char prop1_format[512];
    char prop2_format[512];
    sprintf(prop1_format, "%s/external_leg_p%d_%d_%d_%d_traj%%d.dat", leg_dir, mom1[0], mom1[1], mom1[2], mom1[3]);
    sprintf(prop2_format, "%s/external_leg_p%d_%d_%d_%d_traj%%d.dat", leg_dir, mom2[0], mom2[1], mom2[2], mom2[3]);

    ConfSampleDatabase<WilsonMatrix> prop1_db = ConfSampleDatabase<WilsonMatrix>::LoadBinned(BinningWilsonMatrixLoader(prop1_format), confs, bin_size);
    ConfSampleDatabase<WilsonMatrix> prop2_db = ConfSampleDatabase<WilsonMatrix>::LoadBinned(BinningWilsonMatrixLoader(prop2_format), confs, bin_size);

    // Load ZV and ZA diagrams
    printf("loading ZV/ZA diagrams...\n");
    array<ConfSampleDatabase<WilsonMatrix>, 4> unamputated_V_db;
    array<ConfSampleDatabase<WilsonMatrix>, 4> unamputated_A_db;
    for (int mu = 0; mu < 4; ++mu) {
	char V_format[512];
	char A_format[512];
	sprintf(V_format, "%s/twoq_gamma%d_pa%d_%d_%d_%d_pb%d_%d_%d_%d_traj%%d.dat",
	    ZVZA_dir, mu, mom1[0], mom1[1], mom1[2], mom1[3], mom2[0], mom2[1], mom2[2], mom2[3]);
	sprintf(A_format, "%s/twoq_gamma%d_gamma5_pa%d_%d_%d_%d_pb%d_%d_%d_%d_traj%%d.dat",
	    ZVZA_dir, mu, mom1[0], mom1[1], mom1[2], mom1[3], mom2[0], mom2[1], mom2[2], mom2[3]);

	unamputated_V_db[mu] = ConfSampleDatabase<WilsonMatrix>::LoadBinned(BinningWilsonMatrixLoader(V_format), confs, bin_size);
	unamputated_A_db[mu] = ConfSampleDatabase<WilsonMatrix>::LoadBinned(BinningWilsonMatrixLoader(A_format), confs, bin_size);
    }


    printf("doing jackknife computations...\n");
    vector<vector<int>> jackknife_samples = MakeJackknifeSamples(Nbins, 1);
    const int Njack = jackknife_samples.size();
    printf("Njack = %d\n", Njack);

    JackknifeDatabase<complex<double>> jack_V_proj_gamma_mu;
    JackknifeDatabase<complex<double>> jack_A_proj_gamma_mu;
    JackknifeDatabase<complex<double>> jack_V_proj_qslash;
    JackknifeDatabase<complex<double>> jack_A_proj_qslash;
    JackknifeDatabase<double> jack_Zq_gamma_mu_from_V;
    JackknifeDatabase<double> jack_Zq_gamma_mu_from_A;
    JackknifeDatabase<double> jack_Zq_qslash_from_V;
    JackknifeDatabase<double> jack_Zq_qslash_from_A;
    jack_V_proj_gamma_mu.Resize(Njack);
    jack_A_proj_gamma_mu.Resize(Njack);
    jack_V_proj_qslash.Resize(Njack);
    jack_A_proj_qslash.Resize(Njack);
    jack_Zq_gamma_mu_from_V.Resize(Njack);
    jack_Zq_gamma_mu_from_A.Resize(Njack);
    jack_Zq_qslash_from_V.Resize(Njack);
    jack_Zq_qslash_from_A.Resize(Njack);

//#pragma omp parallel for
    for (unsigned jack = 0; jack < jackknife_samples.size(); ++jack) {
	const vector<int> &sample = jackknife_samples[jack];

	// amputate diagrams
	WilsonMatrix prop1 = prop1_db.MeanOnSample(sample);
	WilsonMatrix prop2 = prop2_db.MeanOnSample(sample);

	array<WilsonMatrix, 4> amputated_V;
	array<WilsonMatrix, 4> amputated_A;
	for (int mu = 0; mu < 4; ++mu) {
	    WilsonMatrix unamputated_V = unamputated_V_db[mu].MeanOnSample(sample);
	    WilsonMatrix unamputated_A = unamputated_A_db[mu].MeanOnSample(sample);
	    amputated_V[mu] = Amputate(unamputated_V, prop1, prop2);
	    amputated_A[mu] = Amputate(unamputated_A, prop1, prop2);
	}

	WilsonMatrix V_dot_gmu;
	WilsonMatrix A_dot_g5gmu;
	WilsonMatrix V_dot_q;
	WilsonMatrix A_dot_q;
	for (int mu = 0; mu < 4; ++mu) {
	    V_dot_gmu += amputated_V[mu] * SpinMatrix::Gamma(mu);
	    A_dot_g5gmu += amputated_A[mu] * SpinMatrix::Gamma5() * SpinMatrix::Gamma(mu);
	    V_dot_q += q[mu] * amputated_V[mu];
	    A_dot_q += q[mu] * amputated_A[mu];
	}

	WilsonMatrix V_dot_q_times_qslash = V_dot_q * qslash;
	WilsonMatrix A_dot_q_times_g5_qslash = A_dot_q * SpinMatrix::Gamma5() * qslash;

	// project diagrams
	complex<double> V_proj_gamma_mu = (1.0 / 48.0) * V_dot_gmu.Trace();
	complex<double> A_proj_gamma_mu = (1.0 / 48.0) * A_dot_g5gmu.Trace();
	complex<double> V_proj_qslash = (1.0 / (12.0 * qsq)) * V_dot_q_times_qslash.Trace();
	complex<double> A_proj_qslash = (1.0 / (12.0 * qsq)) * A_dot_q_times_g5_qslash.Trace();

	double Zq_gamma_mu_from_V = known_ZV * V_proj_gamma_mu.real();
	double Zq_gamma_mu_from_A = known_ZA * A_proj_gamma_mu.real();
	double Zq_qslash_from_V = known_ZV * V_proj_qslash.real();
	double Zq_qslash_from_A = known_ZA * A_proj_qslash.real();

	// store stuff
	jack_V_proj_gamma_mu[jack] = V_proj_gamma_mu;
	jack_A_proj_gamma_mu[jack] = A_proj_gamma_mu;
	jack_V_proj_qslash[jack] = V_proj_qslash;
	jack_A_proj_qslash[jack] = A_proj_qslash;

	jack_Zq_gamma_mu_from_V[jack] = Zq_gamma_mu_from_V;
	jack_Zq_gamma_mu_from_A[jack] = Zq_gamma_mu_from_A;
	jack_Zq_qslash_from_V[jack] = Zq_qslash_from_V;
	jack_Zq_qslash_from_A[jack] = Zq_qslash_from_A;
    }

    printf("known_ZV = %f\n", known_ZV);
    printf("known_ZA = %f\n", known_ZA);
    printf("\n");
    PrintComplexWithError("V_proj_gamma_mu", jack_V_proj_gamma_mu);
    PrintComplexWithError("A_proj_gamma_mu", jack_A_proj_gamma_mu);
    PrintDoubleWithError("Zq_gamma_mu_from_V", jack_Zq_gamma_mu_from_V);
    PrintDoubleWithError("Zq_gamma_mu_from_A", jack_Zq_gamma_mu_from_A);
    printf("\n");
    PrintComplexWithError("V_proj_qslash", jack_V_proj_qslash);
    PrintComplexWithError("A_proj_qslash", jack_A_proj_qslash);
    PrintDoubleWithError("Zq_qslash_from_V", jack_Zq_qslash_from_V);
    PrintDoubleWithError("Zq_qslash_from_A", jack_Zq_qslash_from_A);
}



void DoZq_24cubed_high()
{
    const char* ZVZA_dir = "/home/gregm/fitting/NPR/G1_NPR/data/24nt64_ms0.04_mu0.005/ZVZA_mass0.0100";
    const char* leg_dir = "/home/gregm/fitting/NPR/G1_NPR/data/24nt64_ms0.04_mu0.005/combined_ANL_BNL_fullG1_NPR_mass0.0100";
    vector<int> confs;
    unordered_set<int> excluded{ 1620, 1640, 2020, 2040, 2060, 8960 };
    for (int conf = 1000; conf <= 8961; conf += 20) {
	if (excluded.count(conf) > 0) continue;
	confs.push_back(conf);
    }
    const int mom1[] = { 2, 4, -2, 0 };
    const int mom2[] = { 4, 2, 2, 0 };
    int Ls[] = { 24, 24, 24, 64 };
    const int bin_size = 1;
    const double known_ZV = 1; // I don't know the true value
    const double known_ZA = 1; // I don't know the true value

    DoZq(ZVZA_dir, leg_dir, confs, mom1, mom2, Ls, bin_size, known_ZV, known_ZA);
}

void DoZq_24cubed_low()
{
    const char* ZVZA_dir = "/home/gregm/fitting/NPR/G1_NPR/data/24nt64_ms0.04_mu0.005/ZVZA_mass0.0100";
    const char* leg_dir = "/home/gregm/fitting/NPR/G1_NPR/data/24nt64_ms0.04_mu0.005/combined_ANL_BNL_fullG1_NPR_mass0.0100";
    vector<int> confs;
    unordered_set<int> excluded{ 1620, 1640, 2020, 2040, 2060, 8960 };
    for (int conf = 1000; conf <= 8961; conf += 20) {
	if (excluded.count(conf) > 0) continue;
	confs.push_back(conf);
    }
    const int mom1[] = { 0, 2, 2, 0 };
    const int mom2[] = { 2, 2, 0, 0 };
    int Ls[] = { 24, 24, 24, 64 };
    const int bin_size = 1;
    const double known_ZV = 1; // I don't know the true value
    const double known_ZA = 1; // I don't know the true value

    DoZq(ZVZA_dir, leg_dir, confs, mom1, mom2, Ls, bin_size, known_ZV, known_ZA);
}


void DoZq_32cubed()
{
    const char* ZVZA_dir = "/home/gregm/fitting/NPR/G1_NPR/data/DSDR_b1.75_ms0.045_mu0.001_32x64x32/ZVZA_mass0.0100";
    const char* leg_dir = "/home/gregm/fitting/NPR/G1_NPR/data/DSDR_b1.75_ms0.045_mu0.001_32x64x32/justG1_NPR_mass0.0100";
    vector<int> confs;
    for (int conf = 1066; conf <= 1442; conf += 8) confs.push_back(conf);    
    const int mom1[] = { -4, -2, 2, 0 };
    const int mom2[] = { 0, -2, 4, -4 };
    int Ls[] = { 32, 32, 32, 64 };
    const int bin_size = 1;
    const double known_ZV = 0.6728;  // 1208.4412 gives 0.6728(80)
    const double known_ZA = 0.68778; // 1208.4412 gives 0.68778(34)

    DoZq(ZVZA_dir, leg_dir, confs, mom1, mom2, Ls, bin_size, known_ZV, known_ZA);
}