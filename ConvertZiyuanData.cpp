#include "WilsonMatrix.h"
#include "DoubleWilsonMatrix.h"
#include "util.h"

#include <fstream>
#include <array>
using namespace std;


static WilsonMatrix ReadZiyuanExternalLeg(const char* filename)
{
    WilsonMatrix ret;
    ifstream f(filename);
    if (!f) {
	printf("Couldn't open %s\n", filename);
	exit(-1);
    }
    for (int s1 = 0; s1 < 4; s1++) {
	for (int c1 = 0; c1 < 3; c1++) {
	    for (int s2 = 0; s2 < 4; s2++) {
		for (int c2 = 0; c2 < 3; c2++) {
		    double real;
		    double imag;
		    f >> real;
		    f >> imag;
		    if (f.fail()) {
			printf("Failure reading %s\n", filename);
			exit(-1);
		    }
		    ret(s1, c1, s2, c2) = complex<double>(real, imag);
		}
	    }
	}
    }
    f.close();
    return ret;
}

static DoubleWilsonMatrix ReadZiyuanDoubleWilsonMatrixFromStream(ifstream &f)
{
    DoubleWilsonMatrix ret;
    for (int i = 0; i < 144 * 144; ++i) {
	double real, imag;
	f >> real;
	f >> imag;
	// Qi/Ziyuan/Daiqian have the tensor indices in a different order than I do
	int qi_indices[4];
	int j = i;
	qi_indices[0] = j % 12; j /= 12;
	qi_indices[1] = j % 12; j /= 12;
	qi_indices[2] = j % 12; j /= 12;
	qi_indices[3] = j;
	int my_index =
	    qi_indices[3] +
	    12 * (qi_indices[2] +
	    12 * (qi_indices[1] +
	    12 * qi_indices[0]));

	ret.data[my_index] = complex<double>(real, imag);
    }

    return ret;
}

static WilsonMatrix ReadZiyuanOpLineFromStream(ifstream &f)
{
    WilsonMatrix ret;
    for (int i = 0; i < 12 * 12; ++i) {
	double real, imag;
	f >> real;
	f >> imag;
	// Qi/Ziyuan/Daiqian have the tensor indices in a different order than I do
	int my_index = 12 * (i % 12) + (i / 12);
	ret.data[my_index] = complex<double>(real, imag);
    }

    return ret;
}

// order in which fourq diagrams appear in Ziyuan format bulk files
static const array<Diagram, 12> zbai_diagram_order{{
    FULLY_CONNECTED_COLOR_DIAG, FULLY_CONNECTED_COLOR_MIXED, FULLY_CONNECTED_COLOR_DIAG, FULLY_CONNECTED_COLOR_MIXED,
    DISCONNECTED_LOOP_COLOR_DIAG, DISCONNECTED_LOOP_COLOR_MIXED, DISCONNECTED_LOOP_COLOR_DIAG, DISCONNECTED_LOOP_COLOR_MIXED,
    CONNECTED_LOOP_COLOR_DIAG, CONNECTED_LOOP_COLOR_MIXED, CONNECTED_LOOP_COLOR_DIAG, CONNECTED_LOOP_COLOR_MIXED
}};
// order in which VA structures appear in EVEN PARITY Ziyuan format bulk files
static const array<VAStructure, 12> zbai_even_VA_order{{
    GammaVV, GammaVV, GammaAA, GammaAA, // fully connected
    GammaVV, GammaVV, GammaAA, GammaAA, // disconnected loop
    GammaVV, GammaVV, GammaAA, GammaAA  // connected loop
}};
// order in which VA structures appear in ODD PARITY Ziyuan format bulk files
static const array<VAStructure, 12> zbai_odd_VA_order{{
    GammaAV, GammaAV, GammaVA, GammaVA, // fully connected
    GammaAV, GammaAV, GammaVA, GammaVA, // disconnected loop
    GammaAV, GammaAV, GammaVA, GammaVA  // connected loop
}};

static void ReadZiyuanFourqBulkFile(const char* filename, Parity parity,
    array<array<DoubleWilsonMatrix, 4>, 6> &zbai_unamputated_fourq_VA_diagrams)
{
    ifstream f(filename);
    if (!f) {
	printf("couldn't open bulk file %s\n", filename);
	exit(-1);
    }

    array<DoubleWilsonMatrix, 12> zbai_fourq_diagrams;

    // read fully connected diagrams
    for (int i = 0; i < 4; ++i) {
	zbai_fourq_diagrams[i] = ReadZiyuanDoubleWilsonMatrixFromStream(f);
    }

    // read disconnected diagrams
    WilsonMatrix spectator = ReadZiyuanOpLineFromStream(f);
    for (int i = 4; i < 12; ++i) {
	WilsonMatrix op_line = ReadZiyuanOpLineFromStream(f);
	zbai_fourq_diagrams[i] = DoubleWilsonMatrix::OuterProduct(op_line, spectator);
    }

    f.close();

    const array<VAStructure, 12> &VA_order = (parity == POSITIVE_PARITY ? zbai_even_VA_order : zbai_odd_VA_order);

    // figure out the Diagram and VAStructure of each loaded diagram
    for (int i = 0; i < 12; ++i) {
	Diagram diagram = zbai_diagram_order[i];
	VAStructure va = VA_order[i];
	zbai_unamputated_fourq_VA_diagrams[diagram][va] = zbai_fourq_diagrams[i];
    }
}

static void ReadZiyuanTwoqBulkFile(const char* filename, Parity parity,
    array<array<WilsonMatrix, 4>, 6> &zbai_unamputated_twoq_VA_diagrams)
{
    ifstream f(filename);
    if (!f) {
	printf("couldn't open bulk file %s\n", filename);
	exit(-1);
    }

    array<WilsonMatrix, 8> zbai_twoq_diagrams;

    // read disconnected diagrams
    for (int i = 0; i < 8; ++i) {
	zbai_twoq_diagrams[i] = ReadZiyuanOpLineFromStream(f);
    }

    f.close();

    const array<VAStructure, 12> &VA_order = (parity == POSITIVE_PARITY ? zbai_even_VA_order : zbai_odd_VA_order);

    // figure out the Diagram and VAStructure of each loaded diagram
    for (int i = 0; i < 8; ++i) {
	// the +4's below are because we are only doing disconnected diagrams, skipping the connected ones
	Diagram diagram = zbai_diagram_order[i + 4];
	VAStructure va = VA_order[i + 4];
	zbai_unamputated_twoq_VA_diagrams[diagram][va] = zbai_twoq_diagrams[i];
    }
}

void ConvertZiyuanNPRData(const char* zbai_dir, const char* new_dir, int conf, const int mom1[4], const int mom2[4])
{
    printf("converting conf %d...\n", conf);

    // convert external legs

    char zbai_prop1_filename[512];
    char zbai_prop2_filename[512];
    sprintf(zbai_prop1_filename, "%s/traj_%d_legs_p0.txt", zbai_dir, conf);
    sprintf(zbai_prop2_filename, "%s/traj_%d_legs_p1.txt", zbai_dir, conf);
    WilsonMatrix zbai_prop1 = ReadZiyuanExternalLeg(zbai_prop1_filename);
    WilsonMatrix zbai_prop2 = ReadZiyuanExternalLeg(zbai_prop2_filename);
    char new_prop1_filename[512];
    char new_prop2_filename[512];
    sprintf(new_prop1_filename, "%s/external_leg_p%d_%d_%d_%d_traj%d.dat", new_dir, mom1[0], mom1[1], mom1[2], mom1[3], conf);
    sprintf(new_prop2_filename, "%s/external_leg_p%d_%d_%d_%d_traj%d.dat", new_dir, mom2[0], mom2[1], mom2[2], mom2[3], conf);
    zbai_prop1.WriteToFile(new_prop1_filename);
    zbai_prop2.WriteToFile(new_prop2_filename);

    // convert fourq diagrams

    char even_fourq_bulk_file[512];
    char odd_fourq_bulk_file[512];
    sprintf(even_fourq_bulk_file, "%s/traj_%d_k2pipiNPR_even_p0_p1.txt", zbai_dir, conf);
    sprintf(odd_fourq_bulk_file, "%s/traj_%d_k2pipiNPR_odd_p0_p1.txt", zbai_dir, conf);
    
    array<array<DoubleWilsonMatrix, 4>, 6> zbai_unamputated_fourq_VA_diagrams;
    ReadZiyuanFourqBulkFile(even_fourq_bulk_file, POSITIVE_PARITY, zbai_unamputated_fourq_VA_diagrams);
    ReadZiyuanFourqBulkFile(odd_fourq_bulk_file, NEGATIVE_PARITY, zbai_unamputated_fourq_VA_diagrams);

    for (int diagram = 0; diagram < 6; ++diagram) {
	for (int va = 0; va < 4; ++va) {
	    char stem[512];
	    sprintf(stem, "%s/fourq_%s_%s_pa%d_%d_%d_%d_pb%d_%d_%d_%d_traj%d",
		new_dir, diagram_names[diagram], VA_names[va], 
		mom1[0], mom1[1], mom1[2], mom1[3], mom2[0], mom2[1], mom2[2], mom2[3], conf);
	    char filename_dat[512];
	    char filename_bin[512];
	    sprintf(filename_dat, "%s.dat", stem);
	    sprintf(filename_bin, "%s.bin", stem);

	    DoubleWilsonMatrix &dwm = zbai_unamputated_fourq_VA_diagrams[diagram][va];
	    dwm.WriteToFile(filename_dat);
	    dwm.WriteToBinaryFileReversingEndianness(filename_bin);
	}
    }

    // convert twoq diagrams

    char even_twoq_bulk_file[512];
    char odd_twoq_bulk_file[512];
    sprintf(even_twoq_bulk_file, "%s/traj_%d_subtractionCoeff_even_p0_p1.txt", zbai_dir, conf);
    sprintf(odd_twoq_bulk_file, "%s/traj_%d_subtractionCoeff_odd_p0_p1.txt", zbai_dir, conf);

    array<array<WilsonMatrix, 4>, 6> zbai_unamputated_twoq_VA_diagrams;
    ReadZiyuanTwoqBulkFile(even_twoq_bulk_file, POSITIVE_PARITY, zbai_unamputated_twoq_VA_diagrams);
    ReadZiyuanTwoqBulkFile(odd_twoq_bulk_file, NEGATIVE_PARITY, zbai_unamputated_twoq_VA_diagrams);

    for (int diagram = 2 /* note */; diagram < 6; ++diagram) {
	for (int va = 0; va < 4; ++va) {
	    char filename[512];
	    sprintf(filename, "%s/twoq_%s_%s_pa%d_%d_%d_%d_pb%d_%d_%d_%d_traj%d.dat",
		new_dir, diagram_names[diagram], VA_names[va],
		mom1[0], mom1[1], mom1[2], mom1[3], mom2[0], mom2[1], mom2[2], mom2[3], conf);
	    zbai_unamputated_twoq_VA_diagrams[diagram][va].WriteToFile(filename);
	}
    }
}


void DoZiyuanConversion32cubed()
{
    const int mom1[] = { -4, -2, 2, 0 }; // 32^3
    const int mom2[] = { 0, -2, 4, -4 }; // 32^3

    const char* zbai_dir = "/home/gregm/fitting/NPR/G1_NPR/data/DSDR_b1.75_ms0.045_mu0.0042_32x64x32/zbai_NPR_data";
    const char* new_dir = "/home/gregm/fitting/NPR/G1_NPR/data/DSDR_b1.75_ms0.045_mu0.0042_32x64x32/zbai_converted_NPR_data";

    printf("Doing conversion of Ziyuan NPR data. zbai_dir=\n%s\nnew_dir=\n%s\n\n", zbai_dir, new_dir);

    vector<int> confs;
    for (int conf = 1066; conf <= 2170; conf += 4) confs.push_back(conf);

    for (int conf : confs) {
	ConvertZiyuanNPRData(zbai_dir, new_dir, conf, mom1, mom2);
    }
}
