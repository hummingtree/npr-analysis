#include "tests.h"
#include "SampleLoader.h"
#include "JackknifeDatabase.h"
#include "util.h"
#include "WilsonMatrix.h"
#include "DoubleWilsonMatrix.h"

#include <cassert>
#include <cstdio>
#include <cmath>
#include <array>



void RunAllTests()
{
    TestJackknifeDatabaseDouble();
    TestJackknifeDatabaseComplex();
    TestSampleDatabaseArithmetic();
    ZiyuanCompare();
}



#define assert_agreement(X,Y) printf("compare "#X" =? "#Y"\n"); assert(std::abs(X-Y)<1e-12)

class DoubleLoader : public SampleLoader<double>
{
public:
    double LoadSample(int conf) const
    {
	char filename[512];
	sprintf(filename, "/home/gregm/fitting/NPR/G1_NPR/tests/jackknife_database/double/double.%d", conf);
	FILE *f = OpenFile(filename, "r");
	double ret;
	fscanf(f, "%le", &ret);
	fclose(f);
	printf("loaded %e from %s\n", ret, filename);
	return ret;
    }
};


void TestJackknifeDatabaseDouble()
{
    printf("start TestJackknifeDatabaseDouble\n");
    std::vector<int> double_confs{ 1, 2, 3 };
    DoubleLoader double_loader;
    ConfSampleDatabase<double> double_conf_db = ConfSampleDatabase<double>::Load(double_loader, double_confs);
    std::vector<std::vector<int> > jackknife_samples = MakeJackknifeSamples(double_confs.size(), 1);
    JackknifeDatabase<double> double_jack_db = double_conf_db.ComputeJackknifeMeans(jackknife_samples);
    double central = double_jack_db.CentralValue();
    double err = double_jack_db.Error();
    printf("calculated %e +/- %e\n", central, err);
    double true_central = 12.0;
    double true_err = 1.7320508075688772;
    printf("expected %e +/- %e\n", true_central, true_err);
    assert_agreement(central, true_central);
    assert_agreement(err, true_err);
    printf("Successfully finished TestJackknifeDatabaseDouble\n\n");
}


class ComplexLoader : public SampleLoader<std::complex<double> >
{
public:
    std::complex<double> LoadSample(int conf) const
    {
	char filename[512];
	sprintf(filename, "/home/gregm/fitting/NPR/G1_NPR/tests/jackknife_database/complex/complex.%d", conf);
	FILE *f = OpenFile(filename, "r");
	double re;
	double im;
	char sign;
	fscanf(f, "%le%c%lei ", &re, &sign, &im);
	if (sign == '-') {
	    im *= -1.0;
	} else if (sign != '+') {
	    printf("unknown sign %c\n", sign);
	    exit(-1);
	}
	fclose(f);
	std::complex<double> ret(re, im);
	printf("loaded %e+i%e from %s\n", ret.real(), ret.imag(), filename);
	return ret;
    }
};



void TestJackknifeDatabaseComplex()
{
    printf("start TestJackknifeDatabaseComplex\n");
    std::vector<int> confs{ 1, 2, 3 };
    ComplexLoader loader;
    ConfSampleDatabase<std::complex<double>> conf_db = ConfSampleDatabase<std::complex<double>>::Load(loader, confs);
    std::vector<std::vector<int> > jackknife_samples = MakeJackknifeSamples(confs.size(), 1);
    JackknifeDatabase<std::complex<double>> jack_db = conf_db.ComputeJackknifeMeans(jackknife_samples);
    std::complex<double> central = jack_db.CentralValue();
    std::complex<double> err = jack_db.Error();
    printf("calculated %e+i%e +/- %e+i%e\n", central.real(), central.imag(), err.real(), err.imag());
    std::complex<double> true_central(12.0, 0.0);
    std::complex<double> true_err(1.7320508075688772, 8.660254037844386);
    printf("calculated %e+i%e +/- %e+i%e\n", true_central.real(), true_central.imag(), true_err.real(), true_err.imag());
    assert_agreement(central.real(), true_central.real());
    assert_agreement(central.imag(), true_central.imag());
    assert_agreement(err.real(), true_err.real());
    assert_agreement(err.imag(), true_err.imag());
    printf("Successfully finished TestJackknifeDatabaseComplex\n\n");
}


void TestSampleDatabaseArithmetic()
{
    printf("start TestSampleDatabaseArithmetic\n");
    ConfSampleDatabase<std::complex<double>> a;
    a.Add(2);
    a.Add(3);
    ConfSampleDatabase<std::complex<double>> b;
    b.Add(10);
    b.Add(20);

    ConfSampleDatabase<std::complex<double>> sum = a + b;
    ConfSampleDatabase<std::complex<double>> difference = a - b;
    ConfSampleDatabase<std::complex<double>> five_times_a = 5 * a;
    assert_agreement(sum[0].real(), 12);
    assert_agreement(sum[1].real(), 23);
    assert_agreement(difference[0].real(), -8);
    assert_agreement(difference[1].real(), -17);
    assert_agreement(five_times_a[0].real(), 10);
    assert_agreement(five_times_a[1].real(), 15);

    a += b;
    assert_agreement(sum[0].real(), a[0].real());
    assert_agreement(sum[1].real(), a[1].real());

    a -= b;
    assert_agreement(a[0].real(), 2);
    assert_agreement(a[1].real(), 3);

    a *= 5;
    assert_agreement(a[0].real(), 10);
    assert_agreement(a[1].real(), 15);
    printf("Successfully finished TestSampleDatabaseArithmetic\n\n");
}



void ZiyuanCompare()
{
    printf("\n------------- ZiyuanCompare start ---------------\n\n");

#if 0
    const int mom1[] = { 1, 1, 2, 4 };
    const int mom2[] = { 2, 1, 2, -2 };

    const Parity parity = POSITIVE_PARITY;

    const char* dir = "tests/Ziyuan_compare";
    const int conf = 3620;

    WilsonMatrix prop1 = LoadPropagator(dir, conf, mom1);
    WilsonMatrix prop2 = LoadPropagator(dir, conf, mom2);

    array<array<DoubleWilsonMatrix, 4>, 6> unamputated_diagrams;
    for (int diagram = 0; diagram < 6; ++diagram) {
	unamputated_diagrams[diagram] = LoadLRFourQuarkDiagrams(dir, diagram_names[diagram], parity, conf, mom1, mom2);
    }


    array<array<DoubleWilsonMatrix, 4>, 6> amputated_diagrams;
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

    array<DoubleWilsonMatrix, 7> projector_spin_color_structures = BuildProjectorSpinColorStructures(parity);

    Matrix<std::complex<double>, 7> mixings = 
	ContractionsWithG1::DoContractionsFourQuarkOpFourQuarkExt(amputated_diagrams, projector_spin_color_structures);

    Matrix<std::complex<double>, 7> inv_Zfactors = mixings * tree_level_mixings.Inverse();
    inv_Zfactors.Print("inv_Zfactors");

    const Matrix<std::complex<double>, 7> Ziyuan_Zfactors = {{{
        {{{ 1.2050, -0.0000, -0.0000,  0.0000, 0.0000,  0.0031, -0.0003}}},
        {{{-0.0000,  0.9322,  0.2694, -0.1504, 0.2213, -0.0024, -0.0112}}},
        {{{-0.0000,  0.0216,  1.1962, -0.1252, 0.1232, -0.0010,  0.0035}}},
        {{{-0.0000, -1.1061,  0.0672,  0.6630, 1.0238,  0.0276, -0.0120}}},
        {{{-0.0000, -0.9581, -0.3291,  0.0667, 1.8359, -0.0126,  0.0256}}},
        {{{ 0.0005,  0.0184,  0.0193, -0.0034, 0.0243,  1.0658,  0.2256}}},
        {{{ 0.0002, -0.0266, -0.0630, -0.0094, 0.0001,  0.0877,  1.4526}}}
    }}};

    for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
            assert(abs(Ziyuan_Zfactors[i][j] - inv_Zfactors[i][j]) < 1e-4);
        }
    }
    printf("Ziyuan compare test PASSED!\n");

    printf("\n------------- ZiyuanCompare end ---------------\n\n");
#endif
}



