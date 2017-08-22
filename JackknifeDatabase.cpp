#include "JackknifeDatabase.h"
#include "SpinMatrix.h"

#include <vector>

std::vector<std::vector<int> > MakeJackknifeSamples(int N, int block_size)
{
    std::vector<std::vector<int> > all_samples;

    std::vector<int> full_sample(N);
    for (int i = 0; i < N; ++i) full_sample[i] = i;
    all_samples.push_back(full_sample);

    if (N == 1) return all_samples;

    int num_blocks = N / block_size;
    for (int b = 0; b < num_blocks; ++b) {
	int delete_start = block_size * b;
	int delete_end = block_size * (b + 1);

	std::vector<int> sample;
	for (int i = 0; i < N; ++i) {
	    if (i < delete_start || i >= delete_end) sample.push_back(i);
	}
	all_samples.push_back(sample);
    }

    return all_samples;
}


void PrintDoubleWithError(const char *name, const JackknifeDatabase<double> &jack_db, bool newline)
{
    double central = jack_db.CentralValue();
    double err = jack_db.Error();
    printf("%s = %f     +/-    %f", name, central, err);
    if (newline) printf("\n");
}

void PrintComplexWithError(const char *name, const JackknifeDatabase<std::complex<double> > &jack_db, bool newline)
{
    std::complex<double> central = jack_db.CentralValue();
    std::complex<double> err = jack_db.Error();
    printf("%s = %f + i %f     +/-    %f + i %f", name, central.real(), central.imag(), err.real(), err.imag());
    if (newline) printf("\n");
}

void PrintSpinMatrixWithError(const char *name, const JackknifeDatabase<SpinMatrix> &jack_db)
{
    jack_db.CentralValue().data.Print(name);
    char err_name[512];
    sprintf(err_name, "error on %s", name);
    jack_db.Error().data.Print(err_name);
}
