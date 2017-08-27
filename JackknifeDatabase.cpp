#include "JackknifeDatabase.h"
#include "SpinMatrix.h"

#include <vector>

static std::string format_error(double value, double error);

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
	std::string output = format_error(central, err);
	printf("%s = %s", name, output.c_str());
    if (newline) printf("\n");
}

void PrintComplexWithError(const char *name, const JackknifeDatabase<std::complex<double> > &jack_db, bool newline)
{
//	std::complex<double> central = jack_db.CentralValue();
    std::complex<double> central = central_value_no_bias(jack_db.samples);
    std::complex<double> err = jack_db.Error();
	std::string output_real = format_error(central.real(), err.real());
	std::string output_imag = format_error(central.imag(), err.imag());
    printf("%s = %s + i %s", name, output_real.c_str(), output_imag.c_str());
    if (newline) printf("\n");
}

void PrintSpinMatrixWithError(const char *name, const JackknifeDatabase<SpinMatrix> &jack_db)
{
    jack_db.CentralValue().data.Print(name);
    char err_name[512];
    sprintf(err_name, "error on %s", name);
    jack_db.Error().data.Print(err_name);
}

static std::string format_error(double value, double error){
    assert(error > 0.);
    char output[512];
    if(error > 10.){
        sprintf(output, "%d(%d)", int(value), int(error));
        return std::string(output);
    }else if(error > 1.){
        sprintf(output, "%0.1f(%0.1f)", value, error);
    }else{
        int d = int(2.-1E-12-log(error)/log(10.));
        char format[512];
        sprintf(format, "%%+0.%df(%%d)", d);
        sprintf(output, format, value, int(error*pow(10.,d)));
    }
    return std::string(output);
}
