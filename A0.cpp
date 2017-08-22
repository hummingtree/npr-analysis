#include "JackknifeDatabase.h"
#include "Matrix.h"
#include <complex>
using namespace std;


// Use Wilson coefficients and RI->MS bar conversoin at 2.28 GeV
static array<complex<double>, 10> ComputeA0_2p28GeV(const Matrix<double, 7> &R_lat_to_RI)
{
    array<array<double, 7>, 10> Z_RI_to_MS_2p28_GeV {{
        {{ 0.20094, 0.99765, 0.0070425, 0, 0, 0, 0 }}, 
        {{ 0.20094, 0.0021174, 0.98616, 0.0016416, -0.0049249, 0, 0 }}, 
        {{ 0, 2.9972, 1.9934, 0.0032831, -0.0098498, 0, 0 }}, 
        {{ 0, 1.9803, 2.9233, 0.012754, -0.038093, 0, 0 }}, 
        {{ 0, 0, 0, 1.0010, -0.0028716, 0, 0 }}, 
        {{ 0, -0.036937, -0.086187, -0.023451, 1.0627, 0, 0 }}, 
        {{ 0, 0, 0, 0, 0, 1.0010, -0.0028716 }}, 
        {{ 0, 0, 0, 0, 0, -0.035764, 1.0996 }}, 
        {{ 0.30141, -0.0021174, -0.98616, -0.0016416, 0.0049249, 0, 0 }}, 
        {{ 0.30141, -0.99765, -0.0070425, 0, 0, 0, 0 }}
    }};

    array<double, 7> lattice_A0_matrix_elements {{-0.147, -0.218, 0.295, -0.601, -1.188, 1.33, 4.65 }};

    array<double, 7> RI_A0_matrix_elements = R_lat_to_RI * lattice_A0_matrix_elements;

    array<double, 10> MS_A0_matrix_elements;
    for (int i = 0; i < 10; ++i) {
        MS_A0_matrix_elements[i] = 0;
        for (int j = 0; j < 7; ++j) {
            MS_A0_matrix_elements[i] += Z_RI_to_MS_2p28_GeV[i][j] * RI_A0_matrix_elements[j];
        }
    }

    // All constants from Ziyuan
    double Gf = 1.16637e-5;
    double Vus = 0.2253;
    double Vud = 0.97425;

    array<double, 10> z {{ -0.2879390507, 1.138420148, -0.002847074467, 0.01169943305, -0.001628070946, 
                           0.007936074279, 7.371992154e-05, -9.489150846e-05, 1.513586959e-05, 6.97390239e-05 }};
    array<double, 10> y {{ 0.0, 0.0, 0.02344132124, -0.05789935776, 0.01059643117, -0.06829450251,
                         -0.0003202699843, 0.0006886695743, -0.00994958219, 0.002664737778 }};

    complex<double> tau(0.001543, -0.000635);

    // Full complex Wilson coefficients
    array<complex<double>, 10> C;
    for (int i = 0; i < 10; ++i) {
        C[i] = z[i] + tau * y[i];
    }

    array<complex<double>, 10> A0_amplitudes;
    for (int i = 0; i < 10; ++i) {
        A0_amplitudes[i] = C[i] * MS_A0_matrix_elements[i] * (Gf/sqrt(2.0)) * Vus * Vud;
    }

    return A0_amplitudes;
}

static complex<double> SumA0Amplitudes(const array<complex<double>, 10> &op_amplitudes)
{
    complex<double> total_amplitude(0);
    for (int i = 0; i < 10; ++i) {
        total_amplitude += op_amplitudes[i];
    }
    return total_amplitude;
}

void PrintAmplitude(double value, double error)
{
    if (value < 0) {
        printf("-");
        value = -value;
    }
    if (value == 0) {
        printf("0");
        return;
    }
    int exponent = floor(log10(value));
    double scaled_value = value * pow(10.0, -exponent);
    int error_digits = (int)round(error * pow(10.0, -exponent + 3));
    printf("%0.3f(%d) \\times 10^{%d}", scaled_value, error_digits, exponent);
}

void DoA0_2p28GeV(const JackknifeDatabase<Matrix<complex<double>, 7>> &jack_R_lat_to_RI)
{
    array<JackknifeDatabase<complex<double>>, 10> jack_op_amplitudes;
    JackknifeDatabase<complex<double>> jack_A0;

    for (int jack = 0; jack < jack_R_lat_to_RI.Size(); ++jack) {
        Matrix<double, 7> R_lat_to_RI = RealMatrix(jack_R_lat_to_RI[jack]);
        array<complex<double>, 10> op_amplitudes = ComputeA0_2p28GeV(R_lat_to_RI);
        for (int i = 0; i < 10; ++i) {
            jack_op_amplitudes[i].Add(op_amplitudes[i]);
        }
        complex<double> A0 = SumA0Amplitudes(op_amplitudes);
        jack_A0.Add(A0);
    }

    printf("A0 amplitudes, using Wilson coefficients and RI->MS coefficients at mu=2.28 GeV\n");
    printf("\n");
    complex<double> A0_central = jack_A0.CentralValue();
    complex<double> A0_err = jack_A0.Error();
    printf("Re A0 = %e +/- %e = ", A0_central.real(), A0_err.real());
    PrintAmplitude(A0_central.real(), A0_err.real());
    printf("\n");
    printf("Im A0 = %e +/- %e = ", A0_central.imag(), A0_err.imag());
    PrintAmplitude(A0_central.imag(), A0_err.imag());
    printf("\n");

    printf("\n");
    for (int i = 0; i < 10; ++i) {
        complex<double> amp_central = jack_op_amplitudes[i].CentralValue();
        complex<double> amp_err = jack_op_amplitudes[i].Error();
        printf("contribution from Q_%d: Real part: %e +/- %e    Imag part: %e +/- %e\n", 
                i+1, amp_central.real(), amp_err.real(), amp_central.imag(), amp_err.imag());
    }
    printf("Latex table:\n");
    for (int i = 0; i < 10; ++i) {
        complex<double> amp_central = jack_op_amplitudes[i].CentralValue();
        complex<double> amp_err = jack_op_amplitudes[i].Error();
        printf("%d & $", i+1);
        PrintAmplitude(amp_central.real(), amp_err.real());
        printf("$ & $");
        PrintAmplitude(amp_central.imag(), amp_err.imag());
        printf("$ \\\\\n");
    }
}


void Test_DoA0()
{
    Matrix<double, 7> fake_R;
    fake_R.Identity();
    array<complex<double>, 10> amps = ComputeA0_2p28GeV(fake_R);
    complex<double> A0 = SumA0Amplitudes(amps);
    printf("A0 test with fake data");
    printf("\n");
    printf("A0 = %e + i %e", A0.real(), A0.imag());
    printf("\n");
    for (int i = 0; i < 10; ++i) {
        printf("contribution from Q_%d: %e + i %e\n", i+1, amps[i].real(), amps[i].imag());
    }

    double value = -1.234567e-6;
    double err = 9.87e-8;
    printf("format test: %e +/- %e = ", value, err);
    PrintAmplitude(value, err);
    printf("\n");
}




