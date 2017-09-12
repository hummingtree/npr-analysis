#ifndef NPR_SETTINGS_H
#define NPR_SETTINGS_H

#include "Matrix.h"
#include "DoubleWilsonMatrix.h"
#include "NPR_Utils.h"
#include "util.h"
#include "G1Check.h"
#include <vector>
#include <array>

enum Scheme { SchemeGammaMu, SchemeQslash };

#define Nsub 3

struct NPRSettings 
{
    NPRSettings() {}

    const char* dir;
    const char* sub_dir;
    const char* c1_str;
    std::vector<int> confs;
    std::array<int, 4> mom1;
    std::array<int, 4> mom2;

// added by Jackson
	std::array<double, 4> cont_mom1;
	std::array<double, 4> cont_mom2;
	
    std::array<int, 4> Ls;
    Scheme scheme;
    Parity parity;
    bool do_disconnected;
    bool do_subtractions;
    bool enforce_reality;
    int bin_size;

    std::vector<std::vector<int>> jackknife_samples;
    int Njack;
    void SetJackknifeSamples() {
        const int Nbins = confs.size() / bin_size;
        jackknife_samples = MakeJackknifeSamples(Nbins, 1);
        Njack = jackknife_samples.size();
    }

    double p1[4], p2[4], q[4];
    double cont_p1[4], cont_p2[4], cont_q[4];
    double qsq;
    double cont_qsq;
    void SetMomenta() {
        double p1sq, p2sq;
        p1sq = p2sq = qsq = cont_qsq = 0;
        for (int mu = 0; mu < 4; ++mu) {
	   		p1[mu] = 2 * PI * mom1[mu] / Ls[mu];
	    	p2[mu] = 2 * PI * mom2[mu] / Ls[mu];
	    	cont_p1[mu] = 2 * PI * cont_mom1[mu] / Ls[mu];
	    	cont_p2[mu] = 2 * PI * cont_mom2[mu] / Ls[mu];
	    	q[mu] = p1[mu] - p2[mu];
	    	cont_q[mu] = cont_p1[mu] - cont_p2[mu];
        	p1sq += p1[mu] * p1[mu];
        	p2sq += p2[mu] * p2[mu];
	    	qsq += q[mu] * q[mu];
	    	cont_qsq += cont_q[mu] * cont_q[mu];
        }
        assert(abs(p1sq - p2sq) < 1e-12);
        assert(abs(p1sq - qsq) < 1e-12);
    }
    
    // these matrices are defined in NPRSettings.cpp
    static const Matrix<std::complex<double>, 8> tree_level_greens_funcs_8x8_gammamu;
    static const Matrix<std::complex<double>, 8> tree_level_greens_funcs_8x8_qslash;
    static const Matrix<std::complex<double>, 8> tree_level_greens_funcs_8x8_qslash_daiqian;
    std::array<DoubleWilsonMatrix, 7> projector_spin_color_structures; 
    DoubleWilsonMatrix HF_projector_spin_color_structure;
    Matrix<std::complex<double>, 8> tree_level_greens_funcs_8x8;
    void InitScheme() {
        if (scheme == SchemeGammaMu) {
	    projector_spin_color_structures = BuildProjectorSpinColorStructures(parity);
            tree_level_greens_funcs_8x8 = tree_level_greens_funcs_8x8_gammamu;
        } else { // scheme == SchemeQslash
	    projector_spin_color_structures = BuildQslashProjectorSpinColorStructures(q, qsq, parity);
            tree_level_greens_funcs_8x8 = tree_level_greens_funcs_8x8_qslash;
	    //projector_spin_color_structures = BuildDaiqianQslashProjectorSpinColorStructures(q, qsq, parity);
            //tree_level_greens_funcs_8x8 = tree_level_greens_funcs_8x8_qslash_daiqian;
            //printf("WARNING: USING DAIQIAN PROJECTOR SPIN COLOR STRUCTURES!\n");
        }
        HF_projector_spin_color_structure = BuildHFProjectorSpinColorStructure(parity);

        CheckG1TreeLevelProjectedGreensFunctions(*this);
    }

    std::array<const char*, Nsub> subtraction_op_names;
    std::array<SpinMatrix, Nsub> subtraction_projectors;
    void SetSubtractionOperators() {
        if (parity == POSITIVE_PARITY) {
	    subtraction_op_names = std::array<const char*, Nsub>{{ "scalar", "forward_covariant_dslash", "backward_covariant_dslash" }};
	    subtraction_projectors = std::array<SpinMatrix, Nsub>{{ SpinMatrix::One(), SpinMatrix::Slash(p2), SpinMatrix::Slash(p1) }};
        } else {
            subtraction_op_names = std::array<const char*, Nsub>{{ "pseudoscalar", "forward_covariant_dslash_g5", "backward_covariant_dslash_g5" }};
	    subtraction_projectors = std::array<SpinMatrix, Nsub>{{ SpinMatrix::Gamma5(), SpinMatrix::Slash(p2) * SpinMatrix::Gamma5(), SpinMatrix::Slash(p1) * SpinMatrix::Gamma5() }};
        }
    }

    void Show() {
        printf("--- NPR Settings: ---\n");
        printf("dir = %s\n", dir);
        printf("sub_dir = %s\n", sub_dir);
        printf("c1_str = %s\n", c1_str);
        printf("Nconf = %d\n", (int)confs.size());
        printf("mom1 = (%d, %d, %d, %d)\n", mom1[0], mom1[1], mom1[2], mom1[3]);
        printf("mom2 = (%d, %d, %d, %d)\n", mom2[0], mom2[1], mom2[2], mom2[3]);
        printf("cont_mom1 = (%.3f, %.3f, %.3f, %.3f)\n", cont_mom1[0], cont_mom1[1], cont_mom1[2], cont_mom1[3]);
        printf("cont_mom2 = (%.3f, %.3f, %.3f, %.3f)\n", cont_mom2[0], cont_mom2[1], cont_mom2[2], cont_mom2[3]);
        printf("Ls = (%d, %d, %d, %d)\n", Ls[0], Ls[1], Ls[2], Ls[3]);
        printf("p1 = (%f, %f, %f, %f)\n", p1[0], p1[1], p1[2], p1[3]);
        printf("p2 = (%f, %f, %f, %f)\n", p2[0], p2[1], p2[2], p2[3]);
        printf("q = (%f, %f, %f, %f)\n", q[0], q[1], q[2], q[3]);
        printf("qsq = %f\n", qsq);
        printf("scheme = %s\n", (scheme == SchemeGammaMu ? "GammaMu" : "Qslash"));
        printf("parity = %s\n", (parity == POSITIVE_PARITY ? "POSITIVE_PARITY" : "NEGATIVE_PARITY"));
        printf("do_disconnected = %d\n", do_disconnected);
        printf("do_subtractions = %d\n", do_subtractions);
        printf("enforce_reality = %d\n", enforce_reality);
        printf("bin_size = %d\n", bin_size);
        printf("Njack = %d\n", Njack);
        printf("Nsub = %d\n", Nsub);
        printf("subtraction operators are: ");
        for (int i = 0; i < Nsub; ++i) printf("%s, ", subtraction_op_names[i]);
        printf("\n");
        printf("---------------------\n");
    }

    void Init() {
        SetJackknifeSamples();
        SetMomenta();
		InitScheme();
        SetSubtractionOperators();
        Show();
    }
};

#endif
