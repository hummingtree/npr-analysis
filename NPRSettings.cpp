#include "NPRSettings.h"

// The Z factor matrix is given by
// Z / Zq ^ 2 = F inv(M)
// This is F
// Compare https://rbc.phys.columbia.edu/rbc_ukqcd/individual_postings/qliu/2011/k2pipiNPR/k2pipiNPR_2.pdf
const Matrix<std::complex<double>, 8> NPRSettings::tree_level_greens_funcs_8x8_gammamu = {{{
        {{{ 3072, 3072, 0, 0, 0, 0, 0, 0 }}},
        {{{ 537.6, -230.4, 1152, 0, 0, 0, 0, 0 }}},
        {{{ -230.4, 537.6, 384, 0, 0, 0, 0, 0 }}},
        {{{ 0, 0, 0, 1152, 384, 3456, 1152, 0 }}},
        {{{ 0, 0, 0, 384, 1152, 1152, 3456, 0 }}},
        {{{ 0, 0, 0, 1152, 384, 0, 0, 0 }}},
        {{{ 0, 0, 0, 384, 1152, 0, 0, 0 }}},
        {{{ 0, 768, 1536, 0, 768, 0, 2304, 768 }}} // this row gets checked in the code
}}};

const Matrix<std::complex<double>, 8> NPRSettings::tree_level_greens_funcs_8x8_qslash = {{{
        {{{ 1, 0, 0, 0, 0, 0, 0, 0 }}},
        {{{ 0, 1, 0, 0, 0, 0, 0, 0 }}},
        {{{ 0, 0, 1, 0, 0, 0, 0, 0 }}},
        {{{ 0, 0, 0, 1, 0, 0, 0, 0 }}},
        {{{ 0, 0, 0, 0, 1, 0, 0, 0 }}},
        {{{ 0, 0, 0, 0, 0, 1, 0, 0 }}},
        {{{ 0, 0, 0, 0, 0, 0, 1, 0 }}},
        {{{ 0, 0, 0, 0, 0, 0, 0, 768 }}} 
}}};

const Matrix<std::complex<double>, 8> NPRSettings::tree_level_greens_funcs_8x8_qslash_daiqian = {{{
        {{{ 768,     0,     0,   0,   0,   0,   0,   0 }}},
        {{{   0, 134.4, -57.6,   0,   0,   0,   0,   0 }}},
        {{{   0, -57.6, 134.4,   0,   0,   0,   0,   0 }}},
        {{{   0,     0,     0, 288,  96,   0,   0,   0 }}},
        {{{   0,     0,     0,  96, 288,   0,   0,   0 }}},
        {{{   0,     0,     0,   0,   0, 288,  96,   0 }}},
        {{{   0,     0,     0,   0,   0,  96, 288,   0 }}},
        {{{   0,     0,     0,   0,   0,   0,   0, 768 }}}
}}};

const Matrix<std::complex<double>, 5> NPRSettings::tree_level_greens_funcs_DSeq2_5x5_GammaMu = {{{
        {{{ 3072, 0, 0, 0, 0	}}},
        {{{ 0, 2304, -384, 0, 0	}}},
        {{{ 0, -384, 576, 0, 0	}}},
        {{{ 0, 0, 0, 480, 288	}}},
        {{{ 0, 0, 0, 288, 2016	}}},
}}};

const Matrix<std::complex<double>, 5> NPRSettings::tree_level_greens_funcs_DSeq2_5x5_Qslash = {{{
        {{{ 768, 0, 0, 0, 0		}}},
        {{{ 0, 576, 192, 0, 0	}}},
        {{{ 0, -96, -288, 0, 0	}}},
        {{{ 0, 0, 0, 432, 144 	}}},
        {{{ 0, 0, 0, 720, 1008	}}},
}}};
