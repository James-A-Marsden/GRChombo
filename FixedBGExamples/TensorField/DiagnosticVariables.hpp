/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_chi,    // the conformal factor of the metric, not evolved
    c_rho,    // the energy density of the SF
    c_rhoJ,   // the energy density of the SF
    c_Edot,   // the energy density of the SF
    c_Jdot,   // the energy density of the SF
    c_Kout,
    c_d1Kout1,
    c_d1Kout2,
    c_d1Kout3,
    /*
    c_d2_shift111,
    c_d2_shift112,
    c_d2_shift113,
    c_d2_shift121,
    c_d2_shift122,
    c_d2_shift123,
    c_d2_shift131,
    c_d2_shift132,
    c_d2_shift133,
    c_d1_K_tensor111,
    c_d1_K_tensor112,
    c_d1_K_tensor113,
    c_d1_K_tensor122,
    c_d1_K_tensor123,
    c_d1_K_tensor133,
    c_d1_K_tensor211,
    c_d1_K_tensor212,
    c_d1_K_tensor213,
    c_d1_K_tensor222,
    c_d1_K_tensor223,
    c_d1_K_tensor233,
    c_d1_K_tensor311,
    c_d1_K_tensor312,
    c_d1_K_tensor313,
    c_d1_K_tensor322,
    c_d1_K_tensor323,
    c_d1_K_tensor333,
    c_d1_gamma_UU111,
    c_d1_gamma_UU121,
    c_d1_gamma_UU131,
    c_d1_gamma_UU221,
    c_d1_gamma_UU231,
    c_d1_gamma_UU331,
    c_d1_gamma_UU112,
    c_d1_gamma_UU122,
    c_d1_gamma_UU132,
    c_d1_gamma_UU222,
    c_d1_gamma_UU232,
    c_d1_gamma_UU332,
    c_d1_gamma_UU113,
    c_d1_gamma_UU123,
    c_d1_gamma_UU133,
    c_d1_gamma_UU223,
    c_d1_gamma_UU233,
    c_d1_gamma_UU333,
    c_ricci_phys11,
    c_ricci_phys12,
    c_ricci_phys13,
    c_ricci_phys22,
    c_ricci_phys23,
    c_ricci_phys33,
    */
    /*
    c_d1_K_tensor111,
    c_d1_K_tensor121,
    c_d1_K_tensor131,
    c_d1_K_tensor221,
    c_d1_K_tensor231,
    c_d1_K_tensor331,
    c_d1_K_tensor112,
    c_d1_K_tensor122,
    c_d1_K_tensor132,
    c_d1_K_tensor222,
    c_d1_K_tensor232,
    c_d1_K_tensor332,
    c_d1_K_tensor113,
    c_d1_K_tensor123,
    c_d1_K_tensor133,
    c_d1_K_tensor223,
    c_d1_K_tensor233,
    c_d1_K_tensor333,
    */
    c_traceF,
    c_traceB,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "chi", "rho", "rhoJ", "Edot", "Jdot", "Kout", "d1Kout1", "d1Kout2", "d1Kout3",
    //"ricci_phys11", "ricci_phys12", "ricci_phys13", "ricci_phys22", "ricci_phys23", "ricci_phys33"
    /*
    "d2_shift111", "d2_shift112", "d2_shift113", "d2_shift121", "d2_shift122", "d2_shift123", "d2_shift131","d2_shift132","d2_shift133",
    */
    /*
    "d1_gamma_UU111", "d1_gamma_UU121", "d1_gamma_UU131", "d1_gamma_UU221", "d1_gamma_UU231", "d1_gamma_UU331",
    "d1_gamma_UU112", "d1_gamma_UU122", "d1_gamma_UU132", "d1_gamma_UU222", "d1_gamma_UU232", "d1_gamma_UU332",
    "d1_gamma_UU113", "d1_gamma_UU123", "d1_gamma_UU133", "d1_gamma_UU223", "d1_gamma_UU233", "d1_gamma_UU333"
    */
    /*
    "d1_K_tensor111", "d1_K_tensor121", "d1_K_tensor131", "d1_K_tensor221", "d1_K_tensor231", "d1_K_tensor331",
    "d1_K_tensor112", "d1_K_tensor122", "d1_K_tensor132", "d1_K_tensor222", "d1_K_tensor232", "d1_K_tensor332",
    "d1_K_tensor113", "d1_K_tensor123", "d1_K_tensor133", "d1_K_tensor223", "d1_K_tensor233", "d1_K_tensor333",
    */
    "traceF",
    "traceB"
    };
/*
    
*/
}

#endif /* DIAGNOSTICVARIABLES_HPP */
