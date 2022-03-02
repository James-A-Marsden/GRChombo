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
    c_trace_field,
    c_traceF1,
    c_traceF2,
    c_traceF3,
    c_traceB,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "chi", "rho", "rhoJ", "Edot", "Jdot", "Kout", "d1Kout1", "d1Kout2", "d1Kout3",
    "trace_field",
    "traceF1", "traceF2", "traceF3",
    "traceB"
    };
/*
    
*/
}

#endif /* DIAGNOSTICVARIABLES_HPP */
