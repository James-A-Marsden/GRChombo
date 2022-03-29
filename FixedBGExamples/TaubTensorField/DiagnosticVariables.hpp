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
    c_traceB,
    c_traceRT,
    c_d1fhat1,
    c_d1fhat2,
    c_d1fhat3,
    c_KinvK,
    c_lorenzCont,
    c_lorenzProj1,
    c_lorenzProj2,
    c_lorenzProj3,
    c_primaryConstraintScalar,
    c_primaryConstraintVector1,
    c_primaryConstraintVector2,
    c_primaryConstraintVector3,
    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "chi", "rho", "rhoJ", "Edot", "Jdot", "Kout", "d1Kout1", "d1Kout2", "d1Kout3",
    "trace_field",
    "traceB",
    "traceRT",
    "d1fhat1", "d1fhat2", "d1fhat3",
    "KinvK",
    "lorenzCont", "lorenzProj1", "lorenzProj2", "lorenzProj3",
    "primaryConstraintScalar",
    "primaryConstraintVector1", "primaryConstraintVector2", "primaryConstraintVector3"
    };
/*
    
*/
}

#endif /* DIAGNOSTICVARIABLES_HPP */
