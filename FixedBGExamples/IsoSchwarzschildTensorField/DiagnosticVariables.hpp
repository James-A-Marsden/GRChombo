/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_chi,  // the conformal factor of the metric, not evolved
    c_rho,  // the energy density of the SF
    c_rhoJ, // the energy density of the SF
    c_Edot, // the energy density of the SF
    c_Jdot, // the energy density of the SF
    c_trace_field,
    c_trace_momentum,
    c_transverseScalar,
    c_transverseVector1,
    c_transverseVector2,
    c_transverseVector3,
    c_primaryConstraintScalar,
    c_primaryConstraintVector1,
    c_primaryConstraintVector2,
    c_primaryConstraintVector3,
    c_rho_eff,
    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "chi",
    "rho",
    "rhoJ",
    "Edot",
    "Jdot",
    "trace_field",
    "trace_momentum",
    "transverseScalar",
    "transverseVector1",
    "transverseVector2",
    "transverseVector3",
    "primaryConstraintScalar",
    "primaryConstraintVector1",
    "primaryConstraintVector2",
    "primaryConstraintVector3",
    "rho_eff",
};

}

#endif /* DIAGNOSTICVARIABLES_HPP */
