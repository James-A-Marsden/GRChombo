/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "DiagnosticVariables.hpp"
#include <array>
#include <string>

// assign an enum to each variable
enum
{
    //Spatial components of tensor field
    c_fspatial11, 
    c_fspatial12,
    c_fspatial13,
    c_fspatial22,
    c_fspatial23,
    c_fspatial33,

    //Rank 2 component conjugate field
    c_v11,
    c_v12,
    c_v13,
    c_v22,
    c_v23,
    c_v33,

    c_Xhat,
    c_Xspatial1,
    c_Xspatial2,
    c_Xspatial3,

    NUM_VARS
};

namespace UserVariables
{
    //for checkpoints 
static const std::array<std::string, NUM_VARS> variable_names = {
    "fspatial11", "fspatial12", "fspatial13", "fspatial22", "fspatial23", "fspatial33",
    "v11", "v12", "v13", "v22", "v23", "v33",
    "Xhat", "Xspatial1", "Xspatial2", "Xspatial3"
    };
}

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
