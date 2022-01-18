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
    //c_phi_Re, // matter field added
    //c_Pi_Re,  //(minus) conjugate momentum
    //c_phi_Im, // matter field added
    //c_Pi_Im,  //(minus) conjugate momentum

    //check public - define components 

    //Spatial components of tensor field  //Don't need all as this is symmetric? 
    c_fspatial11, 
    c_fspatial12,
    c_fspatial13,
    c_fspatial22,
    c_fspatial23,
    c_fspatial33,

    //Vector component of tensor field
    c_fbar1, 
    c_fbar2,
    c_fbar3,

    //Scalar component of tensor field
    c_fhat, 

    //First Rank 2 component conjugate field 
    //Recall that this is not symmetric!
    c_u11,
    c_u12,
    c_u13,
    c_u21,
    c_u22,
    c_u23,
    c_u31,
    c_u32,
    c_u33, 

    //Second Rank 2 component conjugate field
    c_v11,
    c_v12,
    c_v13,
    c_v22,
    c_v23,
    c_v33,
    
    //Rank 1 component conjugate field
    c_p1,
    c_p2,
    c_p3,

    //Rank 1 component conjugate field 
    c_q1,
    c_q2,
    c_q3,

    //Scalar component conjugate field 
    c_w, 

    NUM_VARS
};

namespace UserVariables
{
    //for checkpoints 
static const std::array<std::string, NUM_VARS> variable_names = {
    "fspatial11", "fspatial12", "fspatial13", "fspatial22", "fspatial23", "fspatial33", 
    "fbar1", "fbar2", "fbar3", 
    "fhat", 
    "u11", "u12", "u13", "u21", "u22", "u23", "u31", "u32", "u33", 
    "v11", "v12", "v13", "v22", "v23", "v33", 
    "p1", "p2", "p3", 
    "q1", "q2", "q3" 
    "w"};
    //"phi_Re", "Pi_Re", "phi_Im", "Pi_Im"};
}

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
