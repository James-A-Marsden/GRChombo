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

    //Rank 2 component conjugate field
    c_v11,
    c_v12,
    c_v13,
    c_v22,
    c_v23,
    c_v33,

    //Rank 1 component conjugate field 
    c_q1,
    c_q2,
    c_q3,

    //Scalar component conjugate field 
    c_w,

        //Damping variables
    c_X1,
    c_X2,
    c_X3,
    c_theta,
    
    NUM_VARS
};

namespace UserVariables
{
    //for checkpoints 
static const std::array<std::string, NUM_VARS> variable_names = {
    "fspatial11", "fspatial12", "fspatial13", "fspatial22", "fspatial23", "fspatial33", 
    "fbar1", "fbar2", "fbar3", 
    "fhat", 
    "v11", "v12", "v13", "v22", "v23", "v33", 
    "q1", "q2", "q3",
    "w",
    "X1","X2","X3",
    "theta"};
}

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
