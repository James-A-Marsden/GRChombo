/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include <array>
#include <string>

enum
{
    c_d1_gamma_UU11,
    c_d1_gamma_UU12,
    c_d1_gamma_UU13,
    c_d1_gamma_UU22,
    c_d1_gamma_UU23,
    c_d1_gamma_UU33,

    c_d2_gamma11,
    c_d2_gamma12,
    c_d2_gamma13,
    c_d2_gamma22,
    c_d2_gamma23,
    c_d2_gamma33,

    c_d1_K_tensor11,
    c_d1_K_tensor12,
    c_d1_K_tensor13,
    c_d1_K_tensor22,
    c_d1_K_tensor23,
    c_d1_K_tensor33,

    c_d2_lapse, 

    c_d1_shift1,
    c_d1_shift2,
    c_d1_shift3,

    c_d2_shift1,
    c_d2_shift2,
    c_d2_shift3,

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS> variable_names = {
    "d1_gammaUU11",
    "d1_gammaUU12",
    "d1_gammaUU13",
    "d1_gammaUU22",
    "d1_gammaUU23",
    "d1_gammaUU33",

    "d2_gamma11",
    "d2_gamma12",
    "d2_gamma13",
    "d2_gamma22",
    "d2_gamma23",
    "d2_gamma33",

    "d1_K_tensor11",
    "d1_K_tensor12",
    "d1_K_tensor13",
    "d1_K_tensor22",
    "d1_K_tensor23",
    "d1_K_tensor33",

    "d2_lapse",

    "d1_shift1",
    "d1_shift2",
    "d1_shift3",

    "d2_shift1",
    "d2_shift2",
    "d2_shift3",
    };
}

#endif /* USERVARIABLES_HPP */
