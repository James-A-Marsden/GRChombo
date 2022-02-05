/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ADMFIXEDBGVARS_HPP_
#define ADMFIXEDBGVARS_HPP_

#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"

/// Namespace for ADM vars for fixed BG evolution
namespace ADMFixedBGVars
{
/// Vars object for ADM vars used in FixedBG evolution
template <class data_t> struct Vars
{
    // ADM vars needed in matter only rhs (ok for Proca and SF)
    Tensor<2, data_t> gamma;
    Tensor<2, data_t> K_tensor;
    Tensor<2, Tensor<1,data_t>> d1_K_tensor; 
    data_t K;
    Tensor<1, data_t> d1_K;
    data_t lapse;
    Tensor<1, data_t> shift;
    Tensor<2, Tensor<1, data_t>> d1_gamma;
    Tensor<2, Tensor<2, data_t>> d2_gamma;
    Tensor<2, Tensor<1, data_t>> d1_gamma_UU;
    Tensor<3, Tensor<1, data_t>> d1_chris_phys; 
    Tensor<4, data_t> riemann_phys_ULLL; 
    Tensor<1, data_t> d1_lapse;
    Tensor<2, data_t> d1_shift;
    //Extra derivatives needed for tensor field case
    Tensor<2, data_t> d2_lapse;
    Tensor<3, data_t> d2_shift;
};

} // namespace ADMFixedBGVars

#endif /* ADMFIXEDBGVARS_HPP_ */
