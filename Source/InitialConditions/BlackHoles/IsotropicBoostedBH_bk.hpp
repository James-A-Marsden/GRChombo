/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ISOTROPICBOOSTEDBH_HPP_
#define ISOTROPICBOOSTEDBH_HPP_

#include "CCZ4Vars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "parstream.H" //gives pout
#include "simd.hpp"
#include <array>
#include <vector>

//! Class which solves for the initial data for a spherically symmetric boson
//! star with phi^4 coupling
class IsotropicBoostedBH
{

  public:
    //! The constructor
    IsotropicBoostedBH(BoostedBH::params_t a_bh1_params,
                       BoostedBH::params_t a_bh2_params, double a_dx)
        : m_bh1_params(a_bh1_params), m_bh2_params(a_bh2_params), m_dx(a_dx)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    BoostedBH::params_t m_bh1_params;
    BoostedBH::params_t m_bh2_params;
    double m_dx;
};

#include "IsotropicBoostedBH_bk.impl.hpp"

#endif /* ISOTROPICBOOSTEDBH_HPP_ */
