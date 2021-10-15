/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALCONDITIONS_HPP_
#define INITIALCONDITIONS_HPP_

#include "ADMFixedBGVars.hpp"
#include "BoostedIsotropicBHFixedBG.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates the initial conditions

template <class background_t> class InitialConditions
{
  protected:
    const double m_dx;
    const double m_amplitude;
    const double m_omega;
    const std::array<double, CH_SPACEDIM> m_center;
    const background_t m_background;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  public:
    //! The constructor for the class
    InitialConditions(const double a_amplitude, const double a_omega,
                      const std::array<double, CH_SPACEDIM> a_center,
                      const background_t a_background, const double a_dx)
        : m_dx(a_dx), m_amplitude(a_amplitude), m_center(a_center),
          m_omega(a_omega), m_background(a_background)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_center);
        const data_t r = coords.get_radius();
        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;

        // get the metric vars for the background metric
        MetricVars<data_t> metric_vars;
        m_background.compute_metric_background(metric_vars, current_cell);
        const data_t det_gamma =
            TensorAlgebra::compute_determinant_sym(metric_vars.gamma);

        data_t phi = m_amplitude * exp(-r * r);
        data_t Pi = m_amplitude * m_omega * exp(-r * r);

        // Store the initial values of the variables
        current_cell.store_vars(phi, c_phi);
        current_cell.store_vars(Pi, c_Pi);
    }
};

#endif /* INITIALCONDITIONS_HPP_ */
