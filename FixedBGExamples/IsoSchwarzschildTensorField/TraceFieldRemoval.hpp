/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This class enforces A to be trace-free

#ifndef TRACEFIELDREMOVAL_HPP_
#define TRACEFIELDREMOVAL_HPP_

#include "ADMFixedBGVars.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Calculates the energy density rho and angular momentum density rhoJ
//! with matter type matter_t and writes it to the grid
template <class matter_t, class background_t> class TraceFieldRemoval
{
    // Use the variable definition in the matter class
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
                 //!< The matter object
    const background_t m_background; //!< The matter object

  public:
    TraceFieldRemoval(background_t a_background, const double a_dx)
        : m_deriv(a_dx), m_background(a_background)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables, and calc 1st
        const auto vars = current_cell.template load_vars<MatterVars>();

        // get the metric vars
        MetricVars<data_t> metric_vars;
        m_background.compute_metric_background(metric_vars, current_cell);

        using namespace TensorAlgebra;
        const auto gamma_UU =
            TensorAlgebra::compute_inverse_sym(metric_vars.gamma);

        data_t fspatial_trace =
            TensorAlgebra::compute_trace(vars.fspatial, gamma_UU);

        current_cell.store_vars(fspatial_trace, c_fhat);

    }
};
#endif /* TRACEFIELDREMOVAL_HPP_ */
