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
        // derivs
        const auto vars = current_cell.template load_vars<MatterVars>();

        // get the metric vars
        MetricVars<data_t> metric_vars;
        m_background.compute_metric_background(metric_vars, current_cell);

        using namespace TensorAlgebra;
        const auto gamma_UU =
            TensorAlgebra::compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys =
            compute_christoffel(metric_vars.d1_gamma, gamma_UU);

        data_t fspatial_trace =
            TensorAlgebra::compute_trace(vars.fspatial, gamma_UU);

        Tensor<2, data_t> new_fspatial;

        FOR2(i, j)
        {
            new_fspatial[i][j] =
                vars.fspatial[i][j] -
                (1.0 / 3.0) * metric_vars.gamma[i][j] *
                    (fspatial_trace + vars.fhat / metric_vars.lapse);
        }
        /*
        current_cell.store_vars(new_fspatial[0][0], c_fspatial11);
        current_cell.store_vars(new_fspatial[0][1], c_fspatial12);
        current_cell.store_vars(new_fspatial[0][2], c_fspatial13);
        current_cell.store_vars(new_fspatial[1][1], c_fspatial22);
        current_cell.store_vars(new_fspatial[1][2], c_fspatial23);
        current_cell.store_vars(new_fspatial[2][2], c_fspatial33);
        */
    }
};
#endif /* TRACEFIELDREMOVAL_HPP_ */
