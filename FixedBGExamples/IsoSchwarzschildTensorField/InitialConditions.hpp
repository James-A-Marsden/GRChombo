/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALCONDITIONS_HPP_
#define INITIALCONDITIONS_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FixedBGTensorField.hpp"
#include "FourthOrderDerivatives.hpp"
#include "IsoSchwarzschildFixedBG.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include <boost/math/special_functions/bessel.hpp>
//! Class which creates the initial conditions
class InitialConditions
{
  protected:
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;
    const IsoSchwarzschildFixedBG::params_t m_bg_params;
    const FourthOrderDerivatives m_deriv;

    // Define the evolution vars
    template <class data_t>
    using Vars = FixedBGTensorField::template Vars<data_t>;
    // Now the fixed BG vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  public:
    //! The constructor for the class
    InitialConditions(const std::array<double, CH_SPACEDIM> a_center,
                      const IsoSchwarzschildFixedBG::params_t a_bg_params,
                      const double a_dx)
        : m_dx(a_dx), m_center(a_center), m_bg_params(a_bg_params),
          m_deriv(a_dx)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // get the metric vars
        IsoSchwarzschildFixedBG iso_bh(m_bg_params, m_dx);
        MetricVars<data_t> metric_vars;
        iso_bh.compute_metric_background(metric_vars, current_cell);
        const data_t det_gamma =
            TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
        const auto gamma_UU =
            TensorAlgebra::compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys =
            TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU);

        Vars<data_t> vars;
        VarsTools::assign(vars, 0.);
        Tensor<2, data_t> fspatial; // Spatial component of the tensor field
        Tensor<1, data_t> fbar;
        data_t fhat;
        Tensor<2, data_t> v; // Spatial rank 2 v field

        // Damping fields
        Tensor<1, data_t> Xspatial;
        data_t Xhat;

        // Set everything to zero
        vars.fhat = 0.0;
        vars.Xhat = 0.0;
        FOR1(i)
        {
            vars.fbar[i] = 0.0;
            vars.Xspatial[i] = 0.0;
        }

        FOR2(i, j)
        {
            vars.fspatial[i][j] = 0.0;
            vars.v[i][j] = 0.0;
        }
        double radius = coords.get_radius();
        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const data_t x2 = x * x;
        const double y2 = y * y;
        const double z2 = z * z;
        const data_t r = coords.get_radius();
        const data_t r2 = r * r;
        const data_t r3 = r2 * r;
        const data_t r4 = r2 * r2;
        const data_t r5 = r2 * r3;
        const data_t rho = coords.get_rho(); // simd_max(sqrt(x2 + y2), 1e-6);
        const data_t rho2 = rho * rho;

        const data_t costheta = z / r;
        const data_t sintheta = rho / r;

        const data_t cos2theta = costheta * costheta - sintheta * sintheta;
        const data_t sin2theta = 2.0 * sintheta * costheta;

        const data_t sinphi = y / rho;
        const data_t cosphi = x / rho;

        const data_t sin2phi = 2.0 * sinphi * cosphi;
        const data_t cos2phi = cosphi * cosphi - sinphi * sinphi;

        const double M = m_bg_params.mass;
        const double M2 = M * M;

        const double A = 1.0;
        vars.v[0][0] = -A * (2.0 * M + r) * sintheta * sintheta * costheta *
                       sinphi * cosphi / r;
        vars.v[1][1] = A * (2.0 * M + r) * sintheta * sintheta * costheta *
                       sinphi * cosphi / r;
        vars.v[0][1] = A * (2.0 * M + r) * sintheta * sintheta * costheta *
                       cos2phi / (2.0 * r);
        vars.v[0][2] = -A * ((2.0 * M + r) * costheta * costheta + 3.0 * r) *
                       sintheta * sinphi / (2.0 * r);
        vars.v[1][2] = A * ((2.0 * M + r) * costheta * costheta + 3.0 * r) *
                       sintheta * cosphi / (2.0 * r);
        vars.v[1][0] = vars.v[0][1];
        vars.v[2][0] = vars.v[0][2];
        vars.v[2][1] = vars.v[1][2];

        current_cell.store_vars(vars);
    }
};

#endif /* INITIALCONDITIONS_HPP_ */
