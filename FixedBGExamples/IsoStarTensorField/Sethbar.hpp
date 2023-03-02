/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SETHBAR_HPP_
#define SETHBAR_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FixedBGTensorField.hpp"
#include "FourthOrderDerivatives.hpp"
#include "IsoStarFixedBG.hpp"
#include "Tensor.hpp"
#include "TensorPotential.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include <boost/math/special_functions/bessel.hpp>
//! Class which creates the initial constraints
class Sethbar
{
  protected:
    const double m_dx;
    // const double m_amplitude_re, m_amplitude_im;
    // const double m_omega;
    const std::array<double, CH_SPACEDIM> m_center;
    const IsoStarFixedBG::params_t m_bg_params;
    const double m_tensor_mass;
    const double m_initial_constant;
    const FourthOrderDerivatives m_deriv;

    // load in Vars from the field
    //  The evolution vars

    template <class data_t>
    using Vars = FixedBGTensorField<TensorPotential>::template Vars<data_t>;
    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  public:
    //! The constructor for the class
    // const double a_amplitude_re, const double a_amplitude_im, const double
    // a_omega,
    Sethbar(const double tensor_mass,
            const std::array<double, CH_SPACEDIM> a_center,
            const IsoStarFixedBG::params_t a_bg_params, const double a_dx,
            const double a_initial_constant) //, const double a_fhat, const
                                             //Tensor<1,data_t> a_fbar, const
                                             //Tensor<2,data_t> a_fspatial)
        : m_dx(a_dx), m_center(a_center), m_bg_params(a_bg_params),
          m_tensor_mass(tensor_mass), m_initial_constant(a_initial_constant),
          m_deriv(a_dx)
    // m_fhat(a_fhat), m_fbar(a_fbar), m_fspatial(a_fspatial)
    //, m_amplitude_re(a_amplitude_re),
    // m_amplitude_im(a_amplitude_im)
    // m_omega(a_omega),
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // get the metric vars
        IsoStarFixedBG kerr_bh(m_bg_params, m_dx);
        MetricVars<data_t> metric_vars;
        kerr_bh.compute_metric_background(metric_vars, current_cell);
        const data_t det_gamma =
            TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
        const auto gamma_UU =
            TensorAlgebra::compute_inverse_sym(metric_vars.gamma);

        Vars<data_t> vars;

        Tensor<3, data_t> chris_local;
        FOR3(i, j, k)
        {
            chris_local[i][j][k] = 0.0;
            FOR1(l)
            {
                chris_local[i][j][k] += 0.5 * gamma_UU[i][l] *
                                        (metric_vars.d1_gamma[l][k][j] +
                                         metric_vars.d1_gamma[j][l][k] -
                                         metric_vars.d1_gamma[j][k][l]);
            }
        }

        Tensor<3, data_t> cd1_K_tensor;
        FOR3(i, j, k)
        {
            cd1_K_tensor[i][j][k] = metric_vars.d1_K_tensor[i][j][k];
            FOR1(l)
            {
                cd1_K_tensor[i][j][k] +=
                    -chris_local[l][k][i] * metric_vars.K_tensor[l][j] -
                    chris_local[l][k][j] * metric_vars.K_tensor[i][l];
            }
        }

        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const data_t x2 = x * x;
        const double y2 = y * y;
        const double z2 = z * z;
        const data_t r = coords.get_radius();
        const data_t r2 = r * r;
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

        const auto local_vars = current_cell.template load_vars<Vars>();

        // Calculate the derivatives

        vars.fhat = local_vars.fhat;
        FOR1(i)
        {
            vars.fbar[i] = local_vars.fbar[i];
            FOR1(j) { vars.fspatial[i][j] = local_vars.fspatial[i][j]; }
        }

        const auto d1 = m_deriv.template diff1<Vars>(current_cell);

        data_t A = pow(r, 4.0) * pow(0.5 * M + r, -10.0);
        vars.fbar[0] = cosphi * sintheta * A;
        vars.fbar[1] = sinphi * sintheta * A;
        vars.fbar[2] = costheta * A;

        current_cell.store_vars(vars);
    }
};

#endif /* SETHBAR_HPP_ */
