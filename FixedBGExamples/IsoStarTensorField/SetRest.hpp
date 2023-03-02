/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SETREST_HPP_
#define SETREST_HPP_

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
class SetRest
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
    SetRest(const double tensor_mass,
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
        const auto chris_phys =
            TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU);

        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const data_t x2 = x * x;
        const double y2 = y * y;
        const double z2 = z * z;
        const data_t r = coords.get_radius();
        const data_t r2 = r * r;
        const double M = m_bg_params.mass;
        const data_t rho = simd_max(sqrt(x2 + y2), 1e-6);
        const data_t rho2 = rho * rho;

        const data_t costheta = z / r;
        const data_t sintheta = rho / r;

        const data_t sinphi = coords.y / rho;
        const data_t cosphi = coords.x / rho;

        const data_t sin2phi = 2.0 * sinphi * cosphi;
        const data_t cos2phi = cosphi * cosphi - sinphi * sinphi;

        Vars<data_t> vars;

        const auto local_vars = current_cell.template load_vars<Vars>();

        vars.fhat = local_vars.fhat;
        vars.w = local_vars.w;
        FOR1(i)
        {
            vars.fbar[i] = local_vars.fbar[i];
            vars.q[i] = local_vars.q[i];
            FOR1(j)
            {
                vars.fspatial[i][j] = local_vars.fspatial[i][j];
                vars.v[i][j] = local_vars.v[i][j];
            }
        }

        // Calculate the derivatives
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        // Traceless condition

        // Traceless for derivatives

        // test
        Tensor<3, data_t> chris_local;
        FOR3(i, j, k)
        {
            chris_local[i][j][k] = 0.0;
            FOR1(l)
            {
                chris_local[i][j][k] +=
                    gamma_UU[i][l] * (metric_vars.d1_gamma[l][k][j] +
                                      metric_vars.d1_gamma[j][l][k] -
                                      metric_vars.d1_gamma[j][k][l]);
            }
            chris_local[i][j][k] /= 2.0;
        }

        Tensor<2, data_t> i_u;

        FOR2(i, j)
        {
            i_u[i][j] = d1.fbar[j][i] - vars.fhat * metric_vars.K_tensor[i][j];

            FOR1(k)
            {
                i_u[i][j] += -chris_phys.ULL[k][j][i] * vars.fbar[k];

                FOR1(l)
                {
                    i_u[i][j] += -gamma_UU[l][k] * vars.fspatial[j][k] *
                                 metric_vars.K_tensor[i][l];
                }
            }
        }

        Tensor<1, data_t> i_p;
        //  p <-> q
        FOR1(i)
        {
            i_p[i] = d1.fhat[i];

            FOR2(j, k)
            {
                i_p[i] += -2.0 * gamma_UU[j][k] * vars.fbar[j] *
                          metric_vars.K_tensor[i][k];
            }
        }

        Tensor<3, data_t> cd1_K_tensor;
        FOR3(i, j, k)
        {
            cd1_K_tensor[i][j][k] = metric_vars.d1_K_tensor[i][j][k];
            FOR1(l)
            {
                cd1_K_tensor[i][j][k] +=
                    -chris_phys.ULL[l][k][i] * metric_vars.K_tensor[l][j] -
                    chris_phys.ULL[l][k][j] * metric_vars.K_tensor[i][l];
            }
        }

        // NEEDS MASS TERM ADDING
        data_t primaryScalar = 0.0;

        FOR1(i)
        {
            primaryScalar += vars.w * metric_vars.shift[i] *
                             metric_vars.d1_lapse[i] / metric_vars.lapse;
            FOR1(j)
            {
                primaryScalar += gamma_UU[i][j] * (i_p[j] - vars.q[j]) *
                                 metric_vars.d1_lapse[i];
                FOR2(k, l)
                {
                    primaryScalar +=
                        -2.0 * metric_vars.lapse * gamma_UU[i][k] *
                            gamma_UU[j][l] *
                            (metric_vars.K * metric_vars.K_tensor[i][j] +
                             metric_vars.ricci_phys[i][j]) *
                            vars.fspatial[k][l] +
                        metric_vars.lapse * gamma_UU[i][k] * gamma_UU[j][l] *
                            i_u[i][j] * metric_vars.K_tensor[k][l];

                    FOR2(m, n)
                    {
                        primaryScalar +=
                            metric_vars.lapse * gamma_UU[i][l] *
                            gamma_UU[j][m] * gamma_UU[k][n] *
                            (metric_vars.K_tensor[i][m] *
                             metric_vars.K_tensor[j][k] * vars.fspatial[l][n]);
                    }
                }
            }
        }
        Tensor<1, data_t> primaryVector;
        FOR1(i)
        {
            primaryVector[i] =
                vars.w * metric_vars.d1_lapse[i] / metric_vars.lapse;

            FOR1(j)
            {

                FOR1(k)
                {
                    primaryVector[i] +=
                        -gamma_UU[j][k] * metric_vars.K_tensor[k][i] * i_p[j] -
                        2.0 * gamma_UU[j][k] * metric_vars.K_tensor[j][i] *
                            metric_vars.K * vars.fbar[k] -
                        2.0 * gamma_UU[j][k] * metric_vars.ricci_phys[i][j] *
                            vars.fbar[k];

                    FOR2(l, m)
                    {
                        primaryVector[i] +=
                            -2.0 * gamma_UU[j][l] * gamma_UU[k][m] *
                                (cd1_K_tensor[j][k][i] -
                                 cd1_K_tensor[i][j][k]) *
                                vars.fspatial[l][m] +
                            2.0 * gamma_UU[j][l] * gamma_UU[k][m] *
                                metric_vars.K_tensor[k][j] *
                                metric_vars.K_tensor[i][m] * vars.fbar[l];
                    }
                }
            }
        }

        Tensor<1, data_t> transverseVector;

        /*
        FOR1(i)
        {
            vars.q[i] = -metric_vars.K * vars.fbar[i];
            FOR2(j,k)
            {
              vars.q[i] += -gamma_UU[j][k] * metric_vars.K_tensor[i][k] *
        vars.fbar[j]; vars.q[i] += gamma_UU[j][k] * d1.fspatial[j][i][k];
              FOR1(l)
              {
                vars.q[i] += -gamma_UU[j][k] * (chris_local[l][k][i] *
        vars.fspatial[j][l] + chris_local[l][k][j] * vars.fspatial[l][i]);
              }
            }
        }
        */
        FOR2(i, j)
        {
            vars.w += gamma_UU[i][j] * d1.fbar[i][j];
            FOR1(k)
            {
                vars.w +=
                    -gamma_UU[i][j] * chris_phys.ULL[k][i][j] * vars.fbar[k];
            }
        }

        current_cell.store_vars(vars);
    }
};

#endif /* SETREST_HPP_ */
