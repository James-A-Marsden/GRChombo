/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TAUBFIXEDBG_HPP_
#define TAUBFIXEDBG_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

//! Class which computes the initial conditions for a Kerr Schild BH
//! https://arxiv.org/pdf/gr-qc/9805023.pdf
class TaubFixedBG
{
  public:
    //! Struct for the params of the  BH
    struct params_t
    {
        double mass = 1.0;                      //!<< The mass of the BH
        std::array<double, CH_SPACEDIM> center; //!< The center of the BH
        double spin = 0.0;                      //!< The spin param a = J / M
    };

    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;

    const params_t m_params;
    const double m_dx;

    TaubFixedBG(params_t a_params, double a_dx) : m_params(a_params), m_dx(a_dx)
    {
        // check this spin param is sensible
        if ((m_params.spin > m_params.mass) || (m_params.spin < -m_params.mass))
        {
            MayDay::Error(
                "The dimensionless spin parameter must be in the range "
                "-1.0 < spin < 1.0");
        }
    }

    /// This just calculates chi which helps with regridding, debug etc
    /// it is only done once on setup as the BG is fixed
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // get position and set vars
        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        Vars<data_t> metric_vars;
        compute_metric_background(metric_vars, current_cell);

        // calculate and save chi
        data_t chi = TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
        chi = pow(chi, -1.0 / 3.0);
        current_cell.store_vars(chi, c_chi);
    }

    // Kerr Schild solution
    template <class data_t, template <typename> class vars_t>
    void compute_metric_background(vars_t<data_t> &vars,
                                   const Cell<data_t> &current_cell) const
    {
        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

        // black hole params - mass M and spin a
        const double M = m_params.mass;
        const double a = m_params.spin;
        const double a2 = a * a;

        // work out where we are on the grid including effect of spin
        // on x direction (length contraction)
        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z + 200.0;

        vars.lapse = pow(z, -1.0 / 3.0);
        FOR2(i, j) { vars.gamma[i][j] = 0.0; }
        vars.gamma[0][0] = pow(z, 4.0 / 3.0);
        vars.gamma[1][1] = pow(z, 4.0 / 3.0);
        vars.gamma[2][2] = 1.0;
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse_sym(vars.gamma);
        FOR1(i) { vars.shift[i] = 0; }

        // Calculate partial derivative of spatial metric
        FOR3(i, j, k) { vars.d1_gamma[i][j][k] = 0.0; }
        FOR1(k)
        {
            vars.d1_gamma[0][0][k] =
                (4.0 / 3.0) * delta(k, 2) * pow(z, 1.0 / 3.0);
            vars.d1_gamma[1][1][k] =
                (4.0 / 3.0) * delta(k, 2) * pow(z, 1.0 / 3.0);
        }

        vars.gamma_UU = compute_inverse_sym(vars.gamma);

        FOR3(i, j, k) { vars.d1_gamma_UU[i][j][k] = 0.0; }
        FOR1(k)
        {
            vars.d1_gamma_UU[0][0][k] =
                (-4.0 / 3.0) * delta(k, 2) * pow(z, -7.0 / 3.0);
            vars.d1_gamma_UU[1][1][k] =
                (-4.0 / 3.0) * delta(k, 2) * pow(z, -7.0 / 3.0);
        }
        // Second derivative of the spatial metric (can simplify further)
        FOR1(i)
        {
            FOR3(j, k, m) { vars.d2_gamma[i][j][k][m] = 0.0; }
        }
        FOR2(k, m)
        {
            vars.d2_gamma[0][0][k][m] =
                (4.0 / 9.0) * delta(k, 2) * delta(m, 2) * pow(z, -2.0 / 3.0);
            vars.d2_gamma[1][1][k][m] =
                (4.0 / 9.0) * delta(k, 2) * delta(m, 2) * pow(z, -2.0 / 3.0);
        }

        // Calculate derivative of the Christoffel symbol (phys)
        FOR2(i, j)
        {
            FOR2(k, m)
            {
                vars.d1_chris_phys[i][j][k][m] = 0.0;

                FOR1(n)
                {
                    vars.d1_chris_phys[i][j][k][m] +=
                        0.5 * vars.d1_gamma_UU[i][n][m] *
                            (vars.d1_gamma[k][n][j] + vars.d1_gamma[n][j][k] -
                             vars.d1_gamma[j][k][n])

                        + 0.5 * gamma_UU[i][n] *
                              (vars.d2_gamma[k][n][j][m] +
                               vars.d2_gamma[n][j][k][m] -
                               vars.d2_gamma[j][k][n][m]);
                }
            }
        }

        data_t alpha2 = vars.lapse * vars.lapse;
        // calculate derivs of lapse and shift
        FOR1(i)
        {
            vars.d1_lapse[i] = (-1.0 / 3.0) * delta(i, 2) * pow(z, -4.0 / 3.0);
        }

        FOR2(i, j)
        {
            vars.d2_lapse[i][j] =
                (4.0 / 9.0) * delta(i, 2) * delta(j, 2) * pow(z, -7.0 / 3.0);
        }

        // use the fact that shift^i = lapse^2 * shift_i
        FOR2(i, j) { vars.d1_shift[i][j] = 0.0; }

        FOR3(i, j, k) { vars.d2_shift[i][j][k] = 0.0; }
        // calculate the extrinsic curvature, using the fact that
        // 2 * lapse * K_ij = D_i \beta_j + D_j \beta_i - dgamma_ij dt
        // and dgamma_ij dt = 0 in chosen fixed gauge
        // const auto chris_phys = compute_christoffel(vars.d1_gamma, gamma_UU);

        Tensor<3, data_t> chris_local;
        FOR3(i, j, k)
        {
            chris_local[i][j][k] = 0.0;
            FOR1(l)
            {
                chris_local[i][j][k] +=
                    0.5 * gamma_UU[i][l] *
                    (vars.d1_gamma[l][k][j] + vars.d1_gamma[j][l][k] -
                     vars.d1_gamma[j][k][l]);
            }
        }

        FOR2(i, j)
        {
            /*
            vars.K_tensor[i][j] = 0.0;
            FOR1(k)
            {
                vars.K_tensor[i][j] +=
                    vars.gamma[k][j] * vars.d1_shift[k][i] +
                    vars.gamma[k][i] * vars.d1_shift[k][j] +
                    (vars.d1_gamma[k][i][j] + vars.d1_gamma[k][j][i]) *
                        vars.shift[k];
                FOR1(m)
                {
                    vars.K_tensor[i][j] += -2.0 * chris_local[k][i][j] *
                                           vars.gamma[k][m] * vars.shift[m];
                }
            }
            vars.K_tensor[i][j] *= 0.5 / vars.lapse;
            */
            vars.K_tensor[i][j] = 0.0;
        }
        vars.K = compute_trace(gamma_UU, vars.K_tensor);
        FOR3(i, j, k)
        {
            vars.d1_K_tensor[i][j][k] = 0.0;
            /*
            FOR1(m)
            {
                vars.d1_K_tensor[i][j][k] += vars.d1_gamma[m][j][k] *
            vars.d1_shift[m][i] + vars.gamma[m][j] * vars.d2_shift[m][i][k]
                                            + vars.d1_gamma[m][i][k] *
            vars.d1_shift[m][j] + vars.gamma[m][i] * vars.d2_shift[m][j][k]
                + (vars.d2_gamma[m][i][j][k] + vars.d2_gamma[m][j][i][k]) *
            vars.shift[m]
                + (vars.d1_gamma[m][i][j] + vars.d1_gamma[m][j][i]) *
            vars.d1_shift[m][k]; FOR1(n)
                {
                    vars.d1_K_tensor[i][j][k] += -2.0 *
            (vars.d1_chris_phys[m][i][j][k] * vars.gamma[m][n] * vars.shift[n]
                    + chris_local[m][i][j] * vars.d1_gamma[m][n][k] *
            vars.shift[n]
                    + chris_local[m][i][j] * vars.gamma[m][n] *
            vars.d1_shift[n][k]);
                }

            }
            vars.d1_K_tensor[i][j][k] *= 0.5 / vars.lapse;
            vars.d1_K_tensor[i][j][k] += - vars.d1_lapse[k] / vars.lapse *
            vars.K_tensor[i][j];
            */
            // vars.d1_K_tensor[i][j][k] = 0.0;
        }
        // Derivative of the trace, \partial_i K = \partial_i(gamma^jk K_jk)
        FOR1(i)
        {
            vars.d1_K[i] = 0;
            /*
            FOR2(j,k)
            {
                vars.d1_K[i] += vars.d1_gamma_UU[j][k][i] * vars.K_tensor[j][k]
            + gamma_UU[j][k] * vars.d1_K_tensor[j][k][i];
            }
            */
        }

        // spatial riemann curvature tensor

        FOR1(i)
        {
            FOR3(j, k, l)
            {
                vars.riemann_phys_ULLL[i][j][k][l] =
                    vars.d1_chris_phys[i][l][j][k] -
                    vars.d1_chris_phys[i][k][j][l];

                FOR1(m)
                {
                    vars.riemann_phys_ULLL[i][j][k][l] +=
                        chris_local[m][l][j] * chris_local[i][m][k] -
                        chris_local[m][k][j] * chris_local[i][m][l];
                }
            }
        }

        // spatial ricci tensor
        FOR2(i, j)
        {
            vars.ricci_phys[i][j] = 0;
            FOR1(k)
            {
                vars.ricci_phys[i][j] += vars.riemann_phys_ULLL[k][i][k][j];
            }
        }
    }

  public:
    // used to decide when to excise - ie when within the horizon of the BH
    // note that this is not templated over data_t
    double excise(const Cell<double> &current_cell) const
    {
        // black hole params - mass M and spin a
        const double M = m_params.mass;
        const double a = m_params.spin;
        const double a2 = a * a;

        // work out where we are on the grid
        const Coordinates<double> coords(current_cell, m_dx, m_params.center);
        const double x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const double r_plus = M + sqrt(M * M - a2);
        const double r_minus = M - sqrt(M * M - a2);

        const double outer_horizon =
            (x * x + y * y) / (2.0 * M * r_plus) + z * z / r_plus / r_plus;

        const double inner_horizon =
            (x * x + y * y) / (2.0 * M * r_minus) + z * z / r_minus / r_minus;

        // value less than 1 indicates we are within the horizon
        return sqrt(outer_horizon);
    }
};

#endif /* TAUBFIXEDBG_HPP_ */
