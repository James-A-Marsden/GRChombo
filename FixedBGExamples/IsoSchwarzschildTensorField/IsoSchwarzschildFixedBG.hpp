/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ISOSCHWARZSCHILDFIXEDBG_HPP_
#define ISOSCHWARZSCHILDFIXEDBG_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

//! Class which computes the initial conditions for a Schwarzschild BH (a/M=0)
//! in Isotropic coordinates
class IsoSchwarzschildFixedBG
{
  public:
    //! Struct for the params of the  BH
    struct params_t
    {
        double mass = 1.0;                      //!<< The mass of the BH
        std::array<double, CH_SPACEDIM> center; //!< The center of the BH
    };

    template <class data_t> using Vars = ADMFixedBGVars::TensorVars<data_t>;

    const params_t m_params;
    const double m_dx;

    IsoSchwarzschildFixedBG(params_t a_params, double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
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
        current_cell.store_vars(metric_vars.d2_lapse[0][1], c_chi);
    }

    // Schwarzschild solution
    template <class data_t, template <typename> class vars_t>
    void compute_metric_background(vars_t<data_t> &vars,
                                   const Cell<data_t> &current_cell) const
    {
        using namespace TensorAlgebra;
        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

        // black hole params - mass M
        const double M = m_params.mass;

        // work out where we are on the grid including effect of spin
        // on x direction (length contraction)
        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const data_t r = coords.get_radius();
        const data_t r2 = r * r;

        Tensor<1, data_t> d1_r;
        Tensor<2, data_t> d2_r;
        // For neatness, the coordinates in a vector
        Tensor<1, data_t> X;
        X[0] = coords.x;
        X[1] = coords.y;
        X[2] = coords.z;

        // First, get the derivatives of the radius
        FOR1(i) { d1_r[i] = X[i] / r; }

        FOR2(i, j) { d2_r[i][j] = (delta(i, j) - X[i] * X[j] / r2) / r; }

        // The spatial metric
        FOR2(i, j)
        {
            vars.gamma[i][j] = delta(i, j) * pow(1.0 + 0.5 * M / r, 4.0);
        }

        // Calculate partial derivative of spatial metric
        FOR3(i, j, k)
        {
            vars.d1_gamma[i][j][k] = -0.25 * M * pow(M + 2.0 * r, 3.0) *
                                     pow(r, -5.0) * d1_r[k] * delta(i, j);
        }

        FOR2(i, j)
        {
            FOR2(k, l)
            {
                vars.d2_gamma[i][j][k][l] =
                    -delta(i, j) *
                    (M * (M + 2.0 * r) * (M + 2.0 * r) *
                     (-(5.0 * M + 4.0 * r) * d1_r[k] * d1_r[l] +
                      r * (M + 2.0 * r) * d2_r[k][l])) *
                    0.25 * pow(r, -6.0);
            }
        }

        // Now the inverse metric
        const auto gamma_UU = compute_inverse_sym(vars.gamma);
        const auto chris_phys = compute_christoffel(vars.d1_gamma, gamma_UU);
        FOR3(i, j, k)
        {
            vars.d1_gamma_UU[i][j][k] = delta(i, j) * 64.0 * M * d1_r[k] *
                                        pow(r, 3.0) * pow(M + 2.0 * r, -5.0);
        }

        // This is used but is probably not needed, as one can always take the
        // spatial metric out of the covariant derivative before expanding it
        FOR2(i, j)
        {
            FOR2(k, l)
            {
                vars.d2_gamma_UU[i][j][k][l] =
                    delta(i, j) * 64.0 * M * r2 *
                    (d1_r[k] * d1_r[l] * (3.0 * M - 4.0 * r) +
                     r * d2_r[k][l] * (2.0 * r + M)) *
                    pow(2.0 * r + M, -6.0);
            }
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

        vars.lapse = (1.0 - 0.5 * M / r) / (1.0 + 0.5 * M / r);
        // calculate derivs of lapse and shift
        FOR1(i)
        {
            vars.d1_lapse[i] =
                4.0 * M * d1_r[i] / (M + 2.0 * r) / (M + 2.0 * r);

            // vars.d1_ln_lapse[i] = -4.0 * M * d1_r[i] / (M * M - 4.0 * r2);
        }

        FOR2(i, j)
        {
            vars.d2_lapse[i][j] =
                4.0 * M *
                (d2_r[i][j] * (M + 2.0 * r) - 4.0 * d1_r[i] * d1_r[j]) *
                pow(M + 2.0 * r, -3.0);

            // vars.d2_ln_lapse[i][j] =
            //     -(4.0 * M *
            //       (8.0 * r * d1_r[i] * d1_r[j] + M * M * d2_r[i][j] -
            //        4.0 * r2 * d2_r[i][j])) /
            //     (M * M - 4.0 * r2) / (M * M - 4.0 * r2);
        }

        // shift is zero
        FOR1(i) { vars.shift[i] = 0.0; }
        FOR2(i, j) { vars.d1_shift[i][j] = 0.0; }
        FOR3(i, j, k) { vars.d2_shift[i][j][k] = 0.0; }

        // calculate the extrinsic curvature, using the fact that
        // 2 * lapse * K_ij = D_i \beta_j + D_j \beta_i - dgamma_ij dt
        // and dgamma_ij dt = 0 in chosen fixed gauge
        FOR2(i, j) { vars.K_tensor[i][j] = 0.0; }
        vars.K = compute_trace(gamma_UU, vars.K_tensor);
        FOR3(i, j, k) { vars.d1_K_tensor[i][j][k] = 0.0; }

        // Derivative of the trace, \partial_i K = \partial_i(gamma^jk K_jk)
        FOR1(i) { vars.d1_K[i] = 0.0; }

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
                        chris_phys.ULL[m][l][j] * chris_phys.ULL[i][m][k] -
                        chris_phys.ULL[m][k][j] * chris_phys.ULL[i][m][l];
                }
            }
        }

        // spatial ricci tensor
        FOR2(i, j)
        {
            vars.ricci_phys[i][j] = 0.0;
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

        // work out where we are on the grid
        const Coordinates<double> coords(current_cell, m_dx, m_params.center);
        const double x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const double r = sqrt(x * x + y * y + z * z);

        const double horizon_distance = r / (M / 2.0);

        return horizon_distance;
    }
};

#endif /* ISOSCHWARZSCHILDFIXEDBG_HPP_ */
