/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDBGDIAGNOSTICS_HPP_
#define FIXEDBGDIAGNOSTICS_HPP_

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
template <class matter_t, class background_t> class FixedBGDiagnostics
{
    // Use the variable definition in the matter class
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const matter_t m_matter;         //!< The matter object
    const double m_dx;               //!< The matter object
    const background_t m_background; //!< The matter object
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_tensor_field_mass;

  public:
    FixedBGDiagnostics(matter_t a_matter, background_t a_background,
                       double a_dx, std::array<double, CH_SPACEDIM> a_center,
                       double a_tensor_field_mass)
        : m_matter(a_matter), m_deriv(a_dx), m_dx(a_dx),
          m_background(a_background), m_center(a_center),
          m_tensor_field_mass(a_tensor_field_mass)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables, and calc 1st
        // derivs
        const auto vars = current_cell.template load_vars<MatterVars>();
        const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);
        const auto d2 = m_deriv.template diff2<MatterVars>(current_cell);
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // get the metric vars
        MetricVars<data_t> metric_vars;
        m_background.compute_metric_background(metric_vars, current_cell);

        using namespace TensorAlgebra;
        const auto gamma_UU =
            TensorAlgebra::compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys =
            compute_christoffel(metric_vars.d1_gamma, gamma_UU);
        const emtensor_t<data_t> emtensor = m_matter.compute_emtensor(
            vars, metric_vars, d1, gamma_UU, chris_phys.ULL);
        const data_t det_gamma =
            TensorAlgebra::compute_determinant_sym(metric_vars.gamma);

        // Coordinates
        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const data_t x2 = x * x;
        const double y2 = y * y;
        const double z2 = z * z;
        const data_t r = coords.get_radius();
        const data_t r2 = r * r;
        const data_t rho = simd_max(sqrt(x2 + y2), 1e-6);
        const data_t rho2 = rho * rho;

        const data_t costheta = z / r;
        const data_t sintheta = rho / r;

        const data_t sinphi = y / rho;
        const data_t cosphi = x / rho;

        const data_t sin2phi = 2.0 * sinphi * cosphi;
        const data_t cos2phi = cosphi * cosphi - sinphi * sinphi;

        const double horizon = 1.0 / 2.0;

        const double xx = coords.x * coords.x;
        const double yy = coords.y * coords.y;
        const double zz = coords.z * coords.z;
        const double rr = sqrt(xx + yy + zz);

        data_t primaryScalar = 0.0;

        Tensor<1, data_t> primaryVector;
        FOR1(i) { primaryVector[i] = 0.0; }

        data_t trace_field = 0.0;
        data_t trace_momentum = 0.0;

        const double mass = 0.0;
        data_t rho_eff = 0.0;

        if (rr > horizon)
        {
            // replacement of fhat
            data_t fspatial_trace = 0.0;
            FOR2(i, j)
            {
                fspatial_trace +=
                    -metric_vars.lapse * gamma_UU[i][j] * vars.fspatial[i][j];
            }
            // Derivative of the trace of fspatial
            Tensor<1, data_t> d1_fspatial_trace;

            FOR1(i)
            {
                d1_fspatial_trace[i] = 0.0;

                FOR2(j, k)
                {
                    d1_fspatial_trace[i] +=
                        -metric_vars.d1_lapse[i] * gamma_UU[j][k] *
                            vars.fspatial[j][k] -
                        metric_vars.lapse *
                            (metric_vars.d1_gamma_UU[j][k][i] *
                                 vars.fspatial[j][k] +
                             gamma_UU[j][k] * d1.fspatial[j][k][i]);
                }
            }
            Tensor<2, data_t> d2_fspatial_trace;
            FOR2(i, j)
            {
                d2_fspatial_trace[i][j] = 0.0;

                FOR2(k, l)
                {
                    d2_fspatial_trace[i][j] +=
                        -metric_vars.d2_lapse[i][j] * gamma_UU[k][l] *
                            vars.fspatial[k][l] -
                        metric_vars.d1_lapse[i] *
                            metric_vars.d1_gamma_UU[k][l][j] *
                            vars.fspatial[k][l] -
                        metric_vars.d1_lapse[i] * gamma_UU[k][l] *
                            d1.fspatial[k][l][j] -
                        metric_vars.d1_lapse[j] *
                            metric_vars.d1_gamma_UU[k][l][i] *
                            vars.fspatial[k][l] -
                        metric_vars.lapse *
                            metric_vars.d2_gamma_UU[k][l][i][j] *
                            vars.fspatial[k][l] -
                        metric_vars.lapse * metric_vars.d1_gamma_UU[k][l][i] *
                            d1.fspatial[k][l][j] -
                        metric_vars.d1_lapse[j] * gamma_UU[k][l] *
                            d1.fspatial[k][l][i] -
                        metric_vars.lapse * metric_vars.d1_gamma_UU[k][l][j] *
                            d1.fspatial[k][l][i] -
                        metric_vars.lapse * gamma_UU[k][l] *
                            d2.fspatial[k][l][i][j];
                }
            }
            // TRACE DIAGNOSTICS

            trace_field = vars.fhat / metric_vars.lapse;

            FOR2(i, j) { trace_field += gamma_UU[i][j] * vars.fspatial[i][j]; }

            FOR2(i, j) { trace_momentum += gamma_UU[i][j] * vars.v[i][j]; }

            FOR2(i, j)
            {
                trace_momentum += -gamma_UU[i][j] * d1.fbar[i][j];
                FOR1(k)
                {
                    trace_momentum +=
                        gamma_UU[i][j] * chris_phys.ULL[k][i][j] * vars.fbar[k];
                }
            }

            rho_eff = 0.0;
            FOR3(i, j, k)
            {
                FOR1(l)
                {
                    rho_eff += 0.25 * m_tensor_field_mass *
                               m_tensor_field_mass * gamma_UU[i][k] *
                               gamma_UU[j][l] * vars.fspatial[i][j] *
                               vars.fspatial[k][l];
                }
            }

            primaryScalar = metric_vars.lapse * m_tensor_field_mass *
                            m_tensor_field_mass * fspatial_trace;

            FOR2(i, j)
            {
                primaryScalar +=
                    metric_vars.lapse * gamma_UU[i][j] *
                        (fspatial_trace * metric_vars.ricci_phys[i][j] +
                         2.0 * d1_fspatial_trace[i] *
                             metric_vars.d1_ln_lapse[j] -
                         2.0 * fspatial_trace * metric_vars.d1_ln_lapse[i] *
                             metric_vars.d1_ln_lapse[j] -
                         d2_fspatial_trace[i][j]) +
                    gamma_UU[i][j] * fspatial_trace *
                        metric_vars.d2_lapse[i][j];

                FOR1(k)
                {
                    primaryScalar +=
                        metric_vars.lapse * gamma_UU[i][j] *
                        (chris_phys.ULL[k][i][j] * d1_fspatial_trace[k] -
                         fspatial_trace * chris_phys.ULL[k][i][j] *
                             metric_vars.d1_ln_lapse[k]);

                    FOR1(l)
                    {
                        primaryScalar += gamma_UU[i][k] * gamma_UU[j][l] *
                                         metric_vars.lapse * metric_vars.lapse *
                                         metric_vars.ricci_phys[k][l] *
                                         vars.fspatial[i][j];
                        primaryScalar +=
                            gamma_UU[i][k] * gamma_UU[j][l] *
                            metric_vars.lapse *
                            (-metric_vars.lapse * d2.fspatial[i][j][k][l]);

                        FOR1(m)
                        {
                            primaryScalar +=
                                -gamma_UU[i][k] * gamma_UU[j][l] *
                                metric_vars.lapse * metric_vars.lapse *
                                (-chris_phys.ULL[m][l][k] *
                                     d1.fspatial[i][j][m] -
                                 chris_phys.ULL[m][l][i] *
                                     d1.fspatial[m][j][k] -
                                 chris_phys.ULL[m][l][j] *
                                     d1.fspatial[i][m][k] -
                                 metric_vars.d1_chris_phys[m][k][i][l] *
                                     vars.fspatial[m][j] -
                                 chris_phys.ULL[m][k][i] *
                                     d1.fspatial[m][j][l] -
                                 metric_vars.d1_chris_phys[m][j][k][l] *
                                     vars.fspatial[i][m] -
                                 chris_phys.ULL[m][j][k] *
                                     d1.fspatial[i][m][l]);

                            FOR1(n)
                            {
                                primaryScalar +=
                                    -gamma_UU[i][k] * gamma_UU[j][l] *
                                    metric_vars.lapse * metric_vars.lapse *
                                    (chris_phys.ULL[n][l][k] *
                                         chris_phys.ULL[m][n][i] *
                                         vars.fspatial[m][j] +
                                     chris_phys.ULL[n][l][i] *
                                         chris_phys.ULL[m][k][n] *
                                         vars.fspatial[m][j] +
                                     chris_phys.ULL[n][l][j] *
                                         chris_phys.ULL[m][k][i] *
                                         vars.fspatial[m][n] +
                                     chris_phys.ULL[n][l][k] *
                                         chris_phys.ULL[m][n][j] *
                                         vars.fspatial[i][m] +
                                     chris_phys.ULL[n][l][j] *
                                         chris_phys.ULL[m][k][n] *
                                         vars.fspatial[i][m] +
                                     chris_phys.ULL[n][l][i] *
                                         chris_phys.ULL[m][k][j] *
                                         vars.fspatial[n][m]);
                            }
                        }
                    }
                }
            }

            FOR1(i)
            {
                primaryVector[i] = metric_vars.lapse * m_tensor_field_mass *
                                   m_tensor_field_mass * vars.fbar[i];

                FOR2(j, k)
                {
                    primaryVector[i] +=
                        metric_vars.lapse * gamma_UU[j][k] *
                        (d1.v[j][i][k] +
                         vars.fbar[i] * metric_vars.ricci_phys[j][k] -
                         vars.fbar[j] * metric_vars.ricci_phys[i][k] -
                         d2.fbar[i][j][k]);

                    FOR1(l)
                    {
                        primaryVector[i] +=
                            metric_vars.lapse * gamma_UU[j][k] *
                            (-chris_phys.ULL[l][k][i] * vars.v[j][l] -
                             chris_phys.ULL[l][k][j] * vars.v[l][i] +
                             chris_phys.ULL[l][k][j] * d1.fbar[i][l] +
                             chris_phys.ULL[l][k][i] * d1.fbar[l][j] +
                             metric_vars.d1_chris_phys[l][j][i][k] *
                                 vars.fbar[l] +
                             d1.fbar[l][k] * chris_phys.ULL[l][j][i]);

                        FOR1(m)
                        {
                            primaryVector[i] +=
                                metric_vars.lapse * gamma_UU[j][k] *
                                (-chris_phys.ULL[m][k][j] *
                                     chris_phys.ULL[l][m][i] * vars.fbar[l] -
                                 chris_phys.ULL[m][k][i] *
                                     chris_phys.ULL[l][j][m] * vars.fbar[l]);
                        }
                    }
                }
            }
        }
        current_cell.store_vars(trace_field, c_trace_field);
        current_cell.store_vars(trace_momentum, c_trace_momentum);

        current_cell.store_vars(primaryScalar, c_primaryConstraintScalar);

        current_cell.store_vars(primaryVector[0], c_primaryConstraintVector1);
        current_cell.store_vars(primaryVector[1], c_primaryConstraintVector2);
        current_cell.store_vars(primaryVector[2], c_primaryConstraintVector3);

        current_cell.store_vars(rho_eff, c_rho_eff);
    }
};

#endif /* FIXEDBGDIAGNOSTICS_HPP_ */
