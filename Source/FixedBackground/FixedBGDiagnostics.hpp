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

        data_t lapse2 = metric_vars.lapse * metric_vars.lapse;
        data_t lapse3 = metric_vars.lapse * metric_vars.lapse * metric_vars.lapse;

        data_t primaryScalar = 0.0;

        data_t field_amplitude = 0.0;

        Tensor<2,data_t> D2_fhat;
        FOR2(i,j) {D2_fhat[i][j] = 0.0;}

        Tensor<2,data_t> D1_fbar;
        FOR2(i,j) {D1_fbar[i][j] = 0.0;}

        Tensor<3,data_t> D2_fbar;
        FOR3(i,j,k) {D2_fbar[i][j][k] = 0.0;}

        Tensor<3,data_t> D1_fspatial;
        FOR3(i,j,k) {D1_fspatial[i][j][k] = 0.0;}

        Tensor<4,data_t> D2_fspatial;
        FOR2(i,j)
        {
            FOR2(k,l)
            {
                D2_fspatial[i][j][k][l] = 0.0;
            }
        }

        Tensor<1,data_t> shift_L;
        FOR1(i) {shift_L[i] = 0.0;}

        Tensor<3,data_t> D1_K_tensor;
        FOR3(i,j,k) {D1_K_tensor[i][j][k] = 0.0;} 

        Tensor<2, data_t> D2_lapse;
        FOR2(i,j) {D2_lapse[i][j] = 0.0; }

    

        Tensor<1, data_t> primaryVector;
        FOR1(i) { primaryVector[i] = 0.0; }

        data_t fspatial_trace = 0.0;
        Tensor<1, data_t> d1_fspatial_trace;
        FOR1(i) { d1_fspatial_trace[i] = 0.0; }
        Tensor<2, data_t> d2_fspatial_trace;
        FOR2(i, j) { d2_fspatial_trace[i][j] = 0.0; }
        data_t trace_field = 0.0;
        data_t trace_momentum = 0.0;

        const double mass = 0.0;
        data_t rho_eff = 0.0;

        if (rr > horizon)
        {
            // replacement of fhat

            FOR2(i, j)
            {
                fspatial_trace +=
                    -metric_vars.lapse * gamma_UU[i][j] * vars.fspatial[i][j];
            }
            // Derivative of the trace of fspatial

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

            trace_field = vars.fhat;

            FOR2(i, j)
            {
                trace_field +=
                    gamma_UU[i][j] * vars.fspatial[i][j] * metric_vars.lapse;
            }

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

            rho_eff = - 0.75 * m_tensor_field_mass * m_tensor_field_mass * vars.fhat * vars.fhat / metric_vars.lapse / metric_vars.lapse;
            FOR2(i, j)
            {
                rho_eff += 0.5 *  m_tensor_field_mass *  m_tensor_field_mass * gamma_UU[i][j] * vars.fbar[i] * vars.fbar[j] / metric_vars.lapse / metric_vars.lapse;
                FOR1(k)
                {
                    FOR1(l)
                    {
                        rho_eff += 0.25 * m_tensor_field_mass *
                                m_tensor_field_mass * gamma_UU[i][k] *
                                gamma_UU[j][l] * vars.fspatial[i][j] *
                                vars.fspatial[k][l];
                    }
                }
            }


            primaryScalar = metric_vars.lapse * m_tensor_field_mass *
                            m_tensor_field_mass * fspatial_trace;


            field_amplitude = vars.fhat * vars.fhat / metric_vars.lapse / metric_vars.lapse;

            FOR2(i,j)
            {
                field_amplitude += -2.0 * gamma_UU[i][j] * vars.fbar[i] * vars.fbar[j] / metric_vars.lapse / metric_vars.lapse;

                FOR2(k,l)
                {
                    field_amplitude += gamma_UU[i][k] * gamma_UU[j][l] * vars.fspatial[i][j] * vars.fspatial[k][l];
                }
            }

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
            //Modifications to the primary scalar due to extrinsic curvature and shift
            FOR2(i,j)
            {
                FOR2(k,l)
                {
                    primaryScalar+= metric_vars.lapse * vars.fhat * gamma_UU[i][k] * gamma_UU[j][l] *( -2.0 * metric_vars.K_tensor[i][j] * metric_vars.K_tensor[k][l]
                                                                                                        +metric_vars.K_tensor[i][k] * metric_vars.K_tensor[j][l]);

                    primaryScalar += lapse2 * gamma_UU[i][k] * gamma_UU[j][l] * metric_vars.K_tensor[i][j] * vars.v[k][l];     

                    
                    FOR1(m)
                    {
                        FOR1(n)
                        {
                        primaryScalar += lapse2 * gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * metric_vars.K_tensor[i][l] * metric_vars.K_tensor[m][n] * vars.fspatial[j][k];

                        }

                    primaryScalar += gamma_UU[j][l] * gamma_UU[k][m] * (2.0 * metric_vars.K_tensor[l][m] * metric_vars.K_tensor[i][j] * vars.fbar[k] * metric_vars.shift[i]
                                                                        +2.0 * metric_vars.K_tensor[k][m] * metric_vars.K_tensor[i][j] * vars.fbar[l] * metric_vars.shift[i]
                                                                        +2.0 * metric_vars.K_tensor[j][m] * metric_vars.K_tensor[i][k] * vars.fbar[l] * metric_vars.shift[i]);

                  

                                                
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
                        (-metric_vars.lapse * d1.v[j][i][k] +
                         vars.fbar[i] * metric_vars.ricci_phys[j][k] -
                         vars.fbar[j] * metric_vars.ricci_phys[i][k] -
                         d2.fbar[i][j][k]);

                    primaryVector[i] +=
                        2.0 * gamma_UU[j][k] *
                        (metric_vars.d1_lapse[k] * d1.fbar[i][j] -
                         metric_vars.lapse * vars.fbar[i] *
                             metric_vars.d1_ln_lapse[j] *
                             metric_vars.d1_ln_lapse[k]);

                    primaryVector[i] += gamma_UU[j][k] * vars.fbar[i] *
                                        metric_vars.d2_lapse[j][k];

                    FOR1(l)
                    {
                        primaryVector[i] +=
                            metric_vars.lapse * gamma_UU[j][k] *
                            (metric_vars.lapse * chris_phys.ULL[l][k][i] *
                                 vars.v[j][l] +
                             metric_vars.lapse * chris_phys.ULL[l][k][j] *
                                 vars.v[l][i] +
                             chris_phys.ULL[l][k][j] * d1.fbar[i][l] +
                             chris_phys.ULL[l][k][i] * d1.fbar[l][j] +
                             metric_vars.d1_chris_phys[l][j][i][k] *
                                 vars.fbar[l] +
                             d1.fbar[l][k] * chris_phys.ULL[l][j][i]);

                        primaryVector[i] +=
                            -2.0 * gamma_UU[j][k] * metric_vars.d1_lapse[k] *
                            chris_phys.ULL[l][j][i] * vars.fbar[l];

                        primaryVector[i] += -gamma_UU[j][k] * vars.fbar[i] *
                                            chris_phys.ULL[l][j][k] *
                                            metric_vars.d1_lapse[l];

                        //Modification to Bij
                        /*
                        primaryVector[i] += metric_vars.lapse * gamma_UU[j][k] *
                                            (chris_phys.ULL[l][k][i] *
                                            (-d1.fbar[j][l] - d1.fbar[l][j]
                                            +vars.fbar[j] * metric_vars.d1_ln_lapse[l]
                                            +vars.fbar[l] * metric_vars.d1_ln_lapse[j]));

                        primaryVector[i] += metric_vars.lapse * gamma_UU[j][k] *
                                            (chris_phys.ULL[l][k][j] *
                                            (-d1.fbar[i][l] - d1.fbar[l][i]
                                            +vars.fbar[i] * metric_vars.d1_ln_lapse[l]
                                            +vars.fbar[l] * metric_vars.d1_ln_lapse[i]));
                        ///                    
                        */                        
                        FOR1(m)
                        {
                            primaryVector[i] +=
                                metric_vars.lapse * gamma_UU[j][k] *
                                (-chris_phys.ULL[m][k][j] *
                                     chris_phys.ULL[l][m][i] * vars.fbar[l] -
                                 chris_phys.ULL[m][k][i] *
                                     chris_phys.ULL[l][j][m] * vars.fbar[l]);
                            /*
                            //Modification to Bij
                            primaryVector[i] += metric_vars.lapse * gamma_UU[j][k] *
                                            (2.0 * chris_phys.ULL[l][k][i] * chris_phys.ULL[m][j][l] * vars.fbar[m]
                                            +2.0 * chris_phys.ULL[l][k][j] * chris_phys.ULL[m][i][l] * vars.fbar[m]);
                            //                                         
                            */
                        }
                    }
                }
            }
        }

        //Modifications due to extrinsic curvature

                                          

        
        FOR2(i,j)
        {
            FOR2(k,l)
            {
                primaryScalar += -2.0 * metric_vars.d1_lapse[l] * metric_vars.K_tensor[i][j] * vars.fbar[k] * gamma_UU[i][k] * gamma_UU[j][l]
                                     + metric_vars.d1_lapse[i] * metric_vars.K_tensor[l][j] * vars.fbar[k] * gamma_UU[i][k] * gamma_UU[j][l]
                                     + 2.0 * D1_fbar[k][l] * metric_vars.K_tensor[i][j] * gamma_UU[i][k] * gamma_UU[j][l] * metric_vars.lapse
                                     - D1_fbar[k][i] * metric_vars.K_tensor[l][j] * gamma_UU[i][k] * gamma_UU[j][l] * metric_vars.lapse
                                      - 2.0 * D1_K_tensor[i][l][j] * vars.fbar[k] * gamma_UU[i][k] * gamma_UU[j][l] * metric_vars.lapse
                                       - 2.0 * metric_vars.K_tensor[i][j] * metric_vars.K_tensor[k][l] * gamma_UU[i][k] * gamma_UU[j][l] * vars.fhat * metric_vars.lapse
                                       + metric_vars.K_tensor[k][i] * metric_vars.K_tensor[l][j] * gamma_UU[i][k] * gamma_UU[j][l] * vars.fhat * metric_vars.lapse
                                       + vars.v[k][l] * metric_vars.K_tensor[i][j] * gamma_UU[i][k] * gamma_UU[j][l] * lapse2
                                       + D2_fbar[l][i][j] * gamma_UU[i][k] * gamma_UU[j][l] * shift_L[k]
                                       - d1.fhat[l] * metric_vars.K_tensor[i][j] * gamma_UU[i][k] * gamma_UU[j][l] * shift_L[k]
                                       - d1.fhat[i] * metric_vars.K_tensor[l][j] * gamma_UU[i][k] * gamma_UU[j][l] * shift_L[k]
                                        - D1_K_tensor[l][j][i] * gamma_UU[i][k] * gamma_UU[j][l] * vars.fhat * shift_L[k]
                                         - (D1_fbar[l][j] * metric_vars.d1_lapse[i] * gamma_UU[i][k] * gamma_UU[j][l] * shift_L[k])/metric_vars.lapse
                                          -  (D1_fbar[l][i] * metric_vars.d1_lapse[j] * gamma_UU[i][k] * gamma_UU[j][l] * shift_L[k])/metric_vars.lapse
                                          + (metric_vars.d1_lapse[l] * metric_vars.K_tensor[i][j] * gamma_UU[i][k] * gamma_UU[j][l] * vars.fhat * shift_L[k])/metric_vars.lapse
                                          + (metric_vars.d1_lapse[i] * metric_vars.K_tensor[l][j] * gamma_UU[i][k] * gamma_UU[j][l] * vars.fhat * shift_L[k])/metric_vars.lapse
                                          + (2.0 * metric_vars.d1_lapse[i] * metric_vars.d1_lapse[j] * vars.fbar[k] * gamma_UU[i][k] * gamma_UU[j][l] * shift_L[l])/lapse2
                                          -  (D2_lapse[j][i] * vars.fbar[k] * gamma_UU[i][k] * gamma_UU[j][l] * shift_L[l])/metric_vars.lapse
                                           + (2.0 * metric_vars.d1_lapse[i] * metric_vars.d1_lapse[j] * gamma_UU[i][k] * gamma_UU[j][l] * vars.fhat * shift_L[k] * shift_L[l])/lapse3
                                           -  (2.0 * d1.fhat[i] * metric_vars.d1_lapse[j] * gamma_UU[i][k] * gamma_UU[j][l] * shift_L[k] * shift_L[l])/lapse2
                                          -  (D2_lapse[j][i] * gamma_UU[i][k] * gamma_UU[j][l] * vars.fhat * shift_L[k] * shift_L[l])/lapse2
                                          + (D2_fhat[j][i] * gamma_UU[i][k] * gamma_UU[j][l] * shift_L[k] * shift_L[l])/metric_vars.lapse;
                                       

                FOR1(m)
                {
                    FOR1(n)
                    {
                        primaryScalar += 2.0 * D1_K_tensor[i][j][k] * vars.fbar[l] * metric_vars.gamma[m][n] * gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * metric_vars.lapse
                                        + metric_vars.K_tensor[l][i] * metric_vars.K_tensor[m][n] * gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * vars.fspatial[j][k] * lapse2
                                        + 2.0 * D1_fspatial[j][n][k] * metric_vars.K_tensor[i][m] * gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * metric_vars.lapse * shift_L[l]
                                        + D1_fspatial[j][k][i] * metric_vars.K_tensor[m][n] * gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * metric_vars.lapse * shift_L[l]
                                        + D1_K_tensor[j][k][i] * gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * vars.fspatial[m][n] * metric_vars.lapse * shift_L[l]
                                        + 4.0 * metric_vars.K_tensor[i][n] * metric_vars.K_tensor[j][k] * vars.fbar[l] * gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * shift_L[m]
                                        + 2.0 * metric_vars.K_tensor[i][j] * metric_vars.K_tensor[n][k] * vars.fbar[l] * gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * shift_L[m]
                                        - (D1_fbar[n][k] * metric_vars.K_tensor[i][j] * gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * shift_L[l] * shift_L[m])/ metric_vars.lapse
                                        - (4.0 * D1_fbar[n][i] * metric_vars.K_tensor[j][k] * gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * shift_L[l] * shift_L[m])/metric_vars.lapse
                                        + (2.0 * metric_vars.K_tensor[i][n] * metric_vars.K_tensor[j][k] * gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * vars.fhat * shift_L[l] * shift_L[m])/metric_vars.lapse
                                        + (metric_vars.K_tensor[i][j] * metric_vars.K_tensor[n][k] * gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * vars.fhat * shift_L[l] * shift_L[m])/metric_vars.lapse
                                        + (4.0 * metric_vars.d1_lapse[j] * metric_vars.K_tensor[i][k] * vars.fbar[l] * gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * shift_L[m] * shift_L[n])/lapse2 
                                        + (metric_vars.d1_lapse[i] * metric_vars.K_tensor[j][k] * vars.fbar[l] * gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * shift_L[m] * shift_L[n])/lapse2
                                         - (2.0 * D1_K_tensor[i][j][k] * vars.fbar[l] * gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * shift_L[m] * shift_L[n])/metric_vars.lapse;
                        FOR1(o)
                        {
                            FOR1(p)
                            {
                                primaryScalar += - 2.0 * metric_vars.K_tensor[i][o] * metric_vars.K_tensor[j][p] * gamma_UU[i][m] * gamma_UU[j][n] * gamma_UU[k][o] * gamma_UU[l][p] * vars.fspatial[k][l] * shift_L[m] * shift_L[n]
                                                    - metric_vars.K_tensor[i][j] * metric_vars.K_tensor[o][p] * gamma_UU[i][m] * gamma_UU[j][n] * gamma_UU[k][o] * gamma_UU[l][p] * vars.fspatial[k][l] * shift_L[m] * shift_L[n]; 
                            }
                        }
                    }
                }
            }
        }
                                        
        FOR2(i,j)
        {
            FOR1(k)
            {
                primaryVector[i] += - 2.0 * metric_vars.d1_lapse[k] * metric_vars.K_tensor[i][j] * gamma_UU[j][k] * vars.fhat
                                    + 2.0 * d1.fhat[k] * metric_vars.K_tensor[i][j] * gamma_UU[j][k] * metric_vars.lapse;
                                    
                                       

                FOR1(l)
                {  
                    FOR1(m)
                    {
                        primaryVector[i] +=  -  (2.0 * metric_vars.d1_lapse[l] * metric_vars.d1_lapse[m] * vars.fbar[i] * metric_vars.gamma[j][k] * gamma_UU[j][l] * gamma_UU[k][m])/metric_vars.lapse
                                                -  metric_vars.K_tensor[j][k] * metric_vars.K_tensor[l][m] * vars.fbar[i] * gamma_UU[j][l] * gamma_UU[k][m] * metric_vars.lapse
                                                + metric_vars.K_tensor[l][j] * metric_vars.K_tensor[m][k] * vars.fbar[i] * gamma_UU[j][l] * gamma_UU[k][m] * metric_vars.lapse
                                            + D1_K_tensor[i][j][k] * metric_vars.gamma[l][m] * gamma_UU[j][l] * gamma_UU[k][m] * vars.fhat * metric_vars.lapse 
                                            + D1_fspatial[j][m][k] * metric_vars.K_tensor[i][l] * gamma_UU[j][l] * gamma_UU[k][m] * lapse2
                                             - D1_fspatial[i][j][k] * metric_vars.K_tensor[l][m] * gamma_UU[j][l] * gamma_UU[k][m] * lapse2
                                             + D1_K_tensor[j][m][k] * gamma_UU[j][l] * gamma_UU[k][m] * vars.fspatial[i][l] * lapse2
                                             - D1_K_tensor[m][k][j] * gamma_UU[j][l] * gamma_UU[k][m] * vars.fspatial[i][l] * lapse2
                                              + D1_K_tensor[i][j][k] * gamma_UU[j][l] * gamma_UU[k][m] * vars.fspatial[l][m] * lapse2 
                                              -  D1_K_tensor[j][k][i] * gamma_UU[j][l] * gamma_UU[k][m] * vars.fspatial[l][m] * lapse2
                                              -  D1_fbar[m][k] * metric_vars.K_tensor[i][j] * gamma_UU[j][l] * gamma_UU[k][m] * shift_L[l]
                                               - D1_fbar[m][j] * metric_vars.K_tensor[i][k] * gamma_UU[j][l] * gamma_UU[k][m] * shift_L[l]
                                               -  D1_fbar[i][m] * metric_vars.K_tensor[j][k] * gamma_UU[j][l] * gamma_UU[k][m] * shift_L[l]
                                                -  D1_fbar[i][j] * metric_vars.K_tensor[m][k] * gamma_UU[j][l] * gamma_UU[k][m] * shift_L[l]
                                                -  D1_K_tensor[m][k][j] * vars.fbar[i] * gamma_UU[j][l] * gamma_UU[k][m] * shift_L[l]
                                                + metric_vars.K_tensor[i][m] * metric_vars.K_tensor[j][k] * gamma_UU[j][l] * gamma_UU[k][m] * vars.fhat * shift_L[l]
                                                + metric_vars.K_tensor[i][j] * metric_vars.K_tensor[m][k] * gamma_UU[j][l] * gamma_UU[k][m] * vars.fhat * shift_L[l]
                                                + (metric_vars.d1_lapse[m] * metric_vars.K_tensor[j][k] * vars.fbar[i] * gamma_UU[j][l] * gamma_UU[k][m] * shift_L[l])/metric_vars.lapse
                                                + (metric_vars.d1_lapse[j] * metric_vars.K_tensor[m][k] * vars.fbar[i] * gamma_UU[j][l] * gamma_UU[k][m] * shift_L[l])/metric_vars.lapse
                                                - D2_fspatial[i][m][j][k] * gamma_UU[j][l] * gamma_UU[k][m] * metric_vars.lapse * shift_L[l]
                                                + vars.v[i][m] * metric_vars.K_tensor[j][k] * gamma_UU[j][l] * gamma_UU[k][m] * metric_vars.lapse * shift_L[l] 
                                                -  D1_K_tensor[i][j][k] * vars.fbar[l] * gamma_UU[j][l] * gamma_UU[k][m] * shift_L[m]
                                                + (metric_vars.d1_lapse[k] * metric_vars.K_tensor[i][j] * vars.fbar[l] * gamma_UU[j][l] * gamma_UU[k][m] * shift_L[m])/metric_vars.lapse
                                                + (metric_vars.d1_lapse[j] * metric_vars.K_tensor[i][k] * vars.fbar[l] * gamma_UU[j][l] * gamma_UU[k][m] * shift_L[m])/metric_vars.lapse
                                                + (2.0 * metric_vars.d1_lapse[j] * metric_vars.d1_lapse[k] * vars.fbar[i] * gamma_UU[j][l] * gamma_UU[k][m] * shift_L[l] * shift_L[m])/lapse3
                                                -  (2.0 * D1_fbar[i][j] * metric_vars.d1_lapse[k] * gamma_UU[j][l] * gamma_UU[k][m] * shift_L[l] * shift_L[m])/lapse2
                                                -  (D2_lapse[k][j] * vars.fbar[i] * gamma_UU[j][l] * gamma_UU[k][m] * shift_L[l] * shift_L[m])/lapse2
                                                + (2.0 * metric_vars.d1_lapse[j] * metric_vars.K_tensor[i][k] * gamma_UU[j][l] * gamma_UU[k][m] * vars.fhat * shift_L[l] * shift_L[m])/lapse2
                                                + (D2_fbar[i][k][j] * gamma_UU[j][l] * gamma_UU[k][m] * shift_L[l] * shift_L[m])/metric_vars.lapse
                                                -  (2.0 * d1.fhat[j] * metric_vars.K_tensor[i][k] * gamma_UU[j][l] * gamma_UU[k][m] * shift_L[l] * shift_L[m])/metric_vars.lapse
                                                 - (D1_K_tensor[i][j][k] * gamma_UU[j][l] * gamma_UU[k][m] * vars.fhat * shift_L[l] * shift_L[m])/metric_vars.lapse;
                                                 

                                                
                        FOR1(n)
                        {
                            
                            FOR1(o)
                            {
                                primaryVector[i] += - D1_K_tensor[j][l][k] * metric_vars.gamma[m][n] * gamma_UU[j][m] * gamma_UU[k][n] * gamma_UU[l][o] * vars.fspatial[i][o] * lapse2
                                                     - metric_vars.K_tensor[j][n] * metric_vars.K_tensor[k][o] * gamma_UU[j][m] * gamma_UU[k][n] * gamma_UU[l][o] * vars.fspatial[i][l] * metric_vars.lapse * shift_L[m]
                                                     - metric_vars.K_tensor[i][j] * metric_vars.K_tensor[n][o] * gamma_UU[j][m] * gamma_UU[k][n] * gamma_UU[l][o] * vars.fspatial[k][l] * metric_vars.lapse * shift_L[m]
                                                     + D1_fspatial[i][o][l] * metric_vars.K_tensor[j][k] * gamma_UU[j][m] * gamma_UU[k][n] * gamma_UU[l][o] * shift_L[m] * shift_L[n]
                                                      + 2.0 * D1_fspatial[i][l][k] * metric_vars.K_tensor[j][o] * gamma_UU[j][m] * gamma_UU[k][n] * gamma_UU[l][o] * shift_L[m] * shift_L[n]
                                                      + D1_K_tensor[j][l][k] * gamma_UU[j][m] * gamma_UU[k][n] * gamma_UU[l][o] * vars.fspatial[i][o] * shift_L[m] * shift_L[n]
                                                       + (metric_vars.K_tensor[j][o] * metric_vars.K_tensor[k][l] * vars.fbar[i] * gamma_UU[j][m] * gamma_UU[k][n] * gamma_UU[l][o] * shift_L[m] * shift_L[n])/ metric_vars.lapse
                                                       + (metric_vars.K_tensor[j][k] * metric_vars.K_tensor[o][l] * vars.fbar[i] * gamma_UU[j][m] * gamma_UU[k][n] * gamma_UU[l][o] * shift_L[m] * shift_L[m])/metric_vars.lapse
                                                       + (3.0 * metric_vars.K_tensor[i][k] * metric_vars.K_tensor[j][l] * vars.fbar[m] * gamma_UU[j][m] * gamma_UU[k][n] * gamma_UU[l][o] * shift_L[n] * shift_L[o])/metric_vars.lapse
                                                       + (metric_vars.K_tensor[i][j] * metric_vars.K_tensor[k][l] * vars.fbar[m] * gamma_UU[j][m] * gamma_UU[k][n] * gamma_UU[l][o] * shift_L[n] * shift_L[o])/ metric_vars.lapse;
                            }
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

        current_cell.store_vars(field_amplitude, c_field_amplitude);
    }
};

#endif /* FIXEDBGDIAGNOSTICS_HPP_ */
