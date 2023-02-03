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

  public:

  
    FixedBGDiagnostics(matter_t a_matter, background_t a_background, double a_dx,
                     std::array<double, CH_SPACEDIM> a_center)
        : m_matter(a_matter), m_deriv(a_dx), m_dx(a_dx),
          m_background(a_background), m_center(a_center)
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
        const auto gamma_UU = TensorAlgebra::compute_inverse_sym(metric_vars.gamma);
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
 
        const data_t costheta = z/r;
        const data_t sintheta = rho/r;

        const data_t sinphi = y/rho;
        const data_t cosphi = x/rho;

        const data_t sin2phi = 2.0 * sinphi * cosphi; 
        const data_t cos2phi = cosphi*cosphi - sinphi*sinphi;

        //TRACE DIAGNOSTICS
        data_t trace_field = 0.0;

        FOR2(i,j)
        {
            trace_field += gamma_UU[i][j] * vars.fspatial[i][j]; 
        }

        data_t trace_momentum = 0.0;

        FOR2(i,j)
        {
            trace_momentum += gamma_UU[i][j] * vars.v[i][j];
        }
   
        data_t primaryScalar = 0.0;
    
        FOR2(i,j)
        {
            FOR1(k)
            {
                FOR1(l)
                {
                    primaryScalar += gamma_UU[i][k] * gamma_UU[j][l] * (metric_vars.ricci_phys[i][j] * vars.fspatial[k][l] 
                                                                        + metric_vars.d1_lapse[k] * d1.fspatial[i][l][j] / metric_vars.lapse
                                                                        - metric_vars.d1_lapse[k] * metric_vars.d1_lapse[l] * vars.fspatial[i][j] / metric_vars.lapse / metric_vars.lapse
                                                                        + metric_vars.d2_lapse[k][l] * vars.fspatial[i][j] / metric_vars.lapse);
                    FOR1(m)
                    {
                        primaryScalar += - gamma_UU[i][k] * gamma_UU[j][l] * (metric_vars.d1_lapse[k] * (chris_phys.ULL[m][i][j]  * vars.fspatial[m][l] + chris_phys.ULL[m][l][j] * vars.fspatial[i][m]) / metric_vars.lapse
                                        + vars.fspatial[i][j] * chris_phys.ULL[m][l][k] * metric_vars.d1_lapse[m] / metric_vars.lapse);
                    }
                }
            }
        } 
        
        Tensor<1, data_t> primaryVector;
        FOR1(i)
        {
            primaryVector[i] = 0.0;

            FOR2(j,k)
            {
                primaryVector[i] += gamma_UU[j][k] * d1.v[j][i][k];

                FOR1(l)
                {
                    primaryVector[i] += -gamma_UU[j][k] * (chris_phys.ULL[l][k][i] * vars.v[j][l] + chris_phys.ULL[l][k][j] * vars.v[l][i]);

                } 
            }
        }     
        Tensor<1, data_t> transverseVector;

        FOR1(i)
        {
            transverseVector[i] = 0.0;
            FOR2(j,k)
            {
            transverseVector[i] += gamma_UU[j][k] * (metric_vars.lapse * d1.fspatial[j][i][k] + vars.fspatial[i][j] * metric_vars.d1_lapse[k]);
                FOR1(l)
                {
                    transverseVector[i] += -gamma_UU[j][k] * ( metric_vars.lapse * chris_phys.ULL[l][k][i] * vars.fspatial[j][l] + metric_vars.lapse * chris_phys.ULL[l][k][j] * vars.fspatial[l][i]);
                }
            }
        }

        //Effective energy density of the stress energy tensor - just the mass term
        const double mass = 0.1;
        data_t rho_eff;
        rho_eff = 0.0;
        FOR3(i,j,k)
        {
            FOR1(l)
            {
                rho_eff += 0.25 * mass * mass * gamma_UU[i][k] * gamma_UU[j][l] * vars.fspatial[i][j] * vars.fspatial[k][l]; 
            }
        }




        //Store diagnostic variables if outside the event horizon
        //const double horizon = 0.0;//0.8/2.0;
        //const double horizon = 0.8/2.0;
        const double horizon = 1.0/2.0;
      
        const double xx = coords.x * coords.x;
        const double yy = coords.y * coords.y;
        const double zz = coords.z * coords.z;
        const double rr = sqrt(xx + yy + zz);

        if (rr > horizon * 1.50)//((simd_compare_gt(r,horizon)))
        {
      
          current_cell.store_vars(trace_field, c_trace_field);  
          current_cell.store_vars(trace_momentum, c_trace_momentum); 

          current_cell.store_vars(0.0, c_transverseScalar);

          current_cell.store_vars(transverseVector[0], c_transverseVector1);
          current_cell.store_vars(transverseVector[1], c_transverseVector2); 
          current_cell.store_vars(transverseVector[2], c_transverseVector3);  

          current_cell.store_vars(primaryScalar, c_primaryConstraintScalar);

          current_cell.store_vars(primaryVector[0], c_primaryConstraintVector1);
          current_cell.store_vars(primaryVector[1], c_primaryConstraintVector2);
          current_cell.store_vars(primaryVector[2], c_primaryConstraintVector3);

          current_cell.store_vars(rho_eff, c_rho_eff);
        }
        else
        {
          current_cell.store_vars(0.0, c_trace_field);  
          current_cell.store_vars(0.0, c_trace_momentum); 

          current_cell.store_vars(0.0, c_transverseScalar); 

          current_cell.store_vars(0.0, c_transverseVector1);
          current_cell.store_vars(0.0, c_transverseVector2); 
          current_cell.store_vars(0.0, c_transverseVector3);  

          current_cell.store_vars(0.0, c_primaryConstraintScalar);

          current_cell.store_vars(0.0, c_primaryConstraintVector1);
          current_cell.store_vars(0.0, c_primaryConstraintVector2);
          current_cell.store_vars(0.0, c_primaryConstraintVector3);

          current_cell.store_vars(0.0, c_rho_eff);
        }
        
    
      

    }
};

#endif /* FIXEDBGDIAGNOSTICS_HPP_ */
