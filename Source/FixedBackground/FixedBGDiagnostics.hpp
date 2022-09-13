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

        // Defined for convenience 
        Tensor<2, data_t> i_u; 
        FOR2(i,j)
        {
            i_u[i][j] = d1.fbar[j][i] - vars.fhat * metric_vars.K_tensor[i][j];

            FOR1(k)
            {
                i_u[i][j] += -chris_phys.ULL[k][j][i] * vars.fbar[k]; 

                FOR1(l)
                {
                    i_u[i][j] += -gamma_UU[l][k] * vars.fspatial[j][k] * metric_vars.K_tensor[i][l];
                }
            }

        }

        Tensor<3, data_t> d1_i_u; 
        FOR3(i,j,k)
        {
            d1_i_u[i][j][k] = -d1.fhat[k] * metric_vars.K_tensor[i][j] - vars.fhat * metric_vars.d1_K_tensor[i][j][k] + d2.fbar[j][i][k]; 

            FOR1(l)
            {
                d1_i_u[i][j][k] += -chris_phys.ULL[l][j][i] * d1.fbar[l][k] - metric_vars.d1_chris_phys[l][j][i][k] * vars.fbar[l];

                FOR1(m)
                {
                    d1_i_u[i][j][k] += -metric_vars.d1_gamma_UU[l][m][k] * vars.fspatial[j][m] * metric_vars.K_tensor[i][l] - gamma_UU[l][m] * d1.fspatial[j][m][k] * metric_vars.K_tensor[i][l]
                    - gamma_UU[l][m] * vars.fspatial[j][m] * metric_vars.d1_K_tensor[i][l][k];
                }
            }
            
        }

        Tensor<1, data_t> i_p;
        FOR1(i)
        {
            i_p[i] = d1.fhat[i]; 

            FOR2(j,k)
            {
                i_p[i] += -2.0 * gamma_UU[j][k] * vars.fbar[j] * metric_vars.K_tensor[i][k];
            }
            
        }
        Tensor<2, data_t> d1_i_p;
        FOR2(i,j)
        {
            d1_i_p[i][j] = d2.fhat[i][j]; 

            FOR2(k,l)
            {
                d1_i_p[i][j] += -2.0 * (metric_vars.d1_gamma_UU[l][k][j] * vars.fbar[l] * metric_vars.K_tensor[i][k] 
                + gamma_UU[l][k] * d1.fbar[l][j] * metric_vars.K_tensor[i][k] + gamma_UU[l][k] * vars.fbar[l] * metric_vars.d1_K_tensor[i][k][j]); 
            }
            
        }
        
        Tensor<3, data_t> cd1_K_tensor;
        FOR3(i,j,k)
        {
            cd1_K_tensor[i][j][k] = 0.0;//metric_vars.d1_K_tensor[i][j][k];
            FOR1(l)
            {
                cd1_K_tensor[i][j][k] += 0.0;//-chris_phys.ULL[l][k][i] * metric_vars.K_tensor[l][j] -chris_phys.ULL[l][k][j] * metric_vars.K_tensor[i][l];
            }

        }

        //TRACE DIAGNOSTICS
        data_t trace_field = 0.0;

        FOR2(i,j)
        {
            trace_field += gamma_UU[i][j] * vars.fspatial[i][j]; 
        }
        trace_field += -vars.fhat;
        
        Tensor<1, data_t> trace_of_F;

        FOR1(i)
        {
            trace_of_F[i] = -d1.fhat[i];
            FOR2(j,k)
            {
                trace_of_F[i] += 2.0 * gamma_UU[j][k] * metric_vars.K_tensor[i][j] * vars.fbar[k];
            }
        }


        data_t trace_momentum = 0.0;

        FOR2(i,j)
        {
            trace_momentum += gamma_UU[i][j] * vars.v[i][j];
        }
        trace_momentum += -vars.w;
   
        data_t transverseScalar = 0.0;
        
        transverseScalar = -metric_vars.K * vars.fhat - vars.w;
        FOR2(i,j)
        {
          transverseScalar += gamma_UU[i][j] * d1.fbar[i][j];
          FOR1(k)
          {
            transverseScalar += -gamma_UU[i][j] * chris_phys.ULL[k][j][i] * vars.fbar[k];
            FOR1(l)
            {
              transverseScalar += -metric_vars.K_tensor[i][j] * vars.fspatial[k][l] * gamma_UU[i][k] * gamma_UU[j][l]; 
            }
          }
        }
       
        Tensor<1, data_t> transverseVector;
        
        FOR1(i)
        {
          transverseVector[i] = -vars.q[i];//-metric_vars.K * vars.fbar[i] - vars.q[i];
          FOR2(j,k)
          {
            //transverseVector[i] += -gamma_UU[j][k] * metric_vars.K_tensor[i][k] * vars.fbar[j];
            transverseVector[i] += gamma_UU[j][k] * d1.fspatial[j][i][k];
            FOR1(l)
            {
              transverseVector[i] += -gamma_UU[j][k] * (chris_phys.ULL[l][k][i] * vars.fspatial[j][l] + chris_phys.ULL[l][k][j] * vars.fspatial[l][i]);
            }
          }
        }
        //transverseVector[2] = 0.0;//vars.q[2] -gamma_UU[1][1] * d1.fspatial[0][1][1];//- (gamma_UU[0][0] * d1.fspatial[0][0][0] + gamma_UU[1][1] * d1.fspatial[0][1][1] + gamma_UU[2][2] * d1.fspatial[0][2][2]);
        //transverseVector[2] = vars.q[2] - d1.fspatial[2][2][0];
        
        //NEEDS MASS TERM ADDING
        const double temp_mass = 0.0;
        
        data_t primaryScalar = -temp_mass * temp_mass * vars.fhat;
        
        FOR2(i,j)
        {
          primaryScalar += -gamma_UU[i][j] * (d1.q[i][j] - d1_i_p[i][j]);

          FOR1(k)
          {
            primaryScalar += gamma_UU[i][j] * (chris_phys.ULL[k][i][j] * vars.q[k] - chris_phys.ULL[k][i][j] * i_p[k]);

            FOR1(l)
            {
              primaryScalar += gamma_UU[i][k] * gamma_UU[j][l] * metric_vars.ricci_phys[i][j] * vars.fspatial[k][l];
            }
          }
        }        
        /*
        FOR1(i)
        {
          primaryScalar += vars.w * metric_vars.shift[i] * metric_vars.d1_lapse[i] / metric_vars.lapse;
          FOR1(j)
          {
            primaryScalar += gamma_UU[i][j] *(i_p[j] - vars.q[j])* metric_vars.d1_lapse[i];
            FOR2(k,l)
            {

            
              primaryScalar += -2.0 * metric_vars.lapse * gamma_UU[i][k] * gamma_UU[j][l] * (metric_vars.K * metric_vars.K_tensor[i][j] + metric_vars.ricci_phys[i][j]) * vars.fspatial[k][l]
                                + metric_vars.lapse * gamma_UU[i][k] * gamma_UU[j][l] * i_u[i][j] * metric_vars.K_tensor[k][l];
                                
             
              
              FOR2(m,n)
              {
                primaryScalar += metric_vars.lapse * gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * (metric_vars.K_tensor[i][m] * metric_vars.K_tensor[j][k] * vars.fspatial[l][n]);
              }
            }
          }
        }
        */
        Tensor<1, data_t> primaryVector;
        FOR1(i)
        {
          primaryVector[i] = -temp_mass * temp_mass * vars.fbar[i];

          FOR2(j,k)
          {
            primaryVector[i] += -gamma_UU[j][k] * (d1.v[j][i][k] - d1_i_u[j][i][k]) 
                                + gamma_UU[j][k] * metric_vars.ricci_phys[i][j] * vars.fbar[k];

            FOR1(l)
            {
              primaryVector[i] += gamma_UU[j][k] * (chris_phys.ULL[l][k][i] * vars.v[j][l] + chris_phys.ULL[l][k][j] * vars.v[l][i]
                                                   -chris_phys.ULL[l][k][i] * i_u[j][l]    - chris_phys.ULL[l][k][j] * i_u[l][i]);
            } 
          }
        }        


        /*
        FOR1(i)
        {
          primaryVector[i] = vars.w * metric_vars.d1_lapse[i];

          FOR1(j)
          {
            
            FOR1(k)
            {
              primaryVector[i] += 2.0 * metric_vars.lapse * gamma_UU[j][k] * metric_vars.ricci_phys[i][j] * vars.fbar[k];
                                //- gamma_UU[j][k] * metric_vars.K_tensor[k][i] * i_p[j]

                               // -2.0*gamma_UU[j][k] * metric_vars.K_tensor[j][i] * metric_vars.K * vars.fbar[k]

                                //+2.0 * gamma_UU[j][k] * metric_vars.ricci_phys[i][j] * vars.fbar[k];
              
              FOR2(l,m)
              {
                primaryVector[i] += 0.0;//-2.0*gamma_UU[j][l] * gamma_UU[k][m] * (cd1_K_tensor[j][k][i] - cd1_K_tensor[i][j][k]) * vars.fspatial[l][m]
                                    //+2.0*gamma_UU[j][l] * gamma_UU[k][m] * metric_vars.K_tensor[k][j] * metric_vars.K_tensor[i][m] * vars.fbar[l];

              }
            }
          }
        }
        */
        //Store diagnostic variables if outside the event horizon
        //const double horizon = 0.0;
        const double horizon = 0.8/2.0;
        //const double horizon = 0.1;
        const double xx = coords.x * coords.x;
        const double yy = coords.y * coords.y;
        const double zz = coords.z * coords.z;
        const double rr = sqrt(xx + yy + zz);


        /*
        simd_conditional((simd_compare_gt(r,horizon), 1)) 
        {
         //do stuff 
        }
        */
        //const doubledouble_ra d = sqrt(coords.x * coords.x + coords.y * coords.y + coords.z * coords.z);
        if (rr > horizon)//((simd_compare_gt(r,horizon)))
        {
      
          current_cell.store_vars(trace_field, c_trace_field);  
          current_cell.store_vars(trace_momentum, c_trace_momentum); 




          current_cell.store_vars(transverseScalar, c_transverseScalar);
          //current_cell.store_vars(metric_vars.lapse, c_transverseScalar);  

          current_cell.store_vars(transverseVector[0], c_transverseVector1);
          current_cell.store_vars(transverseVector[1], c_transverseVector2); 
          current_cell.store_vars(transverseVector[2], c_transverseVector3);  
          //current_cell.store_vars(metric_vars.d1_lapse[0], c_transverseVector1);
          //current_cell.store_vars(metric_vars.d1_lapse[1], c_transverseVector2); 
          //current_cell.store_vars(metric_vars.d1_lapse[2], c_transverseVector3);     

          current_cell.store_vars(primaryScalar, c_primaryConstraintScalar);
          current_cell.store_vars(primaryVector[0], c_primaryConstraintVector1);
          current_cell.store_vars(primaryVector[1], c_primaryConstraintVector2);
          current_cell.store_vars(primaryVector[2], c_primaryConstraintVector3);
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
        }
        
    
      

    }
};

#endif /* FIXEDBGDIAGNOSTICS_HPP_ */
