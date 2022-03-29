/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDBGDENSITIES_HPP_
#define FIXEDBGDENSITIES_HPP_

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
template <class matter_t, class background_t> class FixedBGDensities
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
    FixedBGDensities(matter_t a_matter, background_t a_background, double a_dx,
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
        const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys =
            compute_christoffel(metric_vars.d1_gamma, gamma_UU);
        const emtensor_t<data_t> emtensor = m_matter.compute_emtensor(
            vars, metric_vars, d1, gamma_UU, chris_phys.ULL);
        const data_t det_gamma =
            TensorAlgebra::compute_determinant_sym(metric_vars.gamma);

        // first rho - note that this is the conserved rho
        // defined using the timelike Killing vector
        // not that of the Eulerian observers
        data_t rho = emtensor.rho * metric_vars.lapse;
        FOR1(i) { rho += -emtensor.Si[i] * metric_vars.shift[i]; }
        rho *= sqrt(det_gamma);

        // now rho J, also the conserved one
        Tensor<1, data_t> dxdphi;
        dxdphi[0] = -coords.y;
        dxdphi[1] = coords.x;
        dxdphi[2] = 0;

        data_t rhoJ = 0;
        FOR1(i) { rhoJ += -emtensor.Si[i] * dxdphi[i]; }
        rhoJ *= sqrt(det_gamma);

        // assign values of conserved density in output box,
        current_cell.store_vars(rho, c_rho);
        current_cell.store_vars(rhoJ, c_rhoJ);

        Tensor<3, data_t> cd1_K_tensor;
        FOR3(i,j,k)
        {
            cd1_K_tensor[i][j][k] = metric_vars.d1_K_tensor[i][j][k];
            FOR1(l)
            {
                cd1_K_tensor[i][j][k] += -chris_phys.ULL[l][k][i] * metric_vars.K_tensor[l][j] -chris_phys.ULL[l][k][j] * metric_vars.K_tensor[i][l];
            }
        }

        //TRACE DIAGNOSTICS
        data_t trace_of_field = 0.0;

        FOR2(i,j)
        {
            trace_of_field += gamma_UU[i][j] * vars.fspatial[i][j]; 
        }
        trace_of_field += -vars.fhat;
        pout()<<"trace here = " << trace_of_field << endl; 

        Tensor<1, data_t> trace_of_F;

        FOR1(i)
        {
            trace_of_F[i] = -d1.fhat[i];
            FOR2(j,k)
            {
                trace_of_F[i] += 2.0 * gamma_UU[j][k] * metric_vars.K_tensor[i][j] * vars.fbar[k];
            }
        }


        data_t trace_of_B = 0.0;

        FOR2(i,j)
        {
            trace_of_B += gamma_UU[i][j] * vars.v[i][j];
        }
        trace_of_B += -vars.w;
        Tensor<1, data_t> fhatderivs;

        FOR1(i)
        {
            fhatderivs[i] = d1.fhat[i];
        }

        data_t KinvK = 0.0;

        const auto K_inverse_UU = TensorAlgebra::compute_inverse_sym(metric_vars.K_tensor);
        FOR2(i,j)
        {
            KinvK += K_inverse_UU[i][j] * metric_vars.K_tensor[i][j];
        }

        /*
        data_t denom = 0.0;

        FOR1(i)
        {
          denom += -gamma_UU[0][i] * metric_vars.d1_K[i]; 
          FOR2(j,k)
          {
            denom += gamma_UU[0][i] * gamma_UU[j][k] * cd1_K_tensor[i][j][k];
          }
        }

        vars.fbar[0] = metric_vars.K * metric_vars.K * vars.fhat;

        FOR2(i,j)
        {
          vars.fbar[0] += gamma_UU[i][j] * vars.fhat * metric_vars.ricci_phys[i][j];
          FOR2(k,l)
          {
            vars.fbar[0] += gamma_UU[k][j] * metric_vars.K_tensor[i][j] * gamma_UU[i][l] * metric_vars.K_tensor[l][k] * vars.fhat;
          }
        }
        vars.fbar[0] /= (2.0 * denom);
        */        
        data_t lorenzCont = 0.0;
        lorenzCont = -metric_vars.K * vars.fhat - vars.w;
        FOR2(i,j)
        {
          lorenzCont += gamma_UU[i][j] * d1.fbar[i][j];
          FOR1(k)
          {
            lorenzCont += -gamma_UU[i][j] * chris_phys.ULL[k][j][i] * vars.fbar[k];
            FOR1(l)
            {
              lorenzCont += -metric_vars.K_tensor[i][j] * vars.fspatial[k][l] * gamma_UU[i][k] * gamma_UU[j][l]; 
            }
          }
        }

        Tensor<1, data_t> lorenzProj;
        FOR1(i)
        {
          lorenzProj[i] = -metric_vars.K * vars.fbar[i] - vars.q[i];
          FOR2(j,k)
          {
            lorenzProj[i] += -gamma_UU[j][k] * metric_vars.K_tensor[i][k] * vars.fbar[j];
            lorenzProj[i] += gamma_UU[j][k] * d1.fspatial[j][i][k];
            FOR1(l)
            {
              lorenzProj[i] += -gamma_UU[j][k] * (chris_phys.ULL[l][k][i] * vars.fspatial[l][j] + chris_phys.ULL[l][k][j] * vars.fspatial[i][l]);
            }
          }
        }

        data_t primaryScalar = 0.0;
        
        FOR2(i,j)
        {

          primaryScalar += d2.fspatial[i][i][j][j] - d2.fspatial[i][j][i][j];
          
        }

        Tensor<1, data_t> primaryVector;
        FOR1(i)
        {
          primaryVector[i] = 0;

          FOR1(j)
          {
            primaryVector[i] += - metric_vars.lapse * d2.fbar[i][j][j] - d1.v[j][j][i] + d1.v[i][j][i] + metric_vars.lapse * d2.fbar[j][j][i]; 
            
            FOR1(k)
            {
              primaryVector[i] += metric_vars.d1_shift[j][i] * d1.fspatial[k][k][i] + metric_vars.shift[j] * d2.fspatial[k][k][j][i]
              + 2.0 * d1.fspatial[j][k][i] * metric_vars.d1_shift[j][k] + vars.fspatial[j][k] * metric_vars.d2_shift[j][k][i]
              - metric_vars.d1_shift[j][k] * d1.fspatial[i][k][j] + metric_vars.shift[k] * d2.fspatial[i][k][j][k]
              + d1.fspatial[j][k][k] * metric_vars.d1_shift[j][i] + vars.fspatial[j][k] * metric_vars.d2_shift[j][i][k]
              + d1.fspatial[i][j][k] * metric_vars.d1_shift[j][k] + vars.fspatial[i][j] * metric_vars.d2_shift[j][k][k];
            }
            
          }
        }
        current_cell.store_vars(trace_of_field, c_trace_field);  
        current_cell.store_vars(trace_of_B, c_traceB); 
        current_cell.store_vars(fhatderivs[0], c_d1fhat1); 
        current_cell.store_vars(fhatderivs[1], c_d1fhat2); 
        current_cell.store_vars(fhatderivs[2], c_d1fhat3); 

        current_cell.store_vars(lorenzCont, c_lorenzCont); 
        current_cell.store_vars(lorenzProj[0], c_lorenzProj1);
        current_cell.store_vars(lorenzProj[1], c_lorenzProj2); 
        current_cell.store_vars(lorenzProj[2], c_lorenzProj3);  
        current_cell.store_vars(primaryScalar, c_primaryConstraintScalar);
        current_cell.store_vars(primaryVector[0], c_primaryConstraintVector1);
        current_cell.store_vars(primaryVector[1], c_primaryConstraintVector2);
        current_cell.store_vars(primaryVector[2], c_primaryConstraintVector3);
        current_cell.store_vars(KinvK, c_KinvK);


    }
};

#endif /* FIXEDBGDENSITIES_HPP_ */
