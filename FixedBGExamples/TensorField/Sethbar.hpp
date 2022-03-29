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
#include "KerrSchildFixedBG.hpp"
#include "Tensor.hpp"
#include "TensorPotential.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "FourthOrderDerivatives.hpp"
//! Class which creates the initial constraints
class Sethbar
{
  protected:
    const double m_dx;
    //const double m_amplitude_re, m_amplitude_im;
    //const double m_omega;
    const std::array<double, CH_SPACEDIM> m_center;
    const KerrSchildFixedBG::params_t m_bg_params;
    const double m_tensor_mass;
    const double m_initial_constant;
    const FourthOrderDerivatives m_deriv;

    //load in Vars from the field
    // The evolution vars
  
    template <class data_t>
    using Vars = FixedBGTensorField<TensorPotential>::template Vars<data_t>;
    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;



  public:
    //! The constructor for the class
    //const double a_amplitude_re, const double a_amplitude_im, const double a_omega,
    Sethbar(const double tensor_mass, const std::array<double, CH_SPACEDIM> a_center,
                      const KerrSchildFixedBG::params_t a_bg_params,
                      const double a_dx, const double a_initial_constant)//, const double a_fhat, const Tensor<1,data_t> a_fbar, const Tensor<2,data_t> a_fspatial)
        : m_dx(a_dx), m_center(a_center), m_bg_params(a_bg_params), m_tensor_mass(tensor_mass),
        m_initial_constant(a_initial_constant), m_deriv(a_dx)
        //m_fhat(a_fhat), m_fbar(a_fbar), m_fspatial(a_fspatial)
        //, m_amplitude_re(a_amplitude_re),   
        //m_amplitude_im(a_amplitude_im)
        //m_omega(a_omega),
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // get the metric vars
        KerrSchildFixedBG kerr_bh(m_bg_params, m_dx);
        MetricVars<data_t> metric_vars;
        kerr_bh.compute_metric_background(metric_vars, current_cell);
        const data_t det_gamma =
            TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
        const auto gamma_UU = TensorAlgebra::compute_inverse_sym(metric_vars.gamma);
        //const auto chris_phys = TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU);

        
        //Populate the field variables with their initial constraints
        //data_t fhat = m_fhat; 
        //Tensor<1, data_t> fbar = m_fbar;
        //Tensor<2, data_t> fspatial = m_fspatial;
        //const double initial_constant = m_initial_constant; 
        //const mass = tensor_mass; 
        //data_t phi_Re = m_amplitude_re;
        //data_t Pi_Im = m_amplitude_im * m_omega;
        //auto vars = current_cell.template load_vars<Vars>();
        Vars<data_t> vars;
        //VarsTools::assign(vars,0.);
        //Tensor<2,data_t> fspatial; //Spatial component of the tensor field
        //Tensor<1,data_t> fbar; //Half-projected component of the tensor field
        //data_t fhat; //Scalar part of the tensor field 
        //Conjugate components
        //Tensor<2,data_t> u; // Spatial rank 2 u field
        //Tensor<2,data_t> v; //Spatial rank 2 v field
        //Tensor<1,data_t> p; //Spatial rank 1 p field
        //Tensor<1,data_t> q; //Spatial rank 1 q field
        //data_t w; //Scalar component

        Tensor<3, data_t> chris_local; 
        FOR3(i,j,k)
        {
          chris_local[i][j][k] = 0.0;
          FOR1(l)
          {
            chris_local[i][j][k] += 0.5 * gamma_UU[i][l] * (metric_vars.d1_gamma[l][k][j] + metric_vars.d1_gamma[j][l][k] - metric_vars.d1_gamma[j][k][l]);
          }
          
        }


        Tensor<3, data_t> cd1_K_tensor;
        FOR3(i,j,k)
        {
            cd1_K_tensor[i][j][k] = metric_vars.d1_K_tensor[i][j][k];
            FOR1(l)
            {
                cd1_K_tensor[i][j][k] += -chris_local[l][k][i] * metric_vars.K_tensor[l][j] -chris_local[l][k][j] * metric_vars.K_tensor[i][l];
            }
        }
        
        const auto local_vars = current_cell.template load_vars<Vars>();

        vars.fhat = local_vars.fhat;
        FOR1(i)
        {
          vars.fbar[i] = local_vars.fbar[i];
          FOR1(j)
          {
            vars.fspatial[i][j] = local_vars.fspatial[i][j];
          }
        }
        
        //Calculate the derivatives
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    
       
        /*
        const auto K_inverse_UU = TensorAlgebra::compute_inverse_sym(metric_vars.K_tensor);
        FOR1(i)
        {
          vars.fbar[i] = 0.0;
          FOR2(j,k)
          {
            vars.fbar[i] += 0.5 * metric_vars.gamma[i][k] * K_inverse_UU[j][k] * d1.fhat[j];
          }
        }
        */
        double denom = 0.0;

        FOR1(i)
        {
          denom += -gamma_UU[0][i] * metric_vars.d1_K[i]; 
          FOR2(j,k)
          {
            denom += gamma_UU[0][i] * gamma_UU[j][k] * cd1_K_tensor[i][j][k];
          }
        }

        //vars.fbar[0] = 0.0;
        vars.fbar[0] = metric_vars.K * metric_vars.K * vars.fhat;
        //vars.fbar[0] = vars.fhat;
        
        FOR2(i,j)
        {
          vars.fbar[0] += gamma_UU[i][j] * vars.fhat * metric_vars.ricci_phys[i][j];
          FOR2(k,l)
          {
            vars.fbar[0] += -gamma_UU[k][j] * metric_vars.K_tensor[i][j] * gamma_UU[i][l] * metric_vars.K_tensor[l][k] * vars.fhat;
          }
        }
        vars.fbar[0] = 0.0;
        double radius = coords.get_radius();

        //vars.fbar[0] = coords.x;
        /*
        FOR1(i)
        {
          vars.fbar[i] = 0;
          FOR1(j)
          {
            vars.fbar[i] += -metric_vars.shift[j] * vars.fspatial[i][j] / metric_vars.lapse;
          }
        }
        */
        //vars.fbar[0] = 1.0;
        //vars.fbar[1] = 1.0;
        //vars.fbar[2] = - (metric_vars.shift[0] * vars.fbar[0] + metric_vars.shift[1] * vars.fbar[1]) / metric_vars.shift[2]; 
        //vars.fbar[0] = 10 * exp(-radius * radius);
        //vars.fbar[0] /= (2.0 * denom);
        /*
        if (-1e-10 < denom < 0)
        {
          denom - -1e-10;
        }
        else if (0 < denom < 1e-10)
        {
          denom = 1e-10;
        }
        */

        //d1fhat -> D1 fspatial^i_i -> gamma_UU D1 f_ _ -> d1 + christoffel * 2  -> christoffel * 2  
        //make sure if its d1 causing the problems 
        /*
        data_t denom = 2.0 * (-2.0 * metric_vars.K_tensor[1][0] * metric_vars.K_tensor[1][2] * metric_vars.K_tensor[2][0] + metric_vars.K_tensor[1][1] * metric_vars.K_tensor[2][0] * metric_vars.K_tensor[2][0]
                                    + metric_vars.K_tensor[1][0] * metric_vars.K_tensor[1][0] * metric_vars.K_tensor[2][2] + metric_vars.K_tensor[1][2] * metric_vars.K_tensor[1][2] * metric_vars.K_tensor[0][0]
                                    - metric_vars.K_tensor[1][1] * metric_vars.K_tensor[2][2] * metric_vars.K_tensor[0][0]);



        vars.fbar[0] = metric_vars.K_tensor[1][2] * metric_vars.K_tensor[1][2] * d1.fhat[0] - metric_vars.K_tensor[1][1] * metric_vars.K_tensor[2][2] * d1.fhat[0]
                      -metric_vars.K_tensor[1][2] * metric_vars.K_tensor[2][0] * d1.fhat[1] + metric_vars.K_tensor[1][0] * metric_vars.K_tensor[2][2] * d1.fhat[1]
                      -metric_vars.K_tensor[1][0] * metric_vars.K_tensor[1][2] * d1.fhat[2] + metric_vars.K_tensor[1][1] * metric_vars.K_tensor[2][0] * d1.fhat[2];
        vars.fbar[0] /= denom;

        vars.fbar[1] = -metric_vars.K_tensor[1][2] * metric_vars.K_tensor[2][0] * d1.fhat[0] + metric_vars.K_tensor[1][0] * metric_vars.K_tensor[2][2] * d1.fhat[0]
                        +metric_vars.K_tensor[2][0] * metric_vars.K_tensor[2][0] * d1.fhat[1] - metric_vars.K_tensor[2][2] * metric_vars.K_tensor[0][0] * d1.fhat[1]
                        -metric_vars.K_tensor[1][0] * metric_vars.K_tensor[2][0] * d1.fhat[2] + metric_vars.K_tensor[1][2] * metric_vars.K_tensor[0][0] * d1.fhat[2];
        vars.fbar[1] /= denom;

        vars.fbar[2] = -metric_vars.K_tensor[1][0]  *metric_vars.K_tensor[1][2] * d1.fhat[0] + metric_vars.K_tensor[1][1] * metric_vars.K_tensor[2][0] * d1.fhat[0]
                        -metric_vars.K_tensor[1][0] * metric_vars.K_tensor[2][0] * d1.fhat[1] + metric_vars.K_tensor[1][2] * metric_vars.K_tensor[0][0] * d1.fhat[1]
                        +metric_vars.K_tensor[1][0] * metric_vars.K_tensor[1][0] * d1.fhat[2] - metric_vars.K_tensor[1][1] * metric_vars.K_tensor[0][0] * d1.fhat[2];
        vars.fbar[2] /= denom;
        */



        current_cell.store_vars(vars);
    }
};

#endif /* SETHBAR_HPP_ */
