/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALCONDITIONS_HPP_
#define INITIALCONDITIONS_HPP_

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
//! Class which creates the initial conditions
class InitialConditions
{
  protected:
    const double m_dx;
    //const double m_amplitude_re, m_amplitude_im;
    //const double m_omega;
    const std::array<double, CH_SPACEDIM> m_center;
    const KerrSchildFixedBG::params_t m_bg_params;
    const double m_tensor_mass;
    const double m_initial_constant;

    //load in Vars from the field
    // The evolution vars
  
    template <class data_t>
    using Vars = FixedBGTensorField<TensorPotential>::template Vars<data_t>;
    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;



  public:
    //! The constructor for the class
    //const double a_amplitude_re, const double a_amplitude_im, const double a_omega,
    InitialConditions(const double tensor_mass, const std::array<double, CH_SPACEDIM> a_center,
                      const KerrSchildFixedBG::params_t a_bg_params,
                      const double a_dx, const double a_initial_constant)//, const double a_fhat, const Tensor<1,data_t> a_fbar, const Tensor<2,data_t> a_fspatial)
        : m_dx(a_dx), m_center(a_center), m_bg_params(a_bg_params), m_tensor_mass(tensor_mass),
        m_initial_constant(a_initial_constant)
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

        //const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        //Populate the field variables with their initial conditions
        //data_t fhat = m_fhat; 
        //Tensor<1, data_t> fbar = m_fbar;
        //Tensor<2, data_t> fspatial = m_fspatial;
        const double initial_constant = m_initial_constant; 
        //const mass = tensor_mass; 
        //data_t phi_Re = m_amplitude_re;
        //data_t Pi_Im = m_amplitude_im * m_omega;

        Vars<data_t> vars;
        VarsTools::assign(vars,0.);
        Tensor<2,data_t> fspatial; //Spatial component of the tensor field
        Tensor<1,data_t> fbar; //Half-projected component of the tensor field
        data_t fhat; //Scalar part of the tensor field 
        //Conjugate components
        //Tensor<2,data_t> u; // Spatial rank 2 u field
        Tensor<2,data_t> v; //Spatial rank 2 v field
        //Tensor<1,data_t> p; //Spatial rank 1 p field
        Tensor<1,data_t> q; //Spatial rank 1 q field
        data_t w; //Scalar component



        /*
        data_t rad = coords.get_radius();
        data_t value = 200*(exp(-rad) * sin(-rad));
        const double frequency = 2 * M_PI /128.0 ;

        data_t factor = 1.0 * exp(-pow((coords.x - 5.0),2)/2 /2. /2.);
        factor = 1;
        data_t amplitude = factor * cos( - frequency * coords.x);
        
        data_t momentum = -frequency * factor * sin(-frequency * coords.x); 
        */
        vars.fhat = 0.0;
       

        FOR(i)
        {
          vars.fbar[i] = 0.0;
          vars.q[i] = 0.0;

          FOR(j)
          {
            vars.fspatial[i][j] = 0.0;
            vars.v[i][j] = 0.0;            
          } 
        }
        vars.fspatial[0][0] = 1.0;
        vars.fspatial[1][1] = -1.0;

        /*
        vars.fspatial[0][0] = initial_constant * amplitude;
        vars.fspatial[1][1] = -initial_constant * amplitude;
        vars.fspatial[0][1] = initial_constant * amplitude;
        vars.fspatial[1][0] = initial_constant * amplitude;

        vars.v[0][0] = initial_constant * momentum;
        vars.v[1][1] = -initial_constant * momentum;
        vars.v[0][1] = initial_constant * momentum;
        vars.v[1][0] = initial_constant * momentum;
        */
 
        //Traceless condition
        vars.fhat = TensorAlgebra::compute_trace(gamma_UU, vars.fspatial);

        //(Riemann constraint?) R = 0 for vacuum
        
        vars.fbar[0] = 0.0;
        FOR2(i,j)
        {
          FOR2(k,l)
          {
            vars.fbar[0] += -gamma_UU[i][k] * gamma_UU[j][l] * metric_vars.K_tensor[i][j] * metric_vars.K_tensor[k][l];
          }
        }
        vars.fbar[0] += metric_vars.K * metric_vars.K;
        vars.fbar[0] *= vars.fhat;
        vars.fbar[0] /= metric_vars.d1_K[0];

        //Spatial projection of Lorentz condition
        FOR1(i)
        {
          vars.q[i] = -metric_vars.K * vars.fbar[i];
          FOR2(j,k)
          {
            vars.q[i] += -gamma_UU[j][k] * metric_vars.K_tensor[i][k] * vars.fbar[j]; //ADD IN DERIVATIVES OF fspatial LATER
          }
        }
        //Normal projection of Lorentz condition
        vars.w = -metric_vars.K * vars.fhat;
        FOR2(i,j)
        {
          vars.w += 0.0;//DERIVATIVE OF fbar
          FOR2(k,l)
          {
            vars.w += -metric_vars.K_tensor[i][j] * vars.fspatial[k][l] * gamma_UU[i][k] * gamma_UU[j][l]; 
          }
        }
        //Momentum traceless condition spatial projection of
        vars.v[0][0] = vars.w / gamma_UU[0][0];



        //vars.fhat = 0.0;
        /*
        FOR2(i,j)
        {
          vars.fhat += vars.fspatial[i][j] * gamma_UU[i][j];
        }
        */
        /*
        const auto chris_phys = TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU);
        FOR1(i)
        {
          vars.q[i] = metric_vars.K;
          FOR2(k,j)
          {
            vars.q[i] += gamma_UU[j][k] * d1.fspatial[i][j][k]
                      - gamma_UU[j][k] * metric_vars.K_tensor[i][k] * vars.fbar[j]; 
            FOR1(l)
            {
              vars.q[i] += -gamma_UU[j][k] * (chris_phys.ULL[l][i][k] * vars.fspatial[l][j] + chris_phys.ULL[l][j][k] * vars.fspatial[i][l]);  
            }
          }
        }
        */

        current_cell.store_vars(vars);
    }
};

#endif /* INITIALCONDITIONS_HPP_ */
