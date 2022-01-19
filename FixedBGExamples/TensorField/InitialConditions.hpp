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
                      const double a_dx)
        : m_dx(a_dx), m_center(a_center), m_bg_params(a_bg_params), m_tensor_mass(tensor_mass)
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

        //data_t phi_Re = m_amplitude_re;
        //data_t Pi_Im = m_amplitude_im * m_omega;

        Vars<data_t> vars;
        VarsTools::assign(vars,0.);
        Tensor<2,data_t> fspatial; //Spatial component of the tensor field
        Tensor<1,data_t> fbar; //Half-projected component of the tensor field
        data_t fhat; //Scalar part of the tensor field 
        //Conjugate components
        Tensor<2,data_t> u; // Spatial rank 2 u field
        Tensor<2,data_t> v; //Spatial rank 2 v field
        Tensor<1,data_t> p; //Spatial rank 1 p field
        Tensor<1,data_t> q; //Spatial rank 1 q field
        data_t w; //Scalar component
        const data_t initial_constant = 1e-5; 
        vars.fhat = initial_constant;
        vars.w = 0;
        for(int i = 0; i < 3; i++)
        {
          vars.fbar[i] = initial_constant;
          vars.p[i] = 0;
          vars.q[i] = 0;
          for(int j = 0; j < 3; j++)
          {
            vars.fspatial[i][j] = initial_constant;
            vars.u[i][j] = 0;
            vars.v[i][j] = 0;
          }
        }
        current_cell.store_vars(vars);
        // Store the initial values of the variables
        //current_cell.store_vars(phi_Re, c_phi_Re);
        //current_cell.store_vars(0.0, c_phi_Im);
        //current_cell.store_vars(0.0, c_Pi_Re);
        //current_cell.store_vars(Pi_Im, c_Pi_Im);
        /*
        current_cell.store_vars(fspatial, c_fspatial);
        current_cell.store_vars(fbar, c_fbar);
        current_cell.store_vars(fhat, c_fhat);
        current_cell.store_vars(u, c_u);
        current_cell.store_vars(v, c_v);
        current_cell.store_vars(p, c_p);
        current_cell.store_vars(q, c_q);
        current_cell.store_vars(0, c_w);
        */
    }
};

#endif /* INITIALCONDITIONS_HPP_ */
