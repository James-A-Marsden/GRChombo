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
#include "IsoSchwarzschildFixedBG.hpp"
#include "Tensor.hpp"
#include "TensorPotential.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "FourthOrderDerivatives.hpp"
 #include <boost/math/special_functions/bessel.hpp> 
//! Class which creates the initial conditions
class InitialConditions
{
  protected:
    const double m_dx;
    //const double m_amplitude_re, m_amplitude_im;
    //const double m_omega;
    const std::array<double, CH_SPACEDIM> m_center;
    const IsoSchwarzschildFixedBG::params_t m_bg_params;
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
    InitialConditions(const double tensor_mass, const std::array<double, CH_SPACEDIM> a_center,
                      const IsoSchwarzschildFixedBG::params_t a_bg_params,
                      const double a_dx, const double a_initial_constant)
        : m_dx(a_dx), m_center(a_center), m_bg_params(a_bg_params), m_tensor_mass(tensor_mass),
        m_initial_constant(a_initial_constant), m_deriv(a_dx)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // get the metric vars
        IsoSchwarzschildFixedBG kerr_bh(m_bg_params, m_dx);
        MetricVars<data_t> metric_vars;
        kerr_bh.compute_metric_background(metric_vars, current_cell);
        const data_t det_gamma =
            TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
        const auto gamma_UU = TensorAlgebra::compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys = TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU);

 
        Vars<data_t> vars;
        VarsTools::assign(vars,0.);
        Tensor<2,data_t> fspatial; //Spatial component of the tensor field
        Tensor<1, data_t> fbar;
        data_t fhat;
        Tensor<2,data_t> v; //Spatial rank 2 v field

        //Damping fields
        Tensor<1,data_t> Xspatial;
        data_t Xhat;


        
       
       
        //Set everything to zero
        vars.fhat = 0.0;
        vars.Xhat = 0.0;
        FOR1(i)
        {
          vars.fbar[i] = 0.0;
          vars.Xspatial[i] = 0.0;
        }

        FOR2(i,j)
        {
          vars.fspatial[i][j] = 0.0;
          vars.v[i][j] = 0.0;
        } 
        double radius = coords.get_radius();
        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const data_t x2 = x * x;
        const double y2 = y * y;
        const double z2 = z * z;
        const data_t r = coords.get_radius();
        const data_t r2 = r * r;
        const data_t r3 = r2 * r;
        const data_t r4 = r2 * r2;
        const data_t r5 = r2 * r3;
        const data_t rho = coords.get_rho();//simd_max(sqrt(x2 + y2), 1e-6);
        const data_t rho2 = rho * rho;
 
        const data_t costheta = z/r;
        const data_t sintheta = rho/r;

        const data_t cos2theta = costheta*costheta - sintheta*sintheta;
        const data_t sin2theta = 2.0 * sintheta * costheta;

        const data_t sinphi = y/rho;
        const data_t cosphi = x/rho;

        const data_t sin2phi = 2.0 * sinphi * cosphi; 
        const data_t cos2phi = cosphi*cosphi - sinphi*sinphi;

        const double M = m_bg_params.mass;
        const double M2 = M * M;

        
          
        data_t A = pow(M + 2.0 * r, -2.0);
        vars.fspatial[0][0] = -A * sintheta * sintheta * sin2phi / r;
        vars.fspatial[1][1] =  A * sintheta * sintheta * sin2phi / r;
        vars.fspatial[2][2] =  0.0;
        vars.fspatial[0][1] =  A * sintheta * sintheta * cos2phi / r;
        vars.fspatial[0][2] = -A * sintheta * sinphi * costheta / r;
        vars.fspatial[1][2] =  A * sintheta * cosphi * costheta / r;

        vars.fspatial[1][0] = vars.fspatial[0][1];
        vars.fspatial[2][0] = vars.fspatial[0][2];
        vars.fspatial[2][1] = vars.fspatial[1][2];  
        
        current_cell.store_vars(vars);

    }
};

#endif /* INITIALCONDITIONS_HPP_ */
