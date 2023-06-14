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
#include "FourthOrderDerivatives.hpp"
#include "IsoSchwarzschildFixedBG.hpp"
#include "Tensor.hpp"
#include "TensorPotential.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include <boost/math/special_functions/bessel.hpp>
//! Class which creates the initial conditions
class InitialConditions
{
  protected:
    const double m_dx;
    // const double m_amplitude_re, m_amplitude_im;
    // const double m_omega;
    const std::array<double, CH_SPACEDIM> m_center;
    const IsoSchwarzschildFixedBG::params_t m_bg_params;
    const FourthOrderDerivatives m_deriv;

    // load in Vars from the field
    //  The evolution vars

    template <class data_t>
    using Vars = FixedBGTensorField<TensorPotential>::template Vars<data_t>;
    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  public:
    //! The constructor for the class
    // const double a_amplitude_re, const double a_amplitude_im, const double
    // a_omega,
    InitialConditions(const std::array<double, CH_SPACEDIM> a_center,
                      const IsoSchwarzschildFixedBG::params_t a_bg_params,
                      const double a_dx)
        : m_dx(a_dx), m_center(a_center), m_bg_params(a_bg_params),
          m_deriv(a_dx)
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
        const auto gamma_UU =
            TensorAlgebra::compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys =
            TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU);

        Vars<data_t> vars;
        VarsTools::assign(vars, 0.);
        Tensor<2, data_t> fspatial; // Spatial component of the tensor field
        Tensor<1, data_t> fbar;
        data_t fhat;
        Tensor<2, data_t> v; // Spatial rank 2 v field

        // Damping fields
        Tensor<1, data_t> Xspatial;
        data_t Xhat;

        // Set everything to zero
        vars.fhat = 0.0;
        vars.Xhat = 0.0;
        FOR1(i)
        {
            vars.fbar[i] = 0.0;
            vars.Xspatial[i] = 0.0;
        }

        FOR2(i, j)
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
        const data_t rho = coords.get_rho(); // simd_max(sqrt(x2 + y2), 1e-6);
        const data_t rho2 = rho * rho;

        const data_t costheta = z / r;
        const data_t sintheta = rho / r;

        const data_t cos2theta = costheta * costheta - sintheta * sintheta;
        const data_t sin2theta = 2.0 * sintheta * costheta;

        const data_t sinphi = y / rho;
        const data_t cosphi = x / rho;

        const data_t sin2phi = 2.0 * sinphi * cosphi;
        const data_t cos2phi = cosphi * cosphi - sinphi * sinphi;

        const double M = m_bg_params.mass;
        const double M2 = M * M;

        data_t A = pow(M + 2.0 * r, -2.0);
        /*
        vars.fspatial[0][0] = -A * sintheta * sintheta * sin2phi / r;
        vars.fspatial[1][1] = A * sintheta * sintheta * sin2phi / r;
        vars.fspatial[2][2] = 0.0;
        vars.fspatial[0][1] = A * sintheta * sintheta * cos2phi / r;
        vars.fspatial[0][2] = -A * sintheta * sinphi * costheta / r;
        vars.fspatial[1][2] = A * sintheta * cosphi * costheta / r;

        vars.fspatial[1][0] = vars.fspatial[0][1];
        vars.fspatial[2][0] = vars.fspatial[0][2];
        vars.fspatial[2][1] = vars.fspatial[1][2];
        */
        
        vars.v[0][0] = -2 * A * cosphi * sinphi * pow(sintheta, 2.0) / r;
        vars.v[1][1] = 2 * A * cosphi * sinphi * pow(sintheta, 2.0) / r;
        vars.v[2][2] = 0;
        vars.v[0][1] = A * cos2phi * pow(sintheta, 2.0) / r;
        vars.v[0][2] = -A * costheta * sintheta * sinphi / r;
        vars.v[1][2] = A * costheta * sintheta * cosphi / r;

        vars.v[1][0] = vars.v[0][1];
        vars.v[2][0] = vars.v[0][2];
        vars.v[2][1] = vars.v[1][2];
      
        current_cell.store_vars(vars);

        // vvv nothing to see below this line

        /*
        Tensor<2,data_t> sph_to_cart_J_UD;
        sph_to_cart_J_UD[0][0] = x/r;
        sph_to_cart_J_UD[0][1] = y/r;
        sph_to_cart_J_UD[0][2] = z/r;
        sph_to_cart_J_UD[1][0] = x * z /(r2 * rho);
        sph_to_cart_J_UD[1][1] = y * z /(r2 * rho);
        sph_to_cart_J_UD[1][2] = -rho / r2;
        sph_to_cart_J_UD[2][0] = -y/rho2;
        sph_to_cart_J_UD[2][1] = x/rho2;
        sph_to_cart_J_UD[2][2] = 0.0;
        const double lam = 1.0;
        const double lam2 = lam * lam;
        const double beta = 20.0;
        //Teukolsky Wave Solutions

        //Teukolsky 'domain wall'
        const data_t tanh1 = tanh((r - beta)/lam);
        const data_t tanh2 = tanh1 * tanh1;
        const data_t sech2 = (1.0/cosh((r - beta)/lam)) * (1.0 / cosh((r -
        beta)/lam));

        const data_t seed = (1e-5) * r5 * (1.0 + tanh1);
        const data_t d1seed =  (1e-5) * -r4 * (r * sech2 + 5.0 * lam * (1.0 +
        tanh1)) / lam; const data_t d2seed = (1e-5) * 2.0 * r3 * (10.0 * (1.0 +
        tanh1) + r * sech2 * (5.0 * lam - r * tanh1) / lam2); const data_t
        d3seed = (1e-5) * pow(lam,-3.0) * 2.0 * r2 * (r3 * sech2 * sech2 - 30.0
        * lam2 * lam * (1.0 + tanh1) + r * sech2 * (-30.0 * lam2 + 15.0 * r *
        lam * tanh1 - 2.0 * r2 * tanh1 * tanh1)); const data_t d4seed = (1e-5) *
        pow(lam, -4.0) * 8.0 * r * (15.0 * lam2 * lam2 * (1.0 + tanh1) + r3 *
        sech2 * sech2 * (-5.0 * lam + 2.0 * r * tanh1) + r * sech2 * (30.0 *
        lam2 * lam - 30.0 * r * lam2 * tanh1 + 10.0 * r2 * lam * tanh1 * tanh1 -
        r3 * tanh1 * tanh1 * tanh1));

        const data_t seed_momentum = -(1e-5) * (r5 * sech2 / lam + 5.0 * r4 *
        (1.0 + tanh1)); const data_t d1seed_momentum = -(1e-5) * 2.0 * r3 *
        (-10.0 * (1.0 + tanh1) + r * sech2 * (-5.0 * lam + r * tanh1) / lam2);
        const data_t d2seed_momentum = -(1e-5) * pow(lam,-3.0) * 2.0 * r2 * (-r3
        * sech2 * sech2 + 30.0 * lam2 * lam * (1.0 + tanh1) + r * sech2 * (30.0
        * lam2 - 15.0 * r * lam * tanh1 + 2.0 * r2 * tanh2)); const data_t
        d3seed_momentum = -(1e-5) * -pow(lam,-4.0) * 8.0 * r * (15.0 * lam2 *
        lam2 * (1.0 + tanh1) + r3 * sech2 * sech2 * (-5.0 * lam + 2.0 * r *
        tanh1) + r * sech2 * (30.0 * lam2 * lam - 30.0 * r * lam2 * tanh1 + 10.0
        * r2 * lam * tanh2 - r3 * tanh2 * tanh1)); const data_t d4seed_momentum
        = -(1e-5) * pow(lam,-5.0) * 8.0 * (2.0 * r5 * pow(sech2,3.0) + 15.0 *
        pow(lam,5.0) * (1.0 + tanh1) + r3 * sech2 * sech2 * (-50.0 * lam2 + 50.0
        * r * lam * tanh1 - 11.0 * r2 * tanh2) + r * sech2 * (75.0 * lam2 * lam2
        - 150.0 * r * lam2 * lam * tanh1 + 100.0 * r2 * lam2 * tanh2 - 25.0 * r3
        * lam * tanh2 * tanh1 + 2.0 * r4 * tanh2 * tanh2));
        */
        /*
        const data_t seed = exp(-pow(r-beta,2.0) * pow(lam,-2.0)/2.0);
        const data_t d1seed = exp(-pow(r-beta,2.0) * pow(lam,-2.0)/2.0) * (r -
        beta) * pow(lam,-2.0); const data_t d2seed = exp(-pow(r-beta,2.0) *
        pow(lam,-2.0)/2.0) * (-lam * lam + pow(beta - r,2.0)) * pow(lam,-4.0);
        const data_t d3seed = exp(-pow(r-beta,2.0) * pow(lam,-2.0)/2.0) * (3.0 *
        lam * lam -  pow(beta-r,2.0)) * (beta-r) * pow(lam,-6.0); const data_t
        d4seed = exp(-pow(r-beta,2.0) * pow(lam,-2.0)/2.0) * (3.0 * pow(lam,4.0)
        - 6.0 * lam * lam * pow(beta-r,2.0) + pow(beta-r,4.0)) * pow(lam,-8.0);

        const data_t seed_momentum = exp(-pow(r-beta,2.0) * pow(lam,-2.0)/2.0) *
        (r - beta); const data_t d1seed_momentum = exp(-pow(r-beta,2.0) *
        pow(lam,-2.0)/2.0) * (-lam * lam * pow(beta-r,2.0)) * pow(lam,-2.0);
        const data_t d2seed_momentum = exp(-pow(r-beta,2.0) * pow(lam,-2.0)/2.0)
        * -1.0 * (-3.0 * lam * lam + pow(beta - r,2.0)) * (beta - r) *
        pow(lam,-4.0); const data_t d3seed_momentum = exp(-pow(r-beta,2.0) *
        pow(lam,-2.0)/2.0) * (3.0 * pow(lam,4.0) - 6.0 * lam * lam * pow(beta -
        r,2.0) + pow(beta-r,4.0)) * pow(lam,-6.0); const data_t d4seed_momentum
        = exp(-pow(r-beta,2.0) * pow(lam,-2.0)/2.0) * -1.0 * (15.0 *
        pow(lam,4.0) - 10.0 * lam * lam * pow(beta-r,2.0) + pow(beta-r,4.0)) *
        (beta-r) * pow(lam,-8.0);
        */

        /*
        const data_t seed = exp(-pow(r-beta,2.0)) * r5;
        const data_t d1seed = exp(-pow(r-beta,2.0)) * r4 * (-5.0 * lam + 2.0 * r
        * (r - beta) ) / lam; const data_t d2seed = exp(-pow(r-beta,2.0)) * 2.0
        * r3 * (10.0 * lam * lam + lam *(10.0 * beta - 11.0 * r) * r + 2.0 *
        pow(beta-r,2.0) * r2) * pow(lam,-2.0); const data_t d3seed =
        exp(-pow(r-beta,2.0)) * -2.0 * r2 * (30.0 * pow(lam,3.0) + 15.0 * lam *
        lam * (4.0 * beta - 5.0 * r) * r + 4.0 * pow(beta-r,3.0) * r3 + 6.0 *
        lam * r2 * (5.0 * beta * beta - 11.0 * beta * r + 6.0 * r2)) *
        pow(lam,-3.0); const data_t d4seed = exp(-pow(r-beta,2.0)) * 4.0 * r *
        (30.0 * pow(lam,4.0) + 60.0 * pow(lam,3.0) * (2.0 * beta - 3.0 * r) * r
        + 4.0 * lam * (10.0 * beta - 13.0 * r) * pow(beta - r,2.0) * r3 + 4.0 *
        pow(beta-r,4.0) * r4 + 3.0 * lam * lam * r2 * (40.0 * beta * beta -
        100.0 * beta * r + 61.0 * r2)) * pow(lam,-4.0);
        */
        /*
        const data_t seed = exp(-pow(r-1.0,2.0)/lam) * (r-1.0);
        const data_t d1seed = exp(-pow(r-1.0,2.0)/lam) * (-lam + 2.0 * pow(r
        - 1.0,2.0) )/ lam; const data_t d2seed = exp(-pow(r-1.0,2.0)/lam) * 2.0
        * (-3.0 * lam + 2.0 * pow(r - 1.0,2.0)) * (r-1.0) * pow(lam,-2.0); const
        data_t d3seed = exp(-pow(r-1.0,2.0)/lam) * (6.0 * lam * lam - 24.0 * lam
        * pow(r-1.0,2.0) + 8.0 * pow(r -1.0,4.0)) * pow(lam,-3.0); const data_t
        d4seed = exp(-pow(r-1.0,2.0)/lam) * 4.0 * (15.0 * lam * lam - 20.0 * lam
        * pow(r-1.0,2.0) + 4.0 * pow(r-1.0,4.0)) * (r-1.0) * pow(lam,-4.0);
        */

        /*
        const data_t A_func = 3.0 * (d2seed/r3 + 3.0 * d1seed/r4 + 3.0 *
        seed/r5); const data_t B_func = - 1.0 * (d3seed/r2 + 3.0 * d2seed/r3
        + 6.0 * d1seed/r4 + 6.0 * seed/r5); const data_t C_func = 0.25 *
        (d4seed/r + 2.0 * d3seed/r2 + 9.0 * d2seed/r3 + 21.0 * d1seed/r4 + 21.0
        * seed/r5);

        const data_t A_func_momentum = 3.0 * (d2seed_momentum/r3 + 3.0 *
        d1seed_momentum/r4 + 3.0 * seed_momentum/r5); const data_t
        B_func_momentum = - 1.0 * (d3seed_momentum/r2 + 3.0 * d2seed_momentum/r3
        + 6.0 * d1seed_momentum/r4 + 6.0 * seed_momentum/r5); const data_t
        C_func_momentum = 0.25 * (d4seed_momentum/r + 2.0 * d3seed_momentum/r2
        + 9.0 * d2seed_momentum/r3 + 21.0 * d1seed_momentum/r4 + 21.0 *
        seed_momentum/r5);



        const data_t frr = 2.0 - 3.0 * sintheta * sintheta;
        const data_t frt = -3.0 * costheta * sintheta;
        const data_t frp = 0.0;
        const data_t ftt_1 = 3.0 * sintheta * sintheta;
        const data_t ftt_2 = -1.0;
        const data_t ftp = 0.0;
        const data_t fpp_1 = -3.0 * sintheta * sintheta;
        const data_t fpp_2 = 3.0 * sintheta * sintheta - 1.0;

        Tensor<2,data_t> spherical_soln;
        spherical_soln[0][0] = A_func * frr;
        spherical_soln[0][1] = r * B_func * frt;
        spherical_soln[1][1] = r2 * (C_func * ftt_1 + A_func * ftt_2);
        spherical_soln[1][2] = 0.0;
        spherical_soln[2][2] = r2 * sintheta * sintheta * (C_func * fpp_1 +
        A_func * fpp_2); spherical_soln[0][2] = 0.0;

        spherical_soln[1][0] = spherical_soln[0][1];
        spherical_soln[2][1] = spherical_soln[1][2];
        spherical_soln[2][0] = spherical_soln[0][2];

        Tensor<2,data_t> spherical_momentum;
        spherical_momentum[0][0] = A_func_momentum * frr;
        spherical_momentum[0][1] = r * B_func_momentum * frt;
        spherical_momentum[1][1] = r2 * (C_func_momentum * ftt_1 +
        A_func_momentum * ftt_2); spherical_momentum[1][2] = 0.0;
        spherical_momentum[2][2] = r2 * sintheta * sintheta * (C_func_momentum *
        fpp_1 + A_func_momentum * fpp_2); spherical_momentum[0][2] = 0.0;

        spherical_momentum[1][0] = spherical_momentum[0][1];
        spherical_momentum[2][1] = spherical_momentum[1][2];
        spherical_momentum[2][0] = spherical_momentum[0][2];

        FOR2(i,j)
        {
          FOR2(k,l)
          {
            vars.fspatial[i][j] += sph_to_cart_J_UD[k][i] *
        sph_to_cart_J_UD[l][j] * spherical_soln[k][l]; vars.v[i][j] +=
        sph_to_cart_J_UD[k][i] * sph_to_cart_J_UD[l][j] *
        spherical_momentum[k][l];
          }
        }
        */

        /*
        namespace bmath = boost::math;

        //const double frequency = 2.0 * M_PI /128.0 ;

        const double omega = 2.0 * M_PI /4.0 ;
        //24!
        const double shift = 20.0;
        //Plane wave initial conditions


        const data_t gaussian = exp(-(coords.z - shift) * (coords.z - shift) /
        (2.0));

        const data_t field_init = sin(omega * coords.z) * gaussian;

        const data_t momentum_init = omega * cos(omega * coords.z) * gaussian
                                    -sin(omega * coords.z) * (coords.z - shift)
        * gaussian;
        */
        /*
        //'Domain Wall' initial conditions
        const data_t wall =       1.0 - tanh(shift - z) * tanh(shift + z);
        const data_t wall_deriv = tanh(shift + z) * pow(cosh(shift - z),-2.0) +
        tanh(shift - z) * pow(cosh(shift + z),-2.0); const data_t field_init =
        wall; const data_t momentum_init = wall_deriv;
        */
        // Wall 2
        /*
        const double shift = 20.0;
        const data_t wall = 0.5 * (tanh(z-shift) + tanh(z+shift));
        const data_t wall_deriv = 0.5 * (pow(cosh(z-shift),-2.0) -
        pow(cosh(z+shift),-2.0)); const data_t field_init = wall; const data_t
        momentum_init = wall_deriv;
        */
        // Wall 3 (single wall)
        /*
        const double shift = 20.0;
        const data_t wall = 1.0 + tanh(z-shift);
        const data_t wall_deriv = -pow(cosh(z-shift),-2.0);
        const data_t field_init = wall;
        const data_t momentum_init = -wall_deriv;
        */
        // PROFILE INITIAL CONDITIONS

        // data_t A = 0.1 * sintheta * sintheta * pow(M + 2.0 * r, -2.0);
        /*
        if (r > M/2.0)
        {
          //data_t A =  sintheta * sintheta * pow(M*M - 4.0 * r*r, -1.0);
          data_t A = 1e10 * sintheta * sintheta * pow(M + 2.0 * r, -2.0);
          vars.fspatial[0][0] = -A * sin2phi / r;
          vars.fspatial[1][1] =  A * sin2phi / r;
          vars.fspatial[2][2] =  0.0;
          vars.fspatial[0][1] =  A * cos2phi / r;
          vars.fspatial[0][2] = -A * sinphi * costheta / (sintheta * r);
          vars.fspatial[1][2] =  A * cosphi * costheta / (sintheta * r);

          vars.fspatial[1][0] = vars.fspatial[0][1];
          vars.fspatial[2][0] = vars.fspatial[0][2];
          vars.fspatial[2][1] = vars.fspatial[1][2];

        }
        else
        {
          FOR2(i,j)
          {
            vars.fspatial[i][j] = 0.0;
          }

        */
        // data_t A =  sintheta * sintheta * pow(M*M - 4.0 * r*r, -1.0);
        /*
        data_t A = sintheta * sintheta * pow(M + 2.0 * r, -2.0);
        vars.fspatial[0][0] = -A * sin2phi / r;
        vars.fspatial[1][1] =  A * sin2phi / r;
        vars.fspatial[2][2] =  0.0;
        vars.fspatial[0][1] =  A * cos2phi / r;
        vars.fspatial[0][2] = -A * sinphi * costheta / (sintheta * r);
        vars.fspatial[1][2] =  A * cosphi * costheta / (sintheta * r);

        vars.fspatial[1][0] = vars.fspatial[0][1];
        vars.fspatial[2][0] = vars.fspatial[0][2];
        vars.fspatial[2][1] = vars.fspatial[1][2];
        */

        /*
        data_t A =  1.0;
        vars.fspatial[0][0] = -A * sintheta * sintheta * sin2phi;
        vars.fspatial[1][1] =  A * sintheta * sintheta * sin2phi;
        vars.fspatial[2][2] =  0.0;
        vars.fspatial[0][1] =  A * sintheta * sintheta * cos2phi;
        vars.fspatial[0][2] = -A * sintheta * sinphi * costheta;
        vars.fspatial[1][2] =  A * sintheta * cosphi * costheta;

        vars.fspatial[1][0] = vars.fspatial[0][1];
        vars.fspatial[2][0] = vars.fspatial[0][2];
        vars.fspatial[2][1] = vars.fspatial[1][2];
        */
        // GW INITIAL CONDITIONS
        /*
        vars.fspatial[0][0] = field_init;
        vars.fspatial[1][1] = -field_init;

        vars.v[0][0] = -momentum_init;
        vars.v[1][1] = momentum_init;
        */
    }
};

#endif /* INITIALCONDITIONS_HPP_ */
