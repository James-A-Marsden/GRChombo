#ifndef SETREST_HPP_
#define SETREST_HPP_

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


//! Class which creates the initial constraints
class SetRest
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
    SetRest(const double tensor_mass, const std::array<double, CH_SPACEDIM> a_center,
                      const IsoSchwarzschildFixedBG::params_t a_bg_params,
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
        IsoSchwarzschildFixedBG kerr_bh(m_bg_params, m_dx);
        MetricVars<data_t> metric_vars;
        kerr_bh.compute_metric_background(metric_vars, current_cell);
        const data_t det_gamma =
            TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
        const auto gamma_UU = TensorAlgebra::compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys = TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU);

        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const data_t x2 = x * x;
        const double y2 = y * y;
        const double z2 = z * z;
        const data_t r = coords.get_radius();
        const data_t r2 = r * r;
        const double M = m_bg_params.mass;
        const data_t rho = simd_max(sqrt(x2 + y2), 1e-6);
        const data_t rho2 = rho * rho;
 
        const data_t costheta = z/r;
        const data_t sintheta = rho/r;

        const data_t sinphi = coords.y/rho;
        const data_t cosphi = coords.x/rho;

        const data_t sin2phi = 2.0 * sinphi * cosphi; 
        const data_t cos2phi = cosphi*cosphi - sinphi*sinphi;

        Vars<data_t> vars;

        const auto local_vars = current_cell.template load_vars<Vars>();

        
        FOR2(i,j)
        {
            vars.fspatial[i][j] = local_vars.fspatial[i][j];
            vars.v[i][j] = local_vars.v[i][j];
        }
        
        
        //Calculate the derivatives
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        //Traceless condition
        
        //Traceless for derivatives 
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
        /*
        vars.Xhat = 0.5 * metric_vars.lapse * primaryScalar;

        FOR1(i)
        {
          vars.Xspatial[i] = metric_vars.lapse * primaryVector[i];
        }
        */
        
        vars.Xhat = 0.0;

        FOR1(i)
        {
            vars.Xspatial[i] = 0.0;
        }
        

        current_cell.store_vars(vars);
    }
};

#endif /* SETREST_HPP_ */

