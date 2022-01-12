/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TENSORPOTENTIAL_HPP_
#define TENSORPOTENTIAL_HPP_

#include "simd.hpp"

class TensorPotential
{
  public:
    struct params_t
    {
        double tensor_mass;
    };

  private:
    const params_t m_params;

  public:
    //! The constructor
    TensorPotential(const params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the tensor field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_F, data_t &dVdF, const vars_t<data_t> &vars) const
    {
        const double m = m_params.tensor_mass;

        // The potential value at a given value of F
        //Indices ! and metrics 

        V_of_F = vars.f_hat * vars.f_hat; 

        FOR(i)
        {

          V_of_F -= 2 * vars.f_bar[i] * vars.f_bar[i];
          
          FOR(j)
          {

            V_of_F += vars.fspatial[i][j] * vars.fspatial[i][j];

          }
        }

        V_of_F *= (0.25 * m * m);
         //          0.5 * m * m * vars.phi_Im * vars.phi_Im;

        // The potential gradient at phi wrt real and im fields
        dVdF = 0; 
       // dVdphi_Re = m * m * vars.phi_Re;
       // dVdphi_Im = m * m * vars.phi_Im;
    }
};

#endif /* TENSORPOTENTIAL_HPP_ */
