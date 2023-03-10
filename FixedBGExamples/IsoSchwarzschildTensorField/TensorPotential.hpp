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
        double tensor_mass = 0.0;
    };

  private:
    const params_t m_params;

  public:
    //! The constructor
    TensorPotential(const params_t a_params) : m_params(a_params) {}

    double return_params() const
    {
        const double m = m_params.tensor_mass;
        return m;
    }
    //! Set the potential function for the tensor field here
    /*
    template <class data_t, template <typename> class vars_t>
    // using MetricVars = typename ADMFixedBGVars::template Vars<data_t>;

    void compute_potential(data_t &V_of_F, data_t &dVdF,
                           const vars_t<data_t> &vars) const
    {
        const double m = m_params.tensor_mass;

    }
    */
};

#endif /* TENSORPOTENTIAL_HPP_ */
