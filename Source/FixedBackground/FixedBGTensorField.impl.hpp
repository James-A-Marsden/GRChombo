/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(FIXEDBGTENSORFIELD_HPP_)
#error "This file should only be included through FixedBGTensorField.hpp"
#endif

#ifndef FIXEDBGTENSORFIELD_IMPL_HPP_
#define FIXEDBGTENSORFIELD_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> FixedBGTensorField<potential_t>::compute_emtensor(
    const vars_t<data_t> &vars, const MetricVars<data_t> &metric_vars,
    const vars_t<Tensor<1, data_t>> &d1, const Tensor<2, data_t> &gamma_UU,
    const Tensor<3, data_t> &chris_phys_ULL) const
{
    emtensor_t<data_t> out;

    // call the function which computes the em tensor excluding the potential
    emtensor_excl_potential(out, vars, metric_vars, d1, gamma_UU,
                            chris_phys_ULL);

    // set the potential values
    data_t V_of_F = 0.0;
    data_t dVdF = 0.0;
    my_potential.compute_potential(V_of_F, dVdF, vars);

    out.rho += 0;//V_of_phi;
    out.S += 0;//-3.0 * V_of_phi;
    FOR2(i, j) { out.Sij[i][j] += 0;}//-metric_vars.gamma[i][j] * V_of_phi; }

    return out;
}

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void FixedBGTensorField<potential_t>::emtensor_excl_potential(
    emtensor_t<data_t> &out, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &gamma_UU, const Tensor<3, data_t> &chris_phys_ULL)
{
    //Initially, set all the matter components to zero.
    // Useful quantity Vt
    data_t Vt = 0;//-vars.Pi * vars.Pi;
    FOR2(i, j) { Vt += 0;}//gamma_UU[i][j] * d1.phi[i] * d1.phi[j]; }

    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR2(i, j)
    {
        out.Sij[i][j] =
            0;//-0.5 * metric_vars.gamma[i][j] * Vt + d1.phi[i] * d1.phi[j];
    }

    // S = Tr_S_ij
    out.S = 0//TensorAlgebra::compute_trace(out.Sij, gamma_UU);

    // S_i (note lower index) = - n^a T_a0
    FOR1(i) { out.Si[i] = 0;}//-d1.phi[i] * vars.Pi; }

    // rho = n^a n^b T_ab
    out.rho = 0;//vars.Pi * vars.Pi + 0.5 * Vt;
}

// Adds in the RHS for the matter vars
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void FixedBGTensorField<potential_t>::matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    // call the function for the rhs excluding the potential
    matter_rhs_excl_potential(total_rhs, vars, metric_vars, d1, d2, advec);

    // set the potential values
    data_t V_of_F = 0.0;
    data_t dVdF = 0.0;

    // compute potential
    my_potential.compute_potential(V_of_F, dVdF, vars);

    // adjust RHS for the potential term
    //total_rhs.Pi += -metric_vars.lapse * dVdphi;
}

// the RHS excluding the potential terms
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void FixedBGTensorField<potential_t>::matter_rhs_excl_potential(
    rhs_vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec)
{
    using namespace TensorAlgebra;

    const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
    const auto chris_phys = compute_christoffel(metric_vars.d1_gamma, gamma_UU);

    // evolution equations for tensor field and (minus) its conjugate momentum
    //rhs.phi = metric_vars.lapse * vars.Pi + advec.phi;
   // rhs.Pi = metric_vars.lapse * metric_vars.K * vars.Pi + advec.Pi;

    rhs.fhat = advec.fhat + vars.lapse * w;

    //1 explicit index
    FOR(i)
    {   
        rhs.fbar[i] = advec.fbar[i] + 0.5 * (fhat * d1.lapse[i] + vars.lapse * d1.fbar[i])
        + 0.5 * (fbar * d1.lapse[i] - vars.lapse * p[i]);

       
       //2 explicit indices 
        FOR(j)
        {   
            rhs.fspatial[i][j] = advec.fspatial[i][j] + fbar[j] * d1.lapse[i] + fbar[i] * d1.lapse[j]
            + vars.lapse * (d1.fbar[i][j] + d1.fbar[j][i]) + 2 * vars.lapse * fhat * K_tensor[i][j] - vars.lapse * v[i][j];
            
            rhs.fhat += 2 * gamma_UU[i][j] * fbar[i] * d1.lapse[j];
           
           //3 explicit indices 
            FOR(k)
            {   
                rhs.fspatial[i][j] += 2 * vars.lapse * chris_phys.ULL[k][i][j] * fbar[k];
                
                rhs.fbar[i] += gamma_UU[j][k] * fbar[j] * d1.shift[i][k]
                + gamma_UU[j][k] * fspatial[i][j] * d1.lapse[k];
                
                //4 explicit indices
                FOR(l)
                {   
                    rhs.fspatial[i][j] += gamma_UU[k][l] * (fspatial[k][j] * d1.shift[i][l] + fspatial[i][k] * d1.shift[j][l]);              
                }
            
            }
         
        }
    
    }

}

#endif /* FIXEDBGTENSORFIELD_IMPL_HPP_ */
