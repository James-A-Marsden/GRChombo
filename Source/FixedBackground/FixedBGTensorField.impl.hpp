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

// emtensor_t<data_t> FixedBGTensorField<potential_t>
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
     //data_t V_of_F = 0.0;
     //data_t dVdF = 0.0;
    //my_potential.compute_potential(V_of_F, dVdF, vars);

    out.rho += 0;                      // V_of_phi;
    out.S += 0;                        //-3.0 * V_of_phi;
    FOR2(i, j) { out.Sij[i][j] += 0; } //-metric_vars.gamma[i][j] * V_of_phi; }

    return out;
}

// Calculate the stress energy tensor elements

// void FixedBGTensorField<potential_t>
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void FixedBGTensorField<potential_t>::emtensor_excl_potential(
    emtensor_t<data_t> &out, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &gamma_UU, const Tensor<3, data_t> &chris_phys_ULL)
{
    // Initially, set all the matter components to zero.
    //  Useful quantity Vt
    data_t Vt = 0;          //-vars.Pi * vars.Pi;
    FOR2(i, j) { Vt += 0; } // gamma_UU[i][j] * d1.phi[i] * d1.phi[j]; }

    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR2(i, j)
    {
        out.Sij[i][j] =
            0; //-0.5 * metric_vars.gamma[i][j] * Vt + d1.phi[i] * d1.phi[j];
    }

    // S = Tr_S_ij
    out.S = 0; // TensorAlgebra::compute_trace(out.Sij, gamma_UU);

    // S_i (note lower index) = - n^a T_a0
    FOR1(i) { out.Si[i] = 0; } //-d1.phi[i] * vars.Pi; }

    // rho = n^a n^b T_ab
    out.rho = 0; // vars.Pi * vars.Pi + 0.5 * Vt;
}

// Adds in the RHS for the matter vars

// void FixedBGTensorField<potential_t>
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
    //double mass = m_tensor_field_mass;
    matter_rhs_excl_potential(total_rhs, vars, metric_vars, d1, d2, advec, m_tensor_field_mass, m_damping_kappa, m_damping_switch);
}

// the RHS excluding the potential terms

template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void FixedBGTensorField<potential_t>::matter_rhs_excl_potential(
    rhs_vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec, const double tensor_field_mass, const double damping_kappa, int damping_switch)
{
    using namespace TensorAlgebra;

    const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
    const auto chris_phys = compute_christoffel(metric_vars.d1_gamma, gamma_UU);

    
    // Evolution equations for the field and the conjugate variables:
    //replacement of fhat
    data_t fspatial_trace = 0.0;
    FOR2(i,j)
    {
        fspatial_trace += - metric_vars.lapse * gamma_UU[i][j] * vars.fspatial[i][j];
    }
    //Derivative of the trace of fspatial
    Tensor<1,data_t> d1_fspatial_trace;
    Tensor<2,data_t> d2_fspatial_trace; 

    FOR1(i)
    {
        d1_fspatial_trace[i] = 0.0;

        FOR2(j,k)
        {
            d1_fspatial_trace[i] += -metric_vars.d1_lapse[i] * gamma_UU[j][k] * vars.fspatial[j][k]
                                    -metric_vars.lapse * (metric_vars.d1_gamma_UU[j][k][i] * vars.fspatial[j][k]
                                                            +gamma_UU[j][k] * d1.fspatial[j][k][i]);
        }
    }

    FOR2(i,j)
    {
        d2_fspatial_trace[i][j] = 0.0;

        FOR2(k,l)
        {
             d2_fspatial_trace[i][j] += -metric_vars.d2_lapse[i][j] * gamma_UU[k][l] * vars.fspatial[k][l]
                                     -metric_vars.d1_lapse[i] * metric_vars.d1_gamma_UU[k][l][j] * vars.fspatial[k][l]
                                     -metric_vars.d1_lapse[i] * gamma_UU[k][l] * d1.fspatial[k][l][j]
                                     -metric_vars.d1_lapse[j] * metric_vars.d1_gamma_UU[k][l][i] * vars.fspatial[k][l]
                                     -metric_vars.lapse * metric_vars.d2_gamma_UU[k][l][i][j] * vars.fspatial[k][l]
                                     -metric_vars.lapse * metric_vars.d1_gamma_UU[k][l][i] * d1.fspatial[k][l][j]
                                     -metric_vars.d1_lapse[j] * gamma_UU[k][l] * d1.fspatial[k][l][i]
                                     -metric_vars.lapse * metric_vars.d1_gamma_UU[k][l][j] * d1.fspatial[k][l][i]
                                     -metric_vars.lapse * gamma_UU[k][l] * d2.fspatial[k][l][i][j];
        }
    }
    
    rhs.fhat = 0.0;
    FOR2(i, j)
    {
        rhs.fhat += metric_vars.lapse * metric_vars.lapse * gamma_UU[i][j] *
                        d1.fbar[i][j] +
                    2.0 * metric_vars.lapse * gamma_UU[i][j] * vars.fbar[i] *
                        metric_vars.d1_lapse[j];
        FOR1(k)
        {
            rhs.fhat += -metric_vars.lapse * metric_vars.lapse *
                        gamma_UU[i][j] * chris_phys.ULL[k][i][j] * vars.fbar[k];
        }
    }

    FOR1(i)
    {

        rhs.fbar[i] = fspatial_trace * metric_vars.d1_ln_lapse[i];
        FOR2(j, k)
        {
            rhs.fbar[i] += gamma_UU[j][k] *
                           (-vars.fspatial[i][j] * metric_vars.d1_lapse[k] -
                            metric_vars.lapse * d1.fspatial[i][j][k]);

            FOR1(l)
            {
                rhs.fbar[i] += -metric_vars.lapse * gamma_UU[j][k] *
                               (-chris_phys.ULL[l][k][i] * vars.fspatial[l][j] -
                                chris_phys.ULL[l][k][j] * vars.fspatial[i][l]);
            }
        }
    }

    FOR2(i, j)
    {
        rhs.fspatial[i][j] = -metric_vars.lapse * vars.v[i][j] -
                             vars.fbar[i] * metric_vars.d1_lapse[j] -
                             vars.fbar[j] * metric_vars.d1_lapse[i];
    }

    
    //Saves having a bunch of indices
    Tensor<2, data_t> tensorRiemannTerm;

    FOR2(i, j)
    {
        tensorRiemannTerm[i][j] = 0.0;
        FOR3(k, l, m)
        {
            FOR2(n, o)
            {
                tensorRiemannTerm[i][j] +=
                    gamma_UU[k][m] * gamma_UU[l][n] * metric_vars.gamma[i][o] *
                    vars.fspatial[m][n] *
                    metric_vars.riemann_phys_ULLL[o][k][j][l];
            }
        }
    }    
    FOR1(i)
    {
        FOR1(j)
        {

            rhs.v[i][j] = metric_vars.lapse * tensor_field_mass * tensor_field_mass *
                              vars.fspatial[i][j] -
                          2.0 * metric_vars.lapse * tensorRiemannTerm[i][j] +
                          2.0 * fspatial_trace * metric_vars.ricci_phys[i][j];
            FOR1(k)
            {
                FOR1(l)
                {
                    rhs.v[i][j] +=
                        -gamma_UU[k][l] *
                            (d1.fspatial[i][j][k] * metric_vars.d1_lapse[l] +
                             metric_vars.lapse * d2.fspatial[i][j][k][l]) -
                        gamma_UU[k][l] *
                            (d1.fspatial[i][k][l] * metric_vars.d1_lapse[j] +
                             d1.fspatial[j][k][l] * metric_vars.d1_lapse[i]);
                    FOR1(m)
                    {
                        rhs.v[i][j] +=
                            gamma_UU[k][l] *
                            (chris_phys.ULL[m][k][i] * vars.fspatial[m][j] +
                             chris_phys.ULL[m][k][j] * vars.fspatial[i][m]) *
                            metric_vars.d1_lapse[l];

                        rhs.v[i][j] +=
                            -gamma_UU[k][l] * metric_vars.lapse *
                            (-chris_phys.ULL[m][l][k] * d1.fspatial[i][j][m] -
                             chris_phys.ULL[m][l][i] * d1.fspatial[m][j][k] -
                             chris_phys.ULL[m][l][j] * d1.fspatial[i][m][k] -
                             metric_vars.d1_chris_phys[m][k][i][l] *
                                 vars.fspatial[m][j] -
                             chris_phys.ULL[m][k][i] * d1.fspatial[m][j][l] -
                             metric_vars.d1_chris_phys[m][j][k][l] *
                                 vars.fspatial[i][m] -
                             chris_phys.ULL[m][j][k] * d1.fspatial[i][m][l]);

                        rhs.v[i][j] +=
                            gamma_UU[k][l] *
                            ((chris_phys.ULL[m][l][i] * vars.fspatial[m][k] +
                              chris_phys.ULL[m][l][k] * vars.fspatial[i][m]) *
                                 metric_vars.d1_lapse[j] +
                             (chris_phys.ULL[m][l][j] * vars.fspatial[m][k] +
                              chris_phys.ULL[m][l][k] * vars.fspatial[j][m]) *
                                 metric_vars.d1_lapse[i]);
                        FOR1(n)
                        {
                            rhs.v[i][j] += -gamma_UU[k][l] * metric_vars.lapse *
                                           (chris_phys.ULL[n][l][k] *
                                                chris_phys.ULL[m][n][i] *
                                                vars.fspatial[m][j] +
                                            chris_phys.ULL[n][l][i] *
                                                chris_phys.ULL[m][k][n] *
                                                vars.fspatial[m][j] +
                                            chris_phys.ULL[n][l][j] *
                                                chris_phys.ULL[m][k][i] *
                                                vars.fspatial[m][n] +
                                            chris_phys.ULL[n][l][k] *
                                                chris_phys.ULL[m][n][j] *
                                                vars.fspatial[i][m] +
                                            chris_phys.ULL[n][l][j] *
                                                chris_phys.ULL[m][k][n] *
                                                vars.fspatial[i][m] +
                                            chris_phys.ULL[n][l][i] *
                                                chris_phys.ULL[m][k][j] *
                                                vars.fspatial[n][m]);
                        }
                    }
                }
            }
        }
    }

    //We need constraints here for the damping variables
    data_t primaryScalar =
        metric_vars.lapse * tensor_field_mass * tensor_field_mass * fspatial_trace;

    FOR2(i, j)
    {
        primaryScalar +=
            metric_vars.lapse * gamma_UU[i][j] *
                (fspatial_trace * metric_vars.ricci_phys[i][j] +
                 2.0 * d1_fspatial_trace[i] * metric_vars.d1_ln_lapse[j] -
                 2.0 * fspatial_trace * metric_vars.d1_ln_lapse[i] *
                     metric_vars.d1_ln_lapse[j] -
                 d2_fspatial_trace[i][j]) +
            gamma_UU[i][j] * fspatial_trace * metric_vars.d2_lapse[i][j];

        FOR1(k)
        {
            primaryScalar += metric_vars.lapse * gamma_UU[i][j] *
                             (chris_phys.ULL[k][i][j] * d1_fspatial_trace[k] -
                              fspatial_trace * chris_phys.ULL[k][i][j] *
                                  metric_vars.d1_ln_lapse[k]);

            FOR1(l)
            {
                primaryScalar += gamma_UU[i][k] * gamma_UU[j][l] *
                                 metric_vars.lapse * metric_vars.lapse *
                                 metric_vars.ricci_phys[k][l] *
                                 vars.fspatial[i][j];
                primaryScalar += gamma_UU[i][k] * gamma_UU[j][l] *
                                 metric_vars.lapse *
                                 (-metric_vars.lapse * d2.fspatial[i][j][k][l]);

                FOR1(m)
                {
                    primaryScalar +=
                        -gamma_UU[i][k] * gamma_UU[j][l] * metric_vars.lapse *
                        metric_vars.lapse *
                        (-chris_phys.ULL[m][l][k] * d1.fspatial[i][j][m] -
                         chris_phys.ULL[m][l][i] * d1.fspatial[m][j][k] -
                         chris_phys.ULL[m][l][j] * d1.fspatial[i][m][k] -
                         metric_vars.d1_chris_phys[m][k][i][l] *
                             vars.fspatial[m][j] -
                         chris_phys.ULL[m][k][i] * d1.fspatial[m][j][l] -
                         metric_vars.d1_chris_phys[m][j][k][l] *
                             vars.fspatial[i][m] -
                         chris_phys.ULL[m][j][k] * d1.fspatial[i][m][l]);

                    FOR1(n)
                    {
                        primaryScalar +=
                            -gamma_UU[i][k] * gamma_UU[j][l] *
                            metric_vars.lapse * metric_vars.lapse *
                            (chris_phys.ULL[n][l][k] * chris_phys.ULL[m][n][i] *
                                 vars.fspatial[m][j] +
                             chris_phys.ULL[n][l][i] * chris_phys.ULL[m][k][n] *
                                 vars.fspatial[m][j] +
                             chris_phys.ULL[n][l][j] * chris_phys.ULL[m][k][i] *
                                 vars.fspatial[m][n] +
                             chris_phys.ULL[n][l][k] * chris_phys.ULL[m][n][j] *
                                 vars.fspatial[i][m] +
                             chris_phys.ULL[n][l][j] * chris_phys.ULL[m][k][n] *
                                 vars.fspatial[i][m] +
                             chris_phys.ULL[n][l][i] * chris_phys.ULL[m][k][j] *
                                 vars.fspatial[n][m]);
                    }
                }
            }
        }
    }

    Tensor<1, data_t> primaryVector;
    FOR1(i)
    {
        primaryVector[i] =
            metric_vars.lapse * tensor_field_mass * tensor_field_mass * vars.fbar[i];

        FOR2(j, k)
        {
            primaryVector[i] +=
                metric_vars.lapse * gamma_UU[j][k] *
                (d1.v[j][i][k] + vars.fbar[i] * metric_vars.ricci_phys[j][k] -
                 vars.fbar[j] * metric_vars.ricci_phys[i][k] -
                 d2.fbar[i][j][k]);

            FOR1(l)
            {
                primaryVector[i] +=
                    metric_vars.lapse * gamma_UU[j][k] *
                    (-chris_phys.ULL[l][k][i] * vars.v[j][l] -
                     chris_phys.ULL[l][k][j] * vars.v[l][i] +
                     chris_phys.ULL[l][k][j] * d1.fbar[i][l] +
                     chris_phys.ULL[l][k][i] * d1.fbar[l][j] +
                     metric_vars.d1_chris_phys[l][j][i][k] * vars.fbar[l] +
                     d1.fbar[l][k] * chris_phys.ULL[l][j][i]);

                FOR1(m)
                {
                    primaryVector[i] +=
                        metric_vars.lapse * gamma_UU[j][k] *
                        (-chris_phys.ULL[m][k][j] * chris_phys.ULL[l][m][i] *
                             vars.fbar[l] -
                         chris_phys.ULL[m][k][i] * chris_phys.ULL[l][j][m] *
                             vars.fbar[l]);
                }
            }
        }
    }

    // Damping evolution


    
    //rhs.Xhat = - 0.5 * damping_kappa * vars.Xhat - 0.5 * primaryScalar /
    //metric_vars.lapse;
    rhs.Xhat = - 0.5 * damping_kappa * vars.Xhat + 0.5 * primaryScalar /
    metric_vars.lapse;    

    FOR1(i)
    {
        rhs.Xspatial[i] =  - damping_kappa * vars.Xspatial[i] +
    //metric_vars.lapse * d1.Xhat[i] - vars.Xhat * metric_vars.d1_lapse[i] + primaryVector[i];
    metric_vars.lapse * d1.Xhat[i] - vars.Xhat * metric_vars.d1_lapse[i] - primaryVector[i];

        FOR1(j)
        {

            rhs.Xhat += gamma_UU[i][j] * ( - vars.Xspatial[i] * metric_vars.d1_lapse[j]);

        }
    }
    
    if(damping_switch > 0.5)
    FOR2(i,j)
    {
        rhs.v[i][j] -=  metric_vars.gamma[i][j] * damping_kappa * vars.Xhat
            - metric_vars.lapse * (d1.Xspatial[i][j] + d1.Xspatial[j][i]); 
    
        FOR1(k)
        {
                rhs.v[i][j] -= metric_vars.lapse * 2.0 * chris_phys.ULL[k][i][j] * vars.Xspatial[k];
        }
    }
        
}

#endif /* FIXEDBGTENSORFIELD_IMPL_HPP_ */
