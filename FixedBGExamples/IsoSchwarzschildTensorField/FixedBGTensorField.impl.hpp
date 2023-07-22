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

template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> FixedBGTensorField::compute_emtensor(
    const vars_t<data_t> &vars, const MetricVars<data_t> &metric_vars,
    const vars_t<Tensor<1, data_t>> &d1, const Tensor<2, data_t> &gamma_UU,
    const Tensor<3, data_t> &chris_phys_ULL) const
{
    emtensor_t<data_t> out;

    // call the function which computes the em tensor excluding the mass terms
    emtensor_excl_potential(out, vars, metric_vars, d1, gamma_UU,
                            chris_phys_ULL);

    double mass = m_tensor_field_mass;

    // add in the mass terms
    out.rho += 0;                      // V_of_phi;
    out.S += 0;                        //-3.0 * V_of_phi;
    FOR2(i, j) { out.Sij[i][j] += 0; } //-metric_vars.gamma[i][j] * V_of_phi; }

    return out;
}

// Calculate the stress energy tensor elements
template <class data_t, template <typename> class vars_t>
void FixedBGTensorField::emtensor_excl_potential(
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
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void FixedBGTensorField::matter_rhs(rhs_vars_t<data_t> &total_rhs,
                                    const vars_t<data_t> &vars,
                                    const MetricVars<data_t> &metric_vars,
                                    const vars_t<Tensor<1, data_t>> &d1,
                                    const diff2_vars_t<Tensor<2, data_t>> &d2,
                                    const vars_t<data_t> &advec) const
{
    // call the function for the rhs excluding the potential
    matter_rhs_excl_potential(total_rhs, vars, metric_vars, d1, d2, advec);
}

// the RHS excluding the potential terms
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void FixedBGTensorField::matter_rhs_excl_potential(
    rhs_vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    using namespace TensorAlgebra;

    const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
    const auto chris_phys = compute_christoffel(metric_vars.d1_gamma, gamma_UU);

    // Evolution equations for the field and the conjugate variables:
    data_t hhat = 0.0;
    FOR2(i, j)
    {
        hhat += -metric_vars.lapse * gamma_UU[i][j] * vars.fspatial[i][j];
    }

    FOR1(i)
    {
        rhs.fbar[i] = -hhat * metric_vars.d1_lapse[i];
        FOR2(j, k)
        {
            rhs.fbar[i] += gamma_UU[j][k] * metric_vars.lapse *
                           (vars.fspatial[i][j] * metric_vars.d1_lapse[k] +
                            metric_vars.lapse * d1.fspatial[i][j][k]);

            FOR1(l)
            {
                rhs.fbar[i] += -gamma_UU[j][k] * metric_vars.lapse *
                               metric_vars.lapse *
                               (chris_phys.ULL[l][k][i] * vars.fspatial[l][j] +
                                chris_phys.ULL[l][k][j] * vars.fspatial[i][l]);
            }
        }
    }

    FOR2(i, j)
    {
        rhs.fspatial[i][j] =
            -metric_vars.lapse * vars.v[i][j] -
            vars.fbar[i] * metric_vars.d1_ln_lapse[j] +
            vars.fbar[j] * metric_vars.d1_ln_lapse[i];
    }

    // No longer really an evolution variable, just for tracking the value of
    // hhat
    rhs.fhat = 0;
    FOR2(i, j)
    {
        rhs.fhat += -metric_vars.lapse * rhs.fspatial[i][j] * gamma_UU[i][j];
    }

    // Saves having a bunch of indices
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
            rhs.v[i][j] = metric_vars.lapse * m_tensor_field_mass *
                              m_tensor_field_mass * vars.fspatial[i][j] -
                          2.0 * metric_vars.lapse * tensorRiemannTerm[i][j] +
                          2.0 * hhat * metric_vars.ricci_phys[i][j];

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
}

#endif /* FIXEDBGTENSORFIELD_IMPL_HPP_ */
