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

//emtensor_t<data_t> FixedBGTensorField<potential_t>
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

    out.rho += 0;//V_of_phi;
    out.S += 0;//-3.0 * V_of_phi;
    FOR2(i, j) { out.Sij[i][j] += 0;}//-metric_vars.gamma[i][j] * V_of_phi; }

    return out;
}

// Calculate the stress energy tensor elements

//void FixedBGTensorField<potential_t>
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
        out.Sij[i][j] = 0;//-0.5 * metric_vars.gamma[i][j] * Vt + d1.phi[i] * d1.phi[j];
    }

    // S = Tr_S_ij
    out.S = 0;//TensorAlgebra::compute_trace(out.Sij, gamma_UU);

    // S_i (note lower index) = - n^a T_a0
    FOR1(i) { out.Si[i] = 0;}//-d1.phi[i] * vars.Pi; }

    // rho = n^a n^b T_ab
    out.rho = 0;//vars.Pi * vars.Pi + 0.5 * Vt;
}

// Adds in the RHS for the matter vars

//void FixedBGTensorField<potential_t>
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
    //data_t V_of_F = 0.0;
    //data_t dVdF = 0.0;

    // compute potential
    //my_potential.compute_potential(V_of_F, dVdF, vars);

    // adjust RHS for the potential term
    //total_rhs.Pi += -metric_vars.lapse * dVdphi;
}

// the RHS excluding the potential terms

//void FixedBGTensorField<potential_t>
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

    const double temp_mass = 0.00;

    //Defining the F variable for convenience

    Tensor<3, data_t> i_F;

    FOR3(i,j,k)
    {
        i_F[i][j][k] = d1.fspatial[j][k][i] - vars.fbar[j] * metric_vars.K_tensor[i][k] - vars.fbar[k] * metric_vars.K_tensor[i][j];
        FOR1(l)
        {
            i_F[i][j][k] += -chris_phys.ULL[l][i][j] * vars.fspatial[l][k] - chris_phys.ULL[l][i][k] * vars.fspatial[j][l]; 
        }
        
    }
    //and its first derivative
    Tensor<4, data_t> d1_i_F;
    FOR3(i,j,k)
    {
        FOR1(l)
        {
            d1_i_F[i][j][k][l] = d2.fspatial[j][k][i][l] - d1.fbar[j][l] * metric_vars.K_tensor[i][k] - vars.fbar[j] * metric_vars.d1_K_tensor[i][k][l]
                                                        - d1.fbar[k][l] * metric_vars.K_tensor[i][j] - vars.fbar[k] * metric_vars.d1_K_tensor[i][j][l];  
            FOR1(m)
            {
                d1_i_F[i][j][k][l] += - metric_vars.d1_chris_phys[m][i][j][l] * vars.fspatial[m][k] - chris_phys.ULL[m][i][j] * d1.fspatial[m][k][l]
                                    - metric_vars.d1_chris_phys[m][i][k][l] * vars.fspatial[j][m] - chris_phys.ULL[m][i][k] * d1.fspatial[j][m][l]; 
            }
            
        }
    }

    //Projections of the Riemann term, R_alpha mu beta nu f^alpha beta
    //First, define the spatial covariant derivative of the extrinsic curvature Kij
    Tensor<3, data_t> cd1_K_tensor;
    FOR3(i,j,k)
    {
        cd1_K_tensor[i][j][k] = metric_vars.d1_K_tensor[i][j][k];
        FOR1(l)
        {
            cd1_K_tensor[i][j][k] += -chris_phys.ULL[l][k][i] * metric_vars.K_tensor[l][j] -chris_phys.ULL[l][k][j] * metric_vars.K_tensor[i][l];
        }
    }

    data_t scalarRiemannTerm; 
    FOR3(i,j,k)
    {
        scalarRiemannTerm = 0.0;
        FOR1(l)
        {
            scalarRiemannTerm += metric_vars.K * gamma_UU[i][k] * gamma_UU[j][l] * metric_vars.K_tensor[k][l] * vars.fspatial[i][j]
                                + gamma_UU[i][k] * gamma_UU[j][l] * vars.fspatial[i][j] * metric_vars.ricci_phys[k][l];
            FOR1(m)
            {
                //scalarRiemannTerm += 0.0;//gamma_UU[i][k] * gamma_UU[j][l] * vars.fspatial[i][j] * metric_vars.riemann_phys_ULLL[m][k][m][l];
                FOR1(n)
                {
                    scalarRiemannTerm += - gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * metric_vars.K_tensor[i][n] * metric_vars.K_tensor[l][m] * vars.fspatial[j][k];
                }
            }
        }   
    }

    Tensor<1, data_t> vectorRiemannTerm;

    FOR1(i)
    {
        vectorRiemannTerm[i] = 0.0;
        FOR2(j,k)
        {
            vectorRiemannTerm[i] += -gamma_UU[j][k] * metric_vars.K_tensor[i][j] * metric_vars.K * vars.fbar[k]
                                    - gamma_UU[j][k] * metric_vars.ricci_phys[i][k] * vars.fbar[j];
            
            FOR1(l)
            {
                //vectorRiemannTerm[i] += 0.0;//-gamma_UU[j][k] * metric_vars.riemann_phys_ULLL[l][i][l][k] * vars.fbar[j];
                
                FOR1(m)
                {
                 vectorRiemannTerm[i] += gamma_UU[j][l] * gamma_UU[k][m] * metric_vars.K_tensor[j][k] * metric_vars.K_tensor[i][m] * vars.fbar[l]
                 + gamma_UU[j][l] * gamma_UU[k][m] * vars.fspatial[l][m] * (cd1_K_tensor[j][k][i] - cd1_K_tensor[i][j][k]);   
                }
            }
        }
    }

    Tensor<2, data_t> tensorRiemannTerm; 

    FOR2(i,j)
    {
        tensorRiemannTerm[i][j] = vars.fhat * metric_vars.K_tensor[i][j] * metric_vars.K
                                    + vars.fhat * metric_vars.ricci_phys[i][j];
        FOR1(k)
        {
            //tensorRiemannTerm[i][j] += 0.0;//vars.fhat * metric_vars.riemann_phys_ULLL[k][i][k][j];
            FOR1(l)
            {
                tensorRiemannTerm[i][j] += gamma_UU[k][l] * vars.fbar[l] * (2.0 * cd1_K_tensor[i][j][k] - cd1_K_tensor[i][k][j] - cd1_K_tensor[j][k][i])
                - vars.fhat * gamma_UU[k][l] * metric_vars.K_tensor[i][l] * metric_vars.K_tensor[j][k];
                FOR2(m,n)
                {
                    tensorRiemannTerm[i][j] += gamma_UU[k][m] * gamma_UU[l][n] * (metric_vars.K_tensor[i][j] * metric_vars.K_tensor[k][l]
                                                                                - metric_vars.K_tensor[i][k] * metric_vars.K_tensor[j][l]) * vars.fspatial[m][n];
                    FOR1(o)
                    {
                    tensorRiemannTerm[i][j] += gamma_UU[k][m] * gamma_UU[l][n] * metric_vars.gamma[i][o] * vars.fspatial[m][n] * metric_vars.riemann_phys_ULLL[o][k][j][l];
                    }
                }
            }
        }
    }
    //Evolution equations for the field and the conjugate variables:
    /*
      Field Variables
      Name        Rank
      fspatial    2
      fbar        1    
      fhat        0
    
      Conjugate Variables
      Name        Rank
      u           2
      v           2     "main"
      p           1
      q           1     "main"
      w           0     "main"
    */
    //No summation indices
    //Need to define the relations between u and v, p and q.

    //  u <-> v
    Tensor<2, data_t> i_u;
    
    FOR2(i,j)
    {
        i_u[i][j] = d1.fbar[j][i] - vars.fhat * metric_vars.K_tensor[i][j];

        FOR1(k)
        {
            i_u[i][j] += -chris_phys.ULL[k][j][i] * vars.fbar[k]; 

            FOR1(l)
            {
                i_u[i][j] += -gamma_UU[l][k] * vars.fspatial[j][k] * metric_vars.K_tensor[i][l];
            }
        }
        //TESTING
        //i_u[i][j] = 0.0;
    }

    // d1_u <-> d1_v
    
    Tensor<3, data_t> d1_i_u; 
    FOR3(i,j,k)
    {
        d1_i_u[i][j][k] = -d1.fhat[k] * metric_vars.K_tensor[i][j] - vars.fhat * metric_vars.d1_K_tensor[i][j][k] + d2.fbar[j][i][k]; 

        FOR1(l)
        {
            d1_i_u[i][j][k] += -chris_phys.ULL[l][j][i] * d1.fbar[l][k] - metric_vars.d1_chris_phys[l][j][i][k] * vars.fbar[l];

            FOR1(m)
            {
                d1_i_u[i][j][k] += -metric_vars.d1_gamma_UU[l][m][k] * vars.fspatial[j][m] * metric_vars.K_tensor[i][l] - gamma_UU[l][m] * d1.fspatial[j][m][k] * metric_vars.K_tensor[i][l]
                - gamma_UU[l][m] * vars.fspatial[j][m] * metric_vars.d1_K_tensor[i][l][k];
            }
        }
        
    }

    Tensor<1, data_t> i_p;
    //  p <-> q
    FOR1(i)
    {
        i_p[i] = d1.fhat[i]; 

        FOR2(j,k)
        {
            i_p[i] += -2.0 * gamma_UU[j][k] * vars.fbar[j] * metric_vars.K_tensor[i][k];
        }
        
    }
    Tensor<2, data_t> d1_i_p;
    // d1_p <-> d1_q
    FOR2(i,j)
    {
        d1_i_p[i][j] = d2.fhat[i][j]; 

        FOR2(k,l)
        {
            d1_i_p[i][j] += -2.0 * (metric_vars.d1_gamma_UU[l][k][j] * vars.fbar[l] * metric_vars.K_tensor[i][k] 
            + gamma_UU[l][k] * d1.fbar[l][j] * metric_vars.K_tensor[i][k] + gamma_UU[l][k] * vars.fbar[l] * metric_vars.d1_K_tensor[i][k][j]); 
        }
        
    }
    


    //EVOLUTION                             
    rhs.fhat = advec.fhat - metric_vars.lapse * vars.w;

    rhs.w = advec.w + metric_vars.lapse * metric_vars.K * vars.w + metric_vars.lapse * temp_mass * temp_mass * vars.fhat - 2.0 * metric_vars.lapse * scalarRiemannTerm;

    //1 summation index
    FOR1(i)
    {   
        rhs.fbar[i] = advec.fbar[i] - vars.fhat * metric_vars.d1_lapse[i] - metric_vars.lapse * vars.q[i];

        rhs.q[i] = advec.q[i] + metric_vars.lapse * vars.q[i] * metric_vars.K - vars.w * metric_vars.d1_lapse[i] + metric_vars.lapse * temp_mass * temp_mass * vars.fbar[i] - 2.0 * metric_vars.lapse * vectorRiemannTerm[i];

        
       //2 summation indices 
        FOR1(j)
        {   
            rhs.fhat += - 2.0 * gamma_UU[i][j] * vars.fbar[j] * metric_vars.d1_lapse[i];

            rhs.fbar[i] += vars.fbar[j] * metric_vars.d1_shift[j][i];

            rhs.fspatial[i][j] = advec.fspatial[i][j] - (vars.fbar[j] * metric_vars.d1_lapse[i] + vars.fbar[i] * metric_vars.d1_lapse[j])
            - metric_vars.lapse * vars.v[i][j]; 

            rhs.w += - 2.0 * gamma_UU[i][j] * vars.q[j] * metric_vars.d1_lapse[i] - metric_vars.lapse * gamma_UU[i][j] * d1_i_p[j][i]
            - gamma_UU[i][j] * i_p[j] * metric_vars.d1_lapse[i];

            rhs.q[i] += vars.q[j] * metric_vars.d1_shift[j][i];

            rhs.v[i][j] = advec.v[i][j] - vars.q[i] * metric_vars.d1_lapse[j] - vars.q[j] * metric_vars.d1_lapse[i] + metric_vars.lapse * temp_mass * temp_mass * vars.fspatial[i][j]
            + metric_vars.lapse * metric_vars.K * vars.v[i][j] - 2.0 * metric_vars.lapse * tensorRiemannTerm[i][j]; 

            //3 summation indices 
            FOR1(k)
            {   
                rhs.fbar[i] += - gamma_UU[j][k] * vars.fspatial[i][j] * metric_vars.d1_lapse[k] - metric_vars.lapse * gamma_UU[j][k] * vars.fbar[j] * metric_vars.K_tensor[k][i];

                rhs.fspatial[i][j] += vars.fspatial[k][j] * metric_vars.d1_shift[k][i] + vars.fspatial[i][k] * metric_vars.d1_shift[k][j];

                rhs.w += metric_vars.lapse * gamma_UU[i][j] * chris_phys.ULL[k][i][j] * i_p[k];

                rhs.q[i] += - gamma_UU[j][k] * (vars.v[j][i] + i_u[j][i]) * metric_vars.d1_lapse[k] + metric_vars.lapse * gamma_UU[j][k] * (i_p[j] - vars.q[j]) * metric_vars.K_tensor[i][k]
                - metric_vars.lapse * gamma_UU[j][k] * d1_i_u[j][i][k];

                rhs.v[i][j] += vars.v[k][j] * metric_vars.d1_shift[k][i] + vars.v[i][k] * metric_vars.d1_shift[k][j];

                //4 summation indices
                FOR1(l)
                {   
                    rhs.fspatial[i][j] += -metric_vars.lapse * gamma_UU[k][l] * (vars.fspatial[k][j] * metric_vars.K_tensor[i][l] + vars.fspatial[k][i] * metric_vars.K_tensor[j][l]);

                    rhs.w += 2.0 * metric_vars.lapse * gamma_UU[i][k] * gamma_UU[j][l] * i_u[i][j] * metric_vars.K_tensor[k][l];

                    rhs.q[i] += metric_vars.lapse * gamma_UU[j][k] * (chris_phys.ULL[l][k][j] * i_u[l][i] + chris_phys.ULL[l][k][i] * i_u[j][l]);

                    rhs.v[i][j] += - gamma_UU[k][l] * i_F[k][i][j] * metric_vars.d1_lapse[l] - metric_vars.lapse * gamma_UU[k][l] * d1_i_F[k][i][j][l]
                    - gamma_UU[k][l] * metric_vars.lapse * (vars.v[k][i] - i_u[k][i]) * metric_vars.K_tensor[l][j] - gamma_UU[k][l] * metric_vars.lapse * (vars.v[k][j] - i_u[k][j]) * metric_vars.K_tensor[l][i];
                    //5 summation indices 
                    FOR1(m)
                    {
                        rhs.q[i] += metric_vars.lapse * gamma_UU[j][l] * gamma_UU[k][m] * i_F[j][k][i] * metric_vars.K_tensor[l][m];

                        rhs.v[i][j] += metric_vars.lapse * gamma_UU[k][l] * (chris_phys.ULL[m][l][k] * i_F[m][i][j] + chris_phys.ULL[m][l][i] * i_F[k][m][j] + chris_phys.ULL[m][l][j] * i_F[k][i][m]);
                    
                        //rhs.fhat = 0;
                        //rhs.fbar[i] = 0.0;
                        //rhs.fspatial[i][j] = 0.0;
                        //rhs.w = 0.0;
                        //rhs.q[i] = 0.0;
                        //rhs.v[i][j] = 0.0;
                    }

                }
            
            }
         
        }
    
    }

}

/*
rhs.fhat = advec.fhat - 2.0 * gamma_UU[i][j] * vars.fbar[j] * metric_vars.d1_lapse[i] - metric_vars.lapse * vars.w; 

rhs.fbar[i] = advec.fbar[i] + vars.fbar[j] * metric_vars.d1_shift[j][i] - gamma_UU[j][k] * vars.fspatial[i][j] * metric_vars.d1_lapse[k]
- metric_vars.lapse * gamma_UU[j][k] * vars.fbar[j] * metric_vars.K_tensor[i][k] - vars.fhat * metric_vars.d1_lapse[i] - metric_vars.lapse * vars.q[i]; 

rhs.fspatial[i][j] = advec.fspatial[i][j] + vars.fspatial[k][j] * metric_vars.d1_shift[k][i] + vars.fspatial[i][k] * metric_vars.d1_shift[k][j]
- metric_vars.lapse * gamma_UU[k][l] * vars.fspatial[k][j] * metric_vars.K_tensor[i][l] - metric_vars.lapse * gamma_UU[k][l] * vars.fspatial[k][i] * metric_vars.K_tensor[j][l]
- (vars.fbar[j] * metric_vars.d1_lapse[i] + vars.fbar[i] * metric_vars.d1_lapse[j])
- metric_vars.lapse * vars.v[i][j]; 

rhs.w = advec.w - 2.0 * gamma_UU[i][j] * vars.q[j] * metric_vars.d1_lapse[i] + 2.0 * metric_vars.lapse * gamma_UU[i][k] * gamma_UU[j][l] * i_u[i][j] * metric_vars.K_tensor[k][l]
- metric_vars.lapse * gamma_UU[i][j] * d1_i_p[j][i] + metric_vars.lapse * gamma_UU[i][j] * chris_phys.ULL[k][i][j] * i_p[k]
- gamma_UU[i][j] * i_p[j] * metric_vars.d1_lapse[i] + metric_vars.lapse * metric_vars.K * vars.w + metric_vars.lapse * temp_mass * temp_mass * vars.fhat;

rhs.q[i] = advec.q[i] + vars.q[j] * metric_vars.d1_shift[j][i] + metric_vars.lapse * gamma_UU[j][l] * gamma_UU[k][m] * i_F[j][k][i] * metric_vars.K_tensor[l][m]
- gamma_UU[j][k] * (vars.v[j][i] + i_u[j][i]) * metric_vars.d1_lapse[k] + metric_vars.lapse * gamma_UU[j][k] * (i_p[j] - vars.q[j]) * metric_vars.K_tensor[i][k]
+ metric_vars.lapse * vars.q[i] * metric_vars.K - vars.w * metric_vars.d1_lapse[i] + metric_vars.lapse * temp_mass * temp_mass * vars.fbar[i]
- gamma_UU[j][k] * d1_i_u[j][i][k] + gamma_UU[j][k] * (chris_phys.ULL[l][k][j] * i_u[k][i] + chris_phys.ULL[l][k][i] * i_u[j][l]);

rhs.v[i][j] = advec.v[i][j] + vars.v[k][j] * metric_vars.d1_shift[k][i] + vars.v[i][k] * metric_vars.d1_shift[k][j]
- gamma_UU[k][l] * i_F[k][i][j] * metric_vars.d1_lapse[l] - metric_vars.lapse * gamma_UU[k][l] * d1_i_F[k][i][j][l]
+ metric_vars.lapse * gamma_UU[k][l] * (chris_phys.ULL[m][l][k] * i_F[m][i][j] + chris_phys.ULL[m][l][i] * i_F[k][m][j] + chris_phys.ULL[m][l][j] * i_F[k][i][m])
- gamma_UU[k][l] * metric_vars.lapse * (vars.v[k][i] - i_u[k][i]) * metric_vars.K_tensor[l][j] - gamma_UU[k][l] * metric_vars.lapse * (vars.v[k][j] - i_u[k][j]) * metric_vars.K_tensor[l][i]
- vars.q[i] * metric_vars.d1_lapse[j] - vars.q[j] * metric_vars.d1_lapse[i] + metric_vars.lapse * temp_mass * temp_mass * vars.fspatial[i][j]
+ metric_vars.lapse * metric_vars.K * vars.v[i][j]; 
*/

#endif /* FIXEDBGTENSORFIELD_IMPL_HPP_ */
