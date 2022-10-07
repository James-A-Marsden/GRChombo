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
    //matter_rhs_damping
    

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

    const double temp_mass = 0.0;

    //Defining the F variable for convenience

    Tensor<3, data_t> i_F;

    FOR3(i,j,k)
    {
        i_F[i][j][k] = d1.fspatial[j][k][i];
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
            d1_i_F[i][j][k][l] = d2.fspatial[j][k][i][l];  
            FOR1(m)
            {
                d1_i_F[i][j][k][l] += - metric_vars.d1_chris_phys[m][i][j][l] * vars.fspatial[m][k] - chris_phys.ULL[m][i][j] * d1.fspatial[m][k][l]
                                    - metric_vars.d1_chris_phys[m][i][k][l] * vars.fspatial[j][m] - chris_phys.ULL[m][i][k] * d1.fspatial[j][m][l]; 
            }
            
        }
    }

    data_t scalarRiemannTerm = 0.0; 
    FOR3(i,j,k)
    {
       
        FOR1(l)
        {
            scalarRiemannTerm += gamma_UU[i][k] * gamma_UU[j][l] * vars.fspatial[i][j] * metric_vars.ricci_phys[k][l];
        }   
    }

    Tensor<1, data_t> vectorRiemannTerm;

    FOR1(i)
    {
        vectorRiemannTerm[i] = 0.0;
        FOR2(j,k)
        {
            vectorRiemannTerm[i] += - gamma_UU[j][k] * metric_vars.ricci_phys[i][k] * vars.fbar[j];
        }
    }

    Tensor<2, data_t> tensorRiemannTerm; 

    FOR2(i,j)
    {
        tensorRiemannTerm[i][j] = - vars.fhat * metric_vars.ricci_phys[i][j] / metric_vars.lapse; 
        FOR3(k,l,m)
        {
            FOR2(n,o)
            {
              tensorRiemannTerm[i][j] += gamma_UU[k][m] * gamma_UU[l][n] * metric_vars.gamma[i][o] * vars.fspatial[m][n] * metric_vars.riemann_phys_ULLL[o][k][j][l];
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
        i_u[i][j] = d1.fbar[j][i];

        FOR1(k)
        {
            i_u[i][j] += -chris_phys.ULL[k][j][i] * vars.fbar[k]; 
        }
        //TESTING
        //i_u[i][j] = 0.0;
    }

    // d1_u <-> d1_v
    
    Tensor<3, data_t> d1_i_u; 
    FOR3(i,j,k)
    {
        d1_i_u[i][j][k] = d2.fbar[j][i][k]; 

        FOR1(l)
        {
            d1_i_u[i][j][k] += -chris_phys.ULL[l][j][i] * d1.fbar[l][k] - metric_vars.d1_chris_phys[l][j][i][k] * vars.fbar[l];
        }
        
    }

    Tensor<1, data_t> i_p;
    //  p <-> q
    FOR1(i)
    {
        i_p[i] = vars.fhat * metric_vars.d1_lapse[i] / metric_vars.lapse / metric_vars.lapse - d1.fhat[i] / metric_vars.lapse;
    }
    Tensor<2, data_t> d1_i_p;
    // d1_p <-> d1_q
    FOR2(i,j)
    {
        //d1_i_p[i][j] = d2.fhat[i][j]; 
        d1_i_p[i][j] = -2.0 * vars.fhat * metric_vars.d1_lapse[i] * metric_vars.d1_lapse[j] * pow(metric_vars.lapse, -3.0)
        + (d1.fhat[j] * metric_vars.d1_lapse[i] + d1.fhat[i] * metric_vars.d1_lapse[j] + vars.fhat * metric_vars.d2_lapse[i][j]) * pow(metric_vars.lapse, -2.0) 
        - d2.fhat[i][j] / metric_vars.lapse; 
        
    }
    
    //For constraint damping.
    //data_t primaryScalar = -temp_mass * temp_mass * vars.fhat;
    data_t primaryScalar = temp_mass * temp_mass * vars.fhat / metric_vars.lapse;
    
    FOR2(i,j)
    {
        primaryScalar += -gamma_UU[i][j] * (d1.q[i][j] - (-2.0 * vars.fhat * metric_vars.d1_lapse[i] * metric_vars.d1_lapse[j] * pow(metric_vars.lapse, -3.0)
        + (d1.fhat[j] * metric_vars.d1_lapse[i] + d1.fhat[i] * metric_vars.d1_lapse[j] + vars.fhat * metric_vars.d2_lapse[i][j]) * pow(metric_vars.lapse, -2.0) 
        - d2.fhat[i][j] / metric_vars.lapse));

        FOR1(k)
        {
        primaryScalar += gamma_UU[i][j] * (chris_phys.ULL[k][i][j] * vars.q[k] - chris_phys.ULL[k][i][j] * (vars.fhat * metric_vars.d1_lapse[k] / metric_vars.lapse / metric_vars.lapse - d1.fhat[k] / metric_vars.lapse));

            FOR1(l)
            {
                primaryScalar += gamma_UU[i][k] * gamma_UU[j][l] * metric_vars.ricci_phys[i][j] * vars.fspatial[k][l];
            }
        }
    } 
     
    Tensor<1, data_t> primaryVector;
    FOR1(i)
    {
        primaryVector[i] = -temp_mass * temp_mass * vars.fbar[i];

        FOR2(j,k)
        {
            primaryVector[i] += -gamma_UU[j][k] * (d1.v[j][i][k] - d1_i_u[j][i][k]) 
                            //+= -gamma_UU[j][k] * (d1.v[j][i][k] - d2.fbar[i][j][k])
                            //d1_i_u[i][j][k] += -chris_phys.ULL[l][j][i] * d1.fbar[l][k] - metric_vars.d1_chris_phys[l][j][i][k] * vars.fbar[l]
                            
                            + gamma_UU[j][k] * metric_vars.ricci_phys[i][j] * vars.fbar[k];

            FOR1(l)
            {
                primaryVector[i] +=  gamma_UU[j][k] * (-chris_phys.ULL[l][i][j] * d1.fbar[l][k] - metric_vars.d1_chris_phys[l][i][j][k] * vars.fbar[l]);
                primaryVector[i] += gamma_UU[j][k] * (chris_phys.ULL[l][k][i] * vars.v[j][l] + chris_phys.ULL[l][k][j] * vars.v[l][i]
                                                    -chris_phys.ULL[l][k][i] * d1.fbar[l][j] - chris_phys.ULL[l][k][j] * d1.fbar[i][l]);
                FOR1(m)
                {
                    primaryVector[i] += gamma_UU[j][k] * (chris_phys.ULL[l][k][i] * chris_phys.ULL[m][j][i] * vars.fbar[m] + chris_phys.ULL[l][k][j] * chris_phys.ULL[m][j][i] * vars.fbar[m]);
                }


                                                    //-chris_phys.ULL[l][k][i] * i_u[j][l]    - chris_phys.ULL[l][k][j] * i_u[l][i]);
            } 
        }
    }     

    data_t transverseScalar = 0.0;
    
    transverseScalar = vars.w;
    FOR2(i,j)
    {
        transverseScalar += -gamma_UU[i][j] * d1.fbar[i][j];
        FOR1(k)
        {
        transverseScalar += gamma_UU[i][j] * chris_phys.ULL[k][j][i] * vars.fbar[k];
        }
    }
    
    Tensor<1, data_t> transverseVector;
    
    FOR1(i)
    {
        transverseVector[i] = vars.q[i];
        FOR2(j,k)
        {
      
        transverseVector[i] += -gamma_UU[j][k] * d1.fspatial[j][i][k];
        FOR1(l)
        {
            transverseVector[i] += gamma_UU[j][k] * (chris_phys.ULL[l][k][i] * vars.fspatial[j][l] + chris_phys.ULL[l][k][j] * vars.fspatial[l][i]);
        }
        }
    }

    //data_t traceFieldConstraint = vars.fhat;
    data_t traceFieldConstraint = -vars.fhat / metric_vars.lapse;
    data_t traceMomentumConstraint = vars.w;
    FOR2(i,j)
    {
        traceFieldConstraint -= gamma_UU[i][j] * vars.fspatial[i][j];
        traceMomentumConstraint -= gamma_UU[i][j] * vars.v[i][j];
    }

    


    
    rhs.fhat = metric_vars.lapse * metric_vars.lapse * vars.w;

    rhs.w = - temp_mass * temp_mass * vars.fhat - 2.0 * metric_vars.lapse * scalarRiemannTerm;

    //1 summation index
    FOR1(i)
    {   
        
        rhs.fbar[i] = vars.fhat * metric_vars.d1_ln_lapse[i] - metric_vars.lapse * vars.q[i];

        rhs.q[i] = - vars.w * metric_vars.d1_lapse[i] + metric_vars.lapse * temp_mass * temp_mass * vars.fbar[i] + 2.0 * metric_vars.lapse * vectorRiemannTerm[i];

       
        
       //2 summation indices 
        FOR1(j)
        {   
            
            //rhs.fhat += - 2.0 * gamma_UU[i][j] * vars.fbar[j] * metric_vars.d1_lapse[i];

            rhs.fhat += + 2.0 * gamma_UU[i][j] * metric_vars.lapse * vars.fbar[j] * metric_vars.d1_lapse[i];

            rhs.fspatial[i][j] = - (vars.fbar[j] * metric_vars.d1_lapse[i] + vars.fbar[i] * metric_vars.d1_lapse[j])
            - metric_vars.lapse * vars.v[i][j]; 

            rhs.w += gamma_UU[i][j] * ( -2.0 * vars.q[j] * metric_vars.d1_lapse[i] - d1.fhat[i] * metric_vars.d1_ln_lapse[j] + d2.fhat[i][j] - vars.fhat * metric_vars.d2_ln_lapse[i][j]);
            
            //rhs.v[i][j] = - vars.q[i] * metric_vars.d1_lapse[j] - vars.q[j] * metric_vars.d1_lapse[i] + metric_vars.lapse * temp_mass * temp_mass * vars.fspatial[i][j]
            // - 2.0 * metric_vars.lapse * tensorRiemannTerm[i][j]; 
            rhs.v[i][j] = - metric_vars.lapse * temp_mass * temp_mass * vars.fspatial[i][j] - 2.0 * metric_vars.lapse * tensorRiemannTerm[i][j]; 
            

            //3 summation indices 
            FOR1(k)
            {   
                
                rhs.fbar[i] += - gamma_UU[j][k] * vars.fspatial[i][j] * metric_vars.d1_lapse[k];

                rhs.w += gamma_UU[i][j] * (-chris_phys.ULL[k][j][i] * d1.fhat[k] + vars.fhat * chris_phys.ULL[k][j][i] * metric_vars.d1_ln_lapse[k]);

                rhs.q[i] += - gamma_UU[j][k] * (vars.v[j][i] + d1.fbar[i][j]) * metric_vars.d1_lapse[k];
                            - metric_vars.lapse * gamma_UU[j][k] * d2.fbar[i][j][k];

                
                //4 summation indices
                FOR1(l)
                {   
                    
                    //rhs.q[i] += metric_vars.lapse * gamma_UU[j][k] * (chris_phys.ULL[l][k][j] * i_u[l][i] + chris_phys.ULL[l][k][i] * i_u[j][l]);
                    rhs.q[i] += gamma_UU[j][k] * chris_phys.ULL[l][i][j] * vars.fbar[l] * metric_vars.d1_lapse[k]
                                -metric_vars.lapse * gamma_UU[j][k] * (-chris_phys.ULL[l][k][j] * d1.fbar[i][l] -chris_phys.ULL[l][k][i] * d1.fbar[l][j]
                                                                        -metric_vars.d1_chris_phys[l][j][i][k] * vars.fbar[l] - chris_phys.ULL[l][j][i] * d1.fbar[l][k]);

                    rhs.v[i][j] += //- gamma_UU[k][l] * i_F[k][i][j] * metric_vars.d1_lapse[l] - metric_vars.lapse * gamma_UU[k][l] * d1_i_F[k][i][j][l];
                                    - gamma_UU[k][l] * d1.fspatial[i][j][k] * metric_vars.d1_lapse[l]
                                    - gamma_UU[k][l] * metric_vars.lapse * d2.fspatial[i][j][k][l];
                     
                    rhs.v[i][j] += -gamma_UU[k][l] * (d1.fspatial[i][k][l] * metric_vars.d1_lapse[j] + d1.fspatial[j][k][l] * metric_vars.d1_lapse[i]);
                     
                    //5 summation indices 
                    FOR1(m)
                    {
                        rhs.q[i] += -metric_vars.lapse * gamma_UU[j][k] * (chris_phys.ULL[m][k][j] * chris_phys.ULL[l][m][i] * vars.fbar[l] + chris_phys.ULL[m][k][i] * chris_phys.ULL[l][m][j] * vars.fbar[l]);

                        rhs.v[i][j] += gamma_UU[k][l] * (chris_phys.ULL[m][k][i] * vars.fspatial[m][j] + chris_phys.ULL[m][k][j] * vars.fspatial[i][m]) * metric_vars.d1_lapse[l];
                        
                        rhs.v[i][j] += -gamma_UU[k][l] * metric_vars.lapse * (- chris_phys.ULL[m][l][k] * d1.fspatial[i][j][m] - chris_phys.ULL[m][l][i] * d1.fspatial[m][j][k] - chris_phys.ULL[m][l][j] * d1.fspatial[i][m][k]
                                                          - metric_vars.d1_chris_phys[m][k][i][l] * vars.fspatial[m][j] - chris_phys.ULL[m][k][i] * d1.fspatial[m][j][l]
                                                          - metric_vars.d1_chris_phys[m][j][k][l] * vars.fspatial[i][m] - chris_phys.ULL[m][j][k] * d1.fspatial[i][m][l]);
                        
                        rhs.v[i][j] += gamma_UU[k][l] * ((chris_phys.ULL[m][l][i] * vars.fspatial[m][k] + chris_phys.ULL[m][l][k] * vars.fspatial[i][m]) * metric_vars.d1_lapse[j]
                                                        +(chris_phys.ULL[m][l][j] * vars.fspatial[m][k] + chris_phys.ULL[m][l][k] * vars.fspatial[j][m]) * metric_vars.d1_lapse[i]);
                        
                        //6 summation indices 
                        FOR1(n)
                        {
                            rhs.v[i][j] += -gamma_UU[k][l] * metric_vars.lapse * (chris_phys.ULL[n][l][k] * chris_phys.ULL[m][n][i] * vars.fspatial[m][j] + chris_phys.ULL[n][l][i] * chris_phys.ULL[m][k][n] * vars.fspatial[m][j] + chris_phys.ULL[n][l][j] * chris_phys.ULL[m][k][i] * vars.fspatial[m][n]
                                                                                + chris_phys.ULL[n][l][k] * chris_phys.ULL[m][n][j] * vars.fspatial[i][m] + chris_phys.ULL[n][l][j] * chris_phys.ULL[m][k][n] * vars.fspatial[i][m] + chris_phys.ULL[n][l][i] * chris_phys.ULL[m][k][j] * vars.fspatial[n][m]);
                        }
                    }

                }
            
            }
         
        }
    
    }
    rhs.w = 0.0;
    rhs.fhat = 0.0;
    FOR1(i)
    {
        rhs.fbar[i] = 0.0;
        rhs.q[i] = 0.0;
    }
}
#endif /* FIXEDBGTENSORFIELD_IMPL_HPP_ */
