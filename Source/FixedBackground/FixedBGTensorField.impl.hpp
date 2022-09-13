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

    data_t scalarRiemannTerm = 0.0; 
    FOR3(i,j,k)
    {
       
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
    
    //For constraint damping.
    data_t primaryScalar = -temp_mass * temp_mass * vars.fhat;
    
    FOR2(i,j)
    {
        primaryScalar += -gamma_UU[i][j] * (d1.q[i][j] - d1_i_p[i][j]);

        FOR1(k)
        {
        primaryScalar += gamma_UU[i][j] * (chris_phys.ULL[k][i][j] * vars.q[k] - chris_phys.ULL[k][i][j] * i_p[k]);

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
                            + gamma_UU[j][k] * metric_vars.ricci_phys[i][j] * vars.fbar[k];

            FOR1(l)
            {
                primaryVector[i] += gamma_UU[j][k] * (chris_phys.ULL[l][k][i] * vars.v[j][l] + chris_phys.ULL[l][k][j] * vars.v[l][i]
                                                    -chris_phys.ULL[l][k][i] * i_u[j][l]    - chris_phys.ULL[l][k][j] * i_u[l][i]);
            } 
        }
    }     

    data_t transverseScalar = 0.0;
    
    transverseScalar = metric_vars.K * vars.fhat + vars.w;
    FOR2(i,j)
    {
        transverseScalar += -gamma_UU[i][j] * d1.fbar[i][j];
        FOR1(k)
        {
        transverseScalar += gamma_UU[i][j] * chris_phys.ULL[k][j][i] * vars.fbar[k];
        FOR1(l)
        {
            transverseScalar += metric_vars.K_tensor[i][j] * vars.fspatial[k][l] * gamma_UU[i][k] * gamma_UU[j][l]; 
        }
        }
    }
    
    Tensor<1, data_t> transverseVector;
    
    FOR1(i)
    {
        transverseVector[i] = vars.q[i];//-metric_vars.K * vars.fbar[i] - vars.q[i];
        FOR2(j,k)
        {
        //transverseVector[i] += -gamma_UU[j][k] * metric_vars.K_tensor[i][k] * vars.fbar[j];
        transverseVector[i] += -gamma_UU[j][k] * d1.fspatial[j][i][k];
        FOR1(l)
        {
            transverseVector[i] += gamma_UU[j][k] * (chris_phys.ULL[l][k][i] * vars.fspatial[j][l] + chris_phys.ULL[l][k][j] * vars.fspatial[l][i]);
        }
        }
    }

    data_t traceFieldConstraint = vars.fhat;
    data_t traceMomentumConstraint = vars.w;
    FOR2(i,j)
    {
        traceFieldConstraint -= gamma_UU[i][j] * vars.fspatial[i][j];
        traceMomentumConstraint -= gamma_UU[i][j] * vars.v[i][j];
    }

    


    rhs.fhat = advec.fhat - metric_vars.lapse * vars.w;

    //rhs.fhat += -metric_vars.lapse * vars.thetahat;

    rhs.w = advec.w + metric_vars.lapse * metric_vars.K * vars.w + metric_vars.lapse * temp_mass * temp_mass * vars.fhat - 2.0 * metric_vars.lapse * scalarRiemannTerm;

    //1 summation index
    FOR1(i)
    {   
        rhs.fbar[i] = advec.fbar[i] - vars.fhat * metric_vars.d1_lapse[i] - metric_vars.lapse * vars.q[i];

        //rhs.fbar[i] += -metric_vars.lapse * vars.thetaspatial[i];

        rhs.q[i] = advec.q[i] + metric_vars.lapse * vars.q[i] * metric_vars.K - vars.w * metric_vars.d1_lapse[i] + metric_vars.lapse * temp_mass * temp_mass * vars.fbar[i] + 2.0 * metric_vars.lapse * vectorRiemannTerm[i];

        //rhs.q[i] += -metric_vars.lapse * d1.thetahat[i];
        
       //2 summation indices 
        FOR1(j)
        {   
            
            
            
            rhs.fhat += - 2.0 * gamma_UU[i][j] * vars.fbar[j] * metric_vars.d1_lapse[i];

            rhs.fbar[i] += vars.fbar[j] * metric_vars.d1_shift[j][i];

            rhs.fspatial[i][j] = advec.fspatial[i][j] - (vars.fbar[j] * metric_vars.d1_lapse[i] + vars.fbar[i] * metric_vars.d1_lapse[j])
            - metric_vars.lapse * vars.v[i][j]; 

            //rhs.fspatial[i][j] += -metric_vars.lapse * (d1.thetaspatial[i][j] + d1.thetaspatial[j][i]);

            rhs.w += - 2.0 * gamma_UU[i][j] * vars.q[j] * metric_vars.d1_lapse[i] - metric_vars.lapse * gamma_UU[i][j] * d1_i_p[j][i]
            - gamma_UU[i][j] * i_p[j] * metric_vars.d1_lapse[i];

            rhs.q[i] += vars.q[j] * metric_vars.d1_shift[j][i];

            rhs.v[i][j] = advec.v[i][j] - vars.q[i] * metric_vars.d1_lapse[j] - vars.q[j] * metric_vars.d1_lapse[i] + metric_vars.lapse * temp_mass * temp_mass * vars.fspatial[i][j]
            + metric_vars.lapse * metric_vars.K * vars.v[i][j] - 2.0 * metric_vars.lapse * tensorRiemannTerm[i][j]; 

            //rhs.v[i][j] += metric_vars.lapse * d1.thetaspatial[i][j];

            //3 summation indices 
            FOR1(k)
            {   
                

                rhs.fbar[i] += - gamma_UU[j][k] * vars.fspatial[i][j] * metric_vars.d1_lapse[k] - metric_vars.lapse * gamma_UU[j][k] * vars.fbar[j] * metric_vars.K_tensor[k][i];

                rhs.fspatial[i][j] += vars.fspatial[k][j] * metric_vars.d1_shift[k][i] + vars.fspatial[i][k] * metric_vars.d1_shift[k][j];

                //rhs.fspatial[i][j] += 2.0 * metric_vars.lapse * chris_phys.ULL[k][i][j] * vars.thetaspatial[k];

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
                    
                     //TRACE SUBTRACTION? 
                    //rhs.fspatial[i][j] -= (1.0 / 3.0) * metric_vars.gamma[i][j] * gamma_UU[k][l] * vars.fspatial[k][l];

                    //rhs.v[i][j] -= (1.0 / 3.0) * metric_vars.gamma[i][j] * gamma_UU[k][l] * vars.v[k][l];

                    
                    
                    //5 summation indices 
                    FOR1(m)
                    {
                        rhs.q[i] += metric_vars.lapse * gamma_UU[j][l] * gamma_UU[k][m] * i_F[j][k][i] * metric_vars.K_tensor[l][m];

                        rhs.v[i][j] += metric_vars.lapse * gamma_UU[k][l] * (chris_phys.ULL[m][l][k] * i_F[m][i][j] + chris_phys.ULL[m][l][i] * i_F[k][m][j] + chris_phys.ULL[m][l][j] * i_F[k][i][m]);
                    
                    }

                }
            
            }
         
        }
    
    }

    
    //EVOLUTION     

    

    //Z4 evolution
    /*
    rhs.thetahat = advec.thetahat - metric_vars.lapse *(damping_param_1 * vars.thetahat - primaryScalar);

    FOR1(i)
    {
        rhs.thetaspatial[i] = advec.thetaspatial[i] - metric_vars.lapse * (damping_param_1 * vars.thetaspatial[i] - primaryVector[i]);
            

        FOR1(j)
        {
            rhs.thetaspatial[i] += vars.thetaspatial[j] * metric_vars.d1_shift[j][i];
        }
    }
    */
    //Damping method 1
    /*
    rhs.w += - metric_vars.lapse * (2.0 * damping_param_1 - damping_param_2) * vars.thetahat;

    FOR1(i)
    {
        rhs.q[i] += metric_vars.lapse * (d1.thetahat[i] + damping_param_2 * vars.thetaspatial[i]);

        FOR1(j)
        {
            rhs.v[i][j] += metric_vars.lapse * (-d1.thetaspatial[i][j] - d1.thetaspatial[j][i] - damping_param_2 * metric_vars.gamma[i][j] * vars.thetahat);

            FOR1(k)
            {
                rhs.q[i] += metric_vars.lapse * gamma_UU[j][k] * metric_vars.K_tensor[i][j] * vars.thetaspatial[k];

                rhs.v[i][j] += metric_vars.lapse * 2.0 * (chris_phys.ULL[k][i][j] * vars.thetaspatial[k]);
            }
        }
    }
    */
    /*
    //Damping method 2 (variable reduction)
    rhs.thetahat = - metric_vars.lapse * damping_param_1 * (2.0 + damping_param_2) * vars.thetahat;
    FOR1(i)
    {
        rhs.thetaspatial[i] = metric_vars.lapse * d1.thetahat[i] - vars.thetahat * metric_vars.d1_lapse[i] - metric_vars.lapse * damping_param_1 * vars.thetaspatial[i];
        
        FOR1(j)
        {
            rhs.thetahat += metric_vars.lapse * gamma_UU[i][j] * d1.thetaspatial[i][j]
                            - gamma_UU[i][j] * vars.thetaspatial[i] * metric_vars.d1_lapse[j];


            FOR1(k)
            {
                rhs.thetahat += -metric_vars.lapse * gamma_UU[i][j] * chris_phys.ULL[k][i][j] * vars.thetaspatial[k];

                rhs.thetaspatial[i] += - gamma_UU[j][k] * vars.v[i][k] * metric_vars.d1_lapse[j];

                FOR1(l)
                {
                    rhs.thetahat += -2.0 * metric_vars.lapse * metric_vars.ricci_phys[i][j] * vars.fspatial[k][l] * gamma_UU[i][k] * gamma_UU[j][l];
                }
            }
        }
    }
    */
    /*
    FOR1(i)
    {
        FOR1(j)
        {
            rhs.v[i][j] += - metric_vars.lapse * (d1.thetaspatial[i][j] + d1.thetaspatial[j][i]) - damping_param_1 * (1.0 + damping_param_2) * metric_vars.lapse * metric_vars.gamma[i][j] * vars.thetahat
            + 2.0 * metric_vars.lapse * metric_vars.K_tensor[i][j] * vars.thetahat;
            FOR1(k)
            {
                rhs.v[i][j] += 2.0 * metric_vars.lapse * chris_phys.ULL[k][i][j] * vars.thetaspatial[k];
            }
        }
        //Switch off the evolution of these fields
        rhs.fbar[i] = 0.0;
        rhs.q[i] = 0.0;
    }
    */
    /*
    FOR1(i)
    {
        rhs.fbar[i] = 0.0;
        rhs.q[i] = 0.0;
    }
    */
    //Damping method 3
    //Damping field evolution

    //data_t damping_param_1 = 1.0;
    //data_t damping_param_2 = 0.0;

    //Fifth damping implementation
    
    /*
    rhs.fhat += -1.0 * metric_vars.lapse * (traceFieldConstraint);
    rhs.w += -1.0 * metric_vars.lapse * (transverseScalar + traceMomentumConstraint);
    FOR1(i)
    {
        rhs.q[i] += -1.0 * metric_vars.lapse * transverseVector[i];
    }
    */
    /*
    rhs.thetahat = 0.5 * metric_vars.lapse * (primaryScalar - damping_param_1 * (2.0 + damping_param_2) * vars.thetahat);
    
    FOR1(i)
    {
        rhs.thetaspatial[i] = 0.5 * metric_vars.lapse * (primaryVector[i] - damping_param_1 * vars.thetaspatial[i]);// - vars.thetahat * metric_vars.d1_lapse[i];

        //FOR1(j)
        //{
        //    //rhs.thetahat += - gamma_UU[i][j] * vars.thetaspatial[i] * metric_vars.d1_lapse[j];
        //}
    }
    */
  
    ///////
    /*
    rhs.w += 0.5 * metric_vars.lapse * primaryScalar;

    FOR1(i)
    {
        rhs.q[i] += 0.5 * metric_vars.lapse * primaryVector[i];
        
        FOR1(j)
        {
            rhs.v[i][j] += metric_vars.lapse * (-d1.thetaspatial[i][j] -d1.thetaspatial[j][i] 
                - metric_vars.gamma[i][j] * damping_param_1 * damping_param_2 * vars.thetahat);

            FOR1(k)
            {
                rhs.v[i][j] += 2.0 * metric_vars.lapse * (chris_phys.ULL[k][i][j] * vars.thetaspatial[k]);
                FOR1(l)
                {
                    rhs.v[i][j] += 2.0 * metric_vars.gamma[i][j] * gamma_UU[k][l] * vars.thetaspatial[k] * metric_vars.d1_lapse[l];
                }
            }
        }
    } 
    */
    /*
    rhs.w += 0.5 * abs(metric_vars.lapse) * primaryScalar;

    FOR1(i)
    {
        rhs.q[i] += 0.5 * abs(metric_vars.lapse) * primaryVector[i];
        
        FOR1(j)
        {
            rhs.v[i][j] += abs(metric_vars.lapse) * (-d1.thetaspatial[i][j] -d1.thetaspatial[j][i] 
                - metric_vars.gamma[i][j] * damping_param_1 * damping_param_2 * vars.thetahat);

            FOR1(k)
            {
                rhs.v[i][j] += 2.0 * abs(metric_vars.lapse) * (chris_phys.ULL[k][i][j] * vars.thetaspatial[k]);
                FOR1(l)
                {
                    rhs.v[i][j] += 2.0 * metric_vars.gamma[i][j] * gamma_UU[k][l] * vars.thetaspatial[k] * metric_vars.d1_lapse[l];
                }
            }
        }
    }
    */
    /*
    data_t lie_n_thetahat;
    Tensor<1, data_t> lie_n_thetaspatial;

    lie_n_thetahat = metric_vars.lapse * (primaryScalar - damping_param_1 * (2.0 + damping_param_2) * vars.thetahat);    

    FOR1(i)
    {
        lie_n_thetaspatial[i] = metric_vars.lapse * (primaryVector[i] - damping_param_1 * vars.thetaspatial[i] + d1.thetahat[i]) - vars.thetahat * metric_vars.d1_lapse[i];// / metric_vars.lapse; 
        FOR1(j)
        {
            lie_n_thetahat += metric_vars.lapse * gamma_UU[i][j] * d1.thetaspatial[i][j] - gamma_UU[i][j] * vars.thetaspatial[i] * metric_vars.d1_lapse[j];
            FOR1(k)
            {
                lie_n_thetahat += - metric_vars.lapse * gamma_UU[i][j] * chris_phys.ULL[k][i][j] * vars.thetaspatial[k];
            }
        }

    }
    
    rhs.thetahat = lie_n_thetahat;// + advec.thetahat;

    FOR1(i)
    {
        rhs.thetaspatial[i] = lie_n_thetaspatial[i];// + advec.thetaspatial[i];

        //FOR1(j)
        //{
           // rhs.thetaspatial[i] += vars.thetaspatial[j] * metric_vars.d1_shift[j][i];
        //}
    }
    */
    
    
    //Modifications to the evolution of the momentum fields.
    
    /*
    rhs.w += metric_vars.lapse  * primaryScalar;

    FOR1(i)
    {
       rhs.q[i] += metric_vars.lapse * primaryVector[i];

       FOR1(j)
       {
           rhs.v[i][j] += metric_vars.lapse * (-d1.thetaspatial[i][j] -d1.thetaspatial[j][i]) + 
           abs(metric_vars.lapse) * metric_vars.gamma[i][j] * (primaryScalar - 2.0 * damping_param_1 * (1.0 + damping_param_2) * vars.thetahat);

           FOR1(k)
           {
               rhs.v[i][j] += 2.0 * metric_vars.lapse * (chris_phys.ULL[k][i][j] * vars.thetaspatial[k]);

               FOR1(l)
               {
                   rhs.v[i][j] += metric_vars.lapse * metric_vars.gamma[i][j] * (2.0 * gamma_UU[k][l] * d1.thetaspatial[k][l]);

                   FOR1(m)
                   {
                       rhs.v[i][j] += - 2.0 * metric_vars.lapse * metric_vars.gamma[i][j] * gamma_UU[k][l] * chris_phys.ULL[m][k][l] * vars.thetaspatial[m];
                   }
               }
           }
       }
    }
    */
    /*
    //rhs.w += metric_vars.lapse  * primaryScalar;
    
    FOR1(i)
    {
        
       //rhs.q[i] += metric_vars.lapse * primaryVector[i];

       FOR1(j)
       {
           rhs.v[i][j] += metric_vars.lapse * (-d1.thetaspatial[i][j] -d1.thetaspatial[j][i]) + 
           metric_vars.lapse * metric_vars.gamma[i][j] * (- 2.0 * damping_param_1 * (1.0 + damping_param_2) * vars.thetahat);

           FOR1(k)
           {
               rhs.v[i][j] += 2.0 * metric_vars.lapse * (chris_phys.ULL[k][i][j] * vars.thetaspatial[k]);

               FOR1(l)
               {
                   rhs.v[i][j] += metric_vars.lapse * metric_vars.gamma[i][j] * (2.0 * gamma_UU[k][l] * d1.thetaspatial[k][l]);

                   FOR1(m)
                   {
                       rhs.v[i][j] += - 2.0 * metric_vars.lapse * metric_vars.gamma[i][j] * gamma_UU[k][l] * chris_phys.ULL[m][k][l] * vars.thetaspatial[m];
                   }
               }
           }
       }  
    }
    */




    //rhs.fhat = 0.0;
    //rhs.w = 0.0;
    

    

    //Evolution (rewritten) Can ignore all terms with the shift in too
    /*
    rhs.fhat = -metric_vars.lapse * vars.w;
    rhs.w = -2.0 * metric_vars.lapse * scalarRiemannTerm;

    FOR1(i)
    {

        rhs.fbar[i] = -metric_vars.lapse * vars.q[i] - vars.fhat * metric_vars.d1_lapse[i];
        rhs.q[i] = -vars.w * metric_vars.d1_lapse[i] - 2.0 * metric_vars.lapse * vectorRiemannTerm[i];

        FOR1(j)
        {
            rhs.fhat += -2.0 * gamma_UU[i][j] * vars.fbar[j] * metric_vars.d1_lapse[i];
            rhs.fspatial[i][j] = -metric_vars.lapse * vars.v[i][j] -vars.fbar[i] * metric_vars.d1_lapse[j] -vars.fbar[j] * metric_vars.d1_lapse[i];

            rhs.w += -2.0 * gamma_UU[i][j] * vars.q[i] * metric_vars.d1_lapse[j] - gamma_UU[i][j] * i_p[i] * metric_vars.d1_lapse[j]
                    -metric_vars.lapse * gamma_UU[i][j] * d1_i_p[i][j];
            rhs.v[i][j] = -vars.q[i] * metric_vars.d1_lapse[j] -vars.q[j] * metric_vars.d1_lapse[i] - 2.0 * metric_vars.lapse * tensorRiemannTerm[i][j];
        
            FOR1(k)
            {
                rhs.fbar[i] += -gamma_UU[j][k] * vars.fspatial[i][j] * metric_vars.d1_lapse[k];
                rhs.w += metric_vars.lapse * gamma_UU[i][j] * chris_phys.ULL[k][i][j] * i_p[k];

                rhs.q[i] += -gamma_UU[j][k] * vars.v[i][j] * metric_vars.d1_lapse[k] - gamma_UU[j][k] * i_u[j][i] * metric_vars.d1_lapse[k]
                - gamma_UU[j][k] * metric_vars.lapse * d1_i_u[j][i][k];

                FOR1(l)
                {
                    rhs.q[i] += gamma_UU[j][k] * metric_vars.lapse * (chris_phys.ULL[l][k][j] * i_u[l][i] + chris_phys.ULL[l][k][i] * i_u[j][l]);
                    rhs.v[i][j] += -gamma_UU[k][l] * i_F[k][i][j] * metric_vars.d1_lapse[l] - metric_vars.lapse * d1_i_F[k][i][j][l];

                    FOR1(m)
                    {
                        rhs.v[i][j] += gamma_UU[k][l] * metric_vars.lapse * (chris_phys.ULL[m][l][k] * i_F[m][i][j] + chris_phys.ULL[m][l][i] * i_F[k][m][i] + chris_phys.ULL[m][l][j] * i_F[k][i][m]);
                    }
                } 
            }
        }
    }
    */

    //  TURNING OFF EVOLUTION OF FHAT, W FOR PURPOSES OF TRACE ENFORCEMENT 
    /*
    rhs.fhat = 0.0;
    rhs.w = 0.0;
    FOR1(i)
    {
        rhs.fbar[i] = 0.0;
    }
    rhs.fspatial[2][2] = 0.0;
    rhs.v[2][2] = 0.0;
    rhs.q[2] = 0.0;
    */


    //matter_rhs_damping(rhs, vars, metric_vars, d1, d2, advec);
    //rhs.fspatial[2][2] = 0.0;
    //rhs.v[2][2] = 0.0;
    //rhs.q[2] = 0.0;
    //rhs.fbar[2] = 0.0;

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





/*
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void FixedBGTensorField<potential_t>::matter_rhs_damping(
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

    const double damping_coeff = 0.1;

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
    Tensor<3, data_t> cd1_K_tensor;
    FOR3(i,j,k)
    {
        cd1_K_tensor[i][j][k] = metric_vars.d1_K_tensor[i][j][k];
        FOR1(l)
        {
            cd1_K_tensor[i][j][k] += -chris_phys.ULL[l][k][i] * metric_vars.K_tensor[l][j] -chris_phys.ULL[l][k][j] * metric_vars.K_tensor[i][l];
        }
    }

    //NEEDS MASS TERM ADDING
    data_t primaryScalar = 0.0;
    
    FOR1(i)
    {
        primaryScalar += vars.w * metric_vars.shift[i] * metric_vars.d1_lapse[i] / metric_vars.lapse;
        FOR1(j)
        {
        primaryScalar += gamma_UU[i][j] *(i_p[j] - vars.q[j])* metric_vars.d1_lapse[i];
        FOR2(k,l)
        {

        
            primaryScalar += -2.0 * metric_vars.lapse * gamma_UU[i][k] * gamma_UU[j][l] * (metric_vars.K * metric_vars.K_tensor[i][j] + metric_vars.ricci_phys[i][j]) * vars.fspatial[k][l]
                            + metric_vars.lapse * gamma_UU[i][k] * gamma_UU[j][l] * i_u[i][j] * metric_vars.K_tensor[k][l];
                            
            
            
            FOR2(m,n)
            {
            primaryScalar += metric_vars.lapse * gamma_UU[i][l] * gamma_UU[j][m] * gamma_UU[k][n] * (metric_vars.K_tensor[i][m] * metric_vars.K_tensor[j][k] * vars.fspatial[l][n]);
            }
        }
        }
    }
    Tensor<1, data_t> primaryVector;
    FOR1(i)
    {
        primaryVector[i] = vars.w * metric_vars.d1_lapse[i] / metric_vars.lapse;

        FOR1(j)
        {
        
        FOR1(k)
        {
            primaryVector[i] += - gamma_UU[j][k] * metric_vars.K_tensor[k][i] * i_p[j]
                            -2.0*gamma_UU[j][k] * metric_vars.K_tensor[j][i] * metric_vars.K * vars.fbar[k]
                            -2.0*gamma_UU[j][k] * metric_vars.ricci_phys[i][j] * vars.fbar[k];
            
            FOR2(l,m)
            {
            primaryVector[i] += -2.0*gamma_UU[j][l] * gamma_UU[k][m] * (cd1_K_tensor[j][k][i] - cd1_K_tensor[i][j][k]) * vars.fspatial[l][m]
                                +2.0*gamma_UU[j][l] * gamma_UU[k][m] * metric_vars.K_tensor[k][j] * metric_vars.K_tensor[i][m] * vars.fbar[l];

            }
        }
        }
    }


    rhs.theta = metric_vars.lapse * (primaryScalar - damping_coeff * vars.theta) + advec.theta;

    FOR1(i)
    {
        rhs.q[i] += metric_vars.lapse * d1.theta[i];

        rhs.X[i] = metric_vars.lapse * (primaryVector[i] - damping_coeff * vars.X[i]) + advec.X[i];
        FOR1(j)
        {
            rhs.v[i][j] += metric_vars.lapse * (d1.X[i][j] + d1.X[j][i]);

            rhs.X[i] += vars.X[j] * metric_vars.d1_shift[j][i];
            
            FOR1(k)
            {
                rhs.v[i][j] += -2.0 * metric_vars.lapse * (chris_phys.ULL[k][i][j] * vars.X[k]);
            }
        }
    }

}
*/
#endif /* FIXEDBGTENSORFIELD_IMPL_HPP_ */
