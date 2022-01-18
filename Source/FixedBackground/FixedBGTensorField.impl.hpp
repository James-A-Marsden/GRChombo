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
    out.S = 0;//TensorAlgebra::compute_trace(out.Sij, gamma_UU);

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

    const double temp_mass = 0;
    const double temp_G = 0; 
    const double temp_pi = 3.14159;
    const double temp_d1_Ktensor = 0;
    //const double temp_K_trace =0;
    const double temp_d1_K_trace =0;
    //const double temp_F 0; 
    const double temp_d1_F =0;
    const double temp_T =0;

    //Define some quantities
    data_t Tr_K;
    FOR2(i,j){Tr_K = gamma_UU[i][j] * metric_vars.K_tensor[i][j];}

    Tensor<3, data_t> F_tensor;

    FOR3(i,j,k)
    {
        F_tensor[i][j][k] = d1.fspatial[j][k][i] - d1.fspatial[i][k][j] - d1.fspatial[j][i][k] - 2 * vars.fbar[i] * metric_vars.K_tensor[j][k];
        FOR1(l)
        {
            F_tensor[i][j][k] += 2 * chris_phys.ULL[l][j][k] * vars.fspatial[i][l]; 
        }
    }

    /*



    //Field evolution equations 
    rhs.fhat = advec.fhat + metric_vars.lapse * vars.w;


    //Conjugate field evolution equations
    rhs.w = advec.w + metric_vars.lapse * temp_mass * temp_mass * vars.fhat + metric_vars.lapse * metric_vars.K * vars.w;
    //+ 8 * metric_vars.lapse * temp_pi * temp_G * (out.S + out.rho) * vars.fhat;

    //1 explicit index
    FOR1(i)
    {   
   //     rhs.fbar[i] = advec.fbar[i] + 0.5 * (vars.fhat * metric_vars.d1_lapse[i] + metric_vars.lapse * d1.fhat[i])
   //     + 0.5 * (vars.fhat * metric_vars.d1_lapse[i] - metric_vars.lapse * vars.p[i]);

        //Conjugate 

      //  rhs.q[i] = advec.q[i] - vars.fhat * metric_vars.d1_lapse[i]
        //-alpha * \partial_i K
      //  - metric_vars.lapse * temp_mass * temp_mass * vars.fbar[i];


       //2 explicit indices 
        FOR1(j)
        {   
      //      rhs.fspatial[i][j] = advec.fspatial[i][j] + vars.fbar[j] * metric_vars.d1_lapse[i] + vars.fbar[i] * metric_vars.d1_lapse[j]
       //     + metric_vars.lapse * (d1.fbar[i][j] + d1.fbar[j][i]) + 2 * metric_vars.lapse * vars.fhat * metric_vars.K_tensor[i][j] - metric_vars.lapse * vars.v[i][j];
            
       //     rhs.fhat += 2 * gamma_UU[i][j] * vars.fbar[i] * metric_vars.d1_lapse[j];

            //Conjugate

         //   rhs.w +=  -gamma_UU[i][j] * (vars.p[j] + 2 * vars.q[j]) * metric_vars.d1_lapse[i];

            
            
           //3 explicit indices 
            FOR1(k)
            {   
               // rhs.fspatial[i][j] += 2 * metric_vars.lapse * chris_phys.ULL[k][i][j] * vars.fbar[k];
                
               // rhs.fbar[i] += gamma_UU[j][k] * vars.fbar[j] * metric_vars.d1_shift[k][i]
            //    +  gamma_UU[j][k] * vars.fspatial[i][j] * metric_vars.d1_lapse[k];

                //Conjugate 
                //rhs.w += - metric_vars.lapse * gamma_UU[i][j] * (d1.p[j][i] - chris_phys.ULL[k][i][j] * vars.p[k]);
                
               // rhs.q[i] += metric_vars.lapse * gamma_UU[j][k] * (d1.u[j][i][k]) //contains A 

               // rhs.q[i] += vars.v[i][j] * gamma_UU[j][k] * metric_vars.d1_lapse[k] 
             //   + gamma_UU[j][k] * metric_vars.d1_lapse[k] * vars.u[j][i]; //Check to make sure this A << is the correct way round with its indices
             //   rhs.q[i] += metric_vars.lapse * vars.q[i] * gamma_UU[j][k] metric_vars.K_tensor[j][k] 
             //   + metric_vars.lapse * (vars.p[k] - vars.q[k]) * gamma_UU[k][j] * metric_vars.K_tensor[j][i];
            //    rhs.q[i] += metric_vars.lapse * gamma_UU[j][k];// * (derivatives of K here + christoffels - add)

                //4 explicit indices
                FOR1(l)
                {   
                    //rhs.fspatial[i][j] += gamma_UU[k][l] * (vars.fspatial[k][j] * metric_vars.d1_shift[i][l] + vars.fspatial[i][k] * metric_vars.d1_shift[j][l]);  


                    //Conjugate
                    //FIX THE d1.K_tensor !
                    //rhs.w += 2 * metric_vars.lapse * gamma_UU[i][k] * gamma_UU[j][l] * ( vars.v[i][j] * metric_vars.K_tensor[k][l])
                    //- 2 * metric_vars.lapse * 
                    //- 2 * metric_vars.lapse * gamma_UU[k][l] * gamma_UU[i][j] * d1.K_tensor[k][l][j] * vars.fbar[i]; 
                    
                    //rhs.q[i] += metric_vars.lapse * gamma_UU[j][k] * (-chris_phys.ULL[l][k][j] vars.u[l][i] - chris_phys.ULL[l][k][i] * vars.u[j][l]);
                    //5 explicit indices 
                    FOR1(m)
                    {
                    //Conjugate 
                    
                    //rhs.w  += -2 * metric_vars.lapse * (- gamma_UU[i][k] * gamma_UU[j][l] * (d1.K_tensor[k][l][i] 
                   // - chris_phys.ULL[m][i][k] * metric_vars.K_tensor[m][l] - chris_phys.ULL[m][i][l] * metric_vars.K_tensor[k][m]) * vars.fbar[j]);
                    }   

                        FOR1(n)
                        {
                            const test_var;
                        }

                    }

                }
            
            }
         
        }
    
    }
    */

}





//PI 



//Leave F as a separate thing?...
/*
rhs.v[i][j] = advec.v[i][j] + (vars.v[k][j] * metric_vars.d1_shift[k][i] + vars.v[i][k] * metric_vars.d1_shift[k][j])
+ Tr_K * vars.v[i][j] - gamma_UU[k][l] * (F_tensor[k][i][j] * metric_vars.d1_lapse[l] + metric_vars.lapse * temp_d1_F)
+ gamma_UU[k][l] * metric_vars.lapse * (chris_phys.ULL[m][l][k] * F_tensor[m][i][j] + chris_phys.ULL[m][l][i] * F_tensor[k][m][j] + chris_phys.ULL[m][l][j] * F_tensor[k][i][m])
- metric_vars.lapse * gamma_UU[k][l] * metric_vars.K_tensor[j][k] * 2 * (vars.v[l][i] + gamma_UU[m][n] * vars.fspatial[l][m] * metric_vars.K_tensor[i][n] - d1.fbar[i][l] + chris_phys.ULL[m][l][i] vars.fbar[m] - 2 * vars.fhat * metric_vars.K_tensor[l][i])
- metric_vars.lapse * gamma_UU[k][l] * metric_vars.K_tensor[i][k] * 2 * (vars.v[l][j] + gamma_UU[m][n] * vars.fspatial[l][m] * metric_vars.K_tensor[j][n] - d1.fbar[j][l] + chris_phys.ULL[m][l][j] vars.fbar[m] - 2 * vars.fhat * metric_vars.K_tensor[l][j])
- vars.q[i] * metric_vars.d1_lapse[j] - vars.q[j] * metric_vars.d1_lapse[i] - 8 * temp_pi * temp_G * metric_vars.lapse * gamma_UU[k][l] * (vars.S_ij[l][i] * vars.fspatial[j][k] + vars.S_ij[l][j] * vars.fspatial[i][k])
+ 8 * temp_G * temp_pi * metric_vars.lapse * temp_T * vars.fspatial[i][j]
- metric_vars.lapse * gamma_UU[k][l] * (temp_d1_Ktensor - chris_phys.ULL[m][k][l] * metric_vars.K_tensor[m][i] - chris_phys.ULL[m][k][i] * metric_vars.K_tensor[l][m]) * vars.fbar[j]
- metric_vars.lapse * gamma_UU[k][l] * (temp_d1_Ktensor - chris_phys.ULL[m][k][l] * metric_vars.K_tensor[m][j] - chris_phys.ULL[m][k][j] * metric_vars.K_tensor[l][m]) * vars.fbar[i]
+ metric_vars.lapse * (vars.fbar[j] temp_d1_K_trace + vars.fbar[i] * temp_d1_K_trace) + metric_vars.lapse * temp_mass * temp_mass * vars.fspatial[i][j];

rhs.q[i] = advec.q[i] + vars.q[i] * gamma_UU[j][k] * metric_vars.d1_shift[k][i]
- 2 * (- gamma_UU[j][k] * vars.fspatial[l][j] * metric_vars.K_tensor[i][k] * gamma_UU[l][m] * metric_vars.d1_lapse[m]
+ d1.fbar[i][j] * gamma_UU[j][k] * metric_vars.d1_lapse[k]
- chris_phys.ULL[j][k][i] * vars.fbar[j] * gamma_UU[k][l] * metric_vars.d1_lapse[l]
+ vars.fhat * metric_vars.K_tensor[j][i] * gamma_UU[j][k] * metric_vars.d1_lapse[k])
+ metric_vars.lapse * gamma_UU[j][k] * (d1.u[j][i][k])
+ metric_vars.lapse * gamma_UU[j][k] * (-chris_phys.ULL[l][k][j] * vars.u[l][i] - chris_phys.ULL[l][k][i] * vars.u[j][l])
+ metric_vars.lapse * vars.q[i] * Tr_K * metric_vars.lapse * gamma_UU[j][k] * (vars.p[k] - vars.q[k]) * metric_vars.K_tensor[j][i]
- vars.w * metric_vars.d1_lapse[i] - gamma_UU[j][k] * gamma_UU[l][m] * F_tensor[j][l] * metric_vars.K_tensor[k][m]
+ 8 * temp_pi * metric_vars.lapse * temp_G * gamma_UU[j][k] * vars.S_ij[k][i] * vars.fbar[j]
- 4 * temp_pi * metric_vars.lapse * temp_G * ( (vars.Tr_S_ij - vars.rho) * vars.fbar[i])
+ metric_vars.lapse * gamma_UU[j][k] * temp_d1_Ktensor - metric_vars.lapse * gamma_UU[j][k] * (chris_phys.ULL[l][j][k] * metric_vars.K_tensor[l][i] + chris_phys.ULL[l][j][i] * metric_vars.K_tensor[k][l])
- metric_vars.lapse * temp_d1_K_trace - gamma_UU[j][k] * gamma_UU[l][m] * (temp_d1_Ktensor) * vars.fspatial[l][i]
+ gamma_UU[j][k] * gamma_UU[l][m] * (chris_phys.ULL[n][j][k] * metric_vars.K_tensor[n][m] + chris_phys.ULL[n][j][m] * metric_vars.K_tensor[k][n]) * vars.fspatial[l][i]
- gamma_UU[j][k] * vars.fspatial[j][i] * temp_d1_K_trace
- 4 * metric_vars.lapse * temp_pi * temp_G * (vars.Tr_S_ij + vars.rho) * vars.fbar[i] - metric_vars.lapse * temp_mass * temp_mass * vars.fbar[i];


rhs.w = advec.w + 2 * metric_vars.lapse * gamma_UU[i][k] * gamma_UU[j][l] * metric_vars.K_tensor[i][j] * (vars.v[k][l]) - gamma_UU[i][m]vars.f[k][i] * metric_vars.K_tensor[l][m]
- (d1.fbar[l][k] - chris_phys.ULL[m][l][k] * vars.fbar[m])
- metric_vars.lapse * gamma_UU[i][j] * (d1.p[j][i] - chris_phys.ULL[k][i][j] * vars.p[k])
- gamma_UU[i][j] * vars.p[j] metric_vars.d1_lapse[i]
- 2 * gamma_UU[i][j] * vars.q[j] * metric_vars.d1_lapse[i]
+ metric_vars.lapse * Tr_K * vars.w
+ 2 * metric_vars.lapse * gamma_UU[i][k] * gamma_UU[j][l] * (temp_d1_Ktensor - chris_phys.ULL[m][i][k] * metric_vars.K_tensor[m][l] - chris_phys.ULL[m][i][l] * metric_vars.K_tensor[k][m]) * vars.fbar[j]
- 2 * metric_vars.lapse * gamma_UU[i][j] * vars.fbar[i] * temp_d1_K_trace
+ 8 * metric_vars.lapse * temp_G * temp_pi * (vars.Tr_S_ij + vars.rho) * vars.fbar
+ metric_vars.lapse * temp_mass * temp_mass * vars.fhat;
*/


#endif /* FIXEDBGTENSORFIELD_IMPL_HPP_ */
