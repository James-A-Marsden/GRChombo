/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef KERRSCHILDFIXEDBG_HPP_
#define KERRSCHILDFIXEDBG_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

//! Class which computes the initial conditions for a Kerr Schild BH
//! https://arxiv.org/pdf/gr-qc/9805023.pdf
class KerrSchildFixedBG
{
  public:
    //! Struct for the params of the  BH
    struct params_t
    {
        double mass = 1.0;                      //!<< The mass of the BH
        std::array<double, CH_SPACEDIM> center; //!< The center of the BH
        double spin = 0.0;                      //!< The spin param a = J / M
    };

    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;

    const params_t m_params;
    const double m_dx;

    KerrSchildFixedBG(params_t a_params, double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
        // check this spin param is sensible
        if ((m_params.spin > m_params.mass) || (m_params.spin < -m_params.mass))
        {
            MayDay::Error(
                "The dimensionless spin parameter must be in the range "
                "-1.0 < spin < 1.0");
        }
    }

    /// This just calculates chi which helps with regridding, debug etc
    /// it is only done once on setup as the BG is fixed
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // get position and set vars
        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        Vars<data_t> metric_vars;
        compute_metric_background(metric_vars, current_cell);

        // calculate and save chi
        data_t chi = TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
        chi = pow(chi, -1.0 / 3.0);
        current_cell.store_vars(chi, c_chi);
        current_cell.store_vars(metric_vars.K, c_Kout);
        current_cell.store_vars(metric_vars.d1_K[0], c_d1Kout1);
        current_cell.store_vars(metric_vars.d1_K[1], c_d1Kout2);
        current_cell.store_vars(metric_vars.d1_K[2], c_d1Kout3);
        
        /*
        current_cell.store_vars(metric_vars.d2_shift[0][0][0], c_d2_shift111);
        current_cell.store_vars(metric_vars.d2_shift[0][0][1], c_d2_shift112);
        current_cell.store_vars(metric_vars.d2_shift[0][0][2], c_d2_shift113);
        current_cell.store_vars(metric_vars.d2_shift[0][1][0], c_d2_shift121);
        current_cell.store_vars(metric_vars.d2_shift[0][1][1], c_d2_shift122);
        current_cell.store_vars(metric_vars.d2_shift[0][1][2], c_d2_shift123);
        current_cell.store_vars(metric_vars.d2_shift[0][2][0], c_d2_shift131);
        current_cell.store_vars(metric_vars.d2_shift[0][2][1], c_d2_shift132);
        current_cell.store_vars(metric_vars.d2_shift[0][2][2], c_d2_shift133);
        */
        /*
        current_cell.store_vars(metric_vars.d1_K_tensor[0][0][0], c_d1_K_tensor111);
        current_cell.store_vars(metric_vars.d1_K_tensor[0][1][0], c_d1_K_tensor121);
        current_cell.store_vars(metric_vars.d1_K_tensor[0][2][0], c_d1_K_tensor131);
        current_cell.store_vars(metric_vars.d1_K_tensor[1][1][0], c_d1_K_tensor221);
        current_cell.store_vars(metric_vars.d1_K_tensor[1][2][0], c_d1_K_tensor231);
        current_cell.store_vars(metric_vars.d1_K_tensor[2][2][0], c_d1_K_tensor331);

        current_cell.store_vars(metric_vars.d1_K_tensor[0][0][1], c_d1_K_tensor112);
        current_cell.store_vars(metric_vars.d1_K_tensor[0][1][1], c_d1_K_tensor122);
        current_cell.store_vars(metric_vars.d1_K_tensor[0][2][1], c_d1_K_tensor132);
        current_cell.store_vars(metric_vars.d1_K_tensor[1][1][1], c_d1_K_tensor222);
        current_cell.store_vars(metric_vars.d1_K_tensor[1][2][1], c_d1_K_tensor232);
        current_cell.store_vars(metric_vars.d1_K_tensor[2][2][1], c_d1_K_tensor332);
        
        current_cell.store_vars(metric_vars.d1_K_tensor[0][0][2], c_d1_K_tensor113);
        current_cell.store_vars(metric_vars.d1_K_tensor[0][1][2], c_d1_K_tensor123);
        current_cell.store_vars(metric_vars.d1_K_tensor[0][2][2], c_d1_K_tensor133);
        current_cell.store_vars(metric_vars.d1_K_tensor[1][1][2], c_d1_K_tensor223);
        current_cell.store_vars(metric_vars.d1_K_tensor[1][2][2], c_d1_K_tensor233);
        current_cell.store_vars(metric_vars.d1_K_tensor[2][2][2], c_d1_K_tensor333);
        */
        /*
        current_cell.store_vars(metric_vars.d1_gamma_UU[0][0][0], c_d1_gamma_UU111);
        current_cell.store_vars(metric_vars.d1_gamma_UU[0][1][0], c_d1_gamma_UU121);
        current_cell.store_vars(metric_vars.d1_gamma_UU[0][2][0], c_d1_gamma_UU131);
        current_cell.store_vars(metric_vars.d1_gamma_UU[1][1][0], c_d1_gamma_UU221);
        current_cell.store_vars(metric_vars.d1_gamma_UU[1][2][0], c_d1_gamma_UU231);
        current_cell.store_vars(metric_vars.d1_gamma_UU[2][2][0], c_d1_gamma_UU331);
        
        current_cell.store_vars(metric_vars.d1_gamma_UU[0][0][2], c_d1_gamma_UU112);
        current_cell.store_vars(metric_vars.d1_gamma_UU[0][1][2], c_d1_gamma_UU122);
        current_cell.store_vars(metric_vars.d1_gamma_UU[0][2][2], c_d1_gamma_UU132);
        current_cell.store_vars(metric_vars.d1_gamma_UU[1][1][2], c_d1_gamma_UU222);
        current_cell.store_vars(metric_vars.d1_gamma_UU[1][2][2], c_d1_gamma_UU232);
        current_cell.store_vars(metric_vars.d1_gamma_UU[2][2][2], c_d1_gamma_UU332);

        current_cell.store_vars(metric_vars.d1_gamma_UU[0][0][3], c_d1_gamma_UU113);
        current_cell.store_vars(metric_vars.d1_gamma_UU[0][1][3], c_d1_gamma_UU123);
        current_cell.store_vars(metric_vars.d1_gamma_UU[0][2][3], c_d1_gamma_UU133);
        current_cell.store_vars(metric_vars.d1_gamma_UU[1][1][3], c_d1_gamma_UU223);
        current_cell.store_vars(metric_vars.d1_gamma_UU[1][2][3], c_d1_gamma_UU233);
        current_cell.store_vars(metric_vars.d1_gamma_UU[2][2][3], c_d1_gamma_UU333);
        *
        /*
        current_cell.store_vars(metric_vars.ricci_phys[0][0], c_ricci_phys11);
        current_cell.store_vars(metric_vars.ricci_phys[0][1], c_ricci_phys12);
        current_cell.store_vars(metric_vars.ricci_phys[0][2], c_ricci_phys13);
        current_cell.store_vars(metric_vars.ricci_phys[1][1], c_ricci_phys22);
        current_cell.store_vars(metric_vars.ricci_phys[1][2], c_ricci_phys23);
        current_cell.store_vars(metric_vars.ricci_phys[2][2], c_ricci_phys33);
        
        */
    }

    // Kerr Schild solution
    template <class data_t, template <typename> class vars_t>
    void compute_metric_background(vars_t<data_t> &vars,
                                   const Cell<data_t> &current_cell) const
    {
        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

        // black hole params - mass M and spin a
        const double M = m_params.mass;
        const double a = m_params.spin;
        const double a2 = a * a;

        // work out where we are on the grid including effect of spin
        // on x direction (length contraction)
        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const data_t rho = coords.get_radius();
        const data_t rho2 = rho * rho;

        // the Kerr Schild radius r
        const data_t r2 = 0.5 * (rho2 - a2) +
                          sqrt(0.25 * (rho2 - a2) * (rho2 - a2) + a2 * z * z);
        const data_t r = sqrt(r2);
        const data_t cos_theta = z / r;

        // find the H and el quantities (el decomposed into space and time)
        data_t H = M * r / (r2 + a2 * cos_theta * cos_theta);
        const Tensor<1, data_t> el = {(r * x + a * y) / (r2 + a2),
                                      (r * y - a * x) / (r2 + a2), z / r};
        const data_t el_t = 1.0;

        // Calculate the gradients in el and H
        Tensor<1, data_t> dHdx;
        Tensor<1, data_t> dltdx;
        Tensor<2, data_t> dldx;
        Tensor<3, data_t> d2ldx2;
        Tensor<2, data_t> d2ltdx2;
        Tensor<2, data_t> d2Hdx2; 
        get_KS_derivs(dHdx, d2Hdx2, dldx, d2ldx2, dltdx, d2ltdx2, H, coords);

        // populate ADM vars
        vars.lapse = pow(1.0 + 2.0 * H * el_t * el_t, -0.5);
        FOR2(i, j)
        {
            vars.gamma[i][j] =
                TensorAlgebra::delta(i, j) + 2.0 * H * el[i] * el[j];
        }
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse_sym(vars.gamma);
        FOR1(i)
        {
            vars.shift[i] = 0;
            FOR1(j)
            {
                vars.shift[i] += gamma_UU[i][j] * 2.0 * H * el[j] * el_t;
            }
        }
        FOR2(i, j)


        // Calculate partial derivative of spatial metric
        FOR3(i, j, k)
        {
            vars.d1_gamma[i][j][k] =
                2.0 * (el[i] * el[j] * dHdx[k] + H * el[i] * dldx[j][k] +
                       H * el[j] * dldx[i][k]);
        }

                // First define some useful quantities:

         
        //ADDED
        vars.gamma_UU = compute_inverse_sym(vars.gamma);

        Tensor<2, data_t> lambda_UU; 
        //Manually define its contents as there isn't any nice symmetry we can use 
        lambda_UU[0][0] = el[1] * el[1] + el[2] * el[2]; 
        lambda_UU[0][1] = -el[0] * el[1];
        lambda_UU[0][2] = -el[0] * el[2];
        lambda_UU[1][1] = el[0] * el[0] + el[2] * el[2];
        lambda_UU[1][2] = -el[1] * el[2];
        lambda_UU[2][2] = el[0] * el[0] + el[1] * el[1];
        //Symmetries
        lambda_UU[1][0] = lambda_UU[0][1];
        lambda_UU[2][0] = lambda_UU[0][2];
        lambda_UU[2][1] = lambda_UU[1][2];
        //Deivatives of lambda_UU
        Tensor<2, Tensor<1, data_t>> d1_lambda_UU;
        FOR1(i)
        {
            d1_lambda_UU[0][0][i] = 2.0 * dldx[1][i] * el[1] + 2.0 * dldx[2][i] * el[2];
            d1_lambda_UU[0][1][i] = -dldx[0][i] * el[1] - el[0] * dldx[1][i]; 
            d1_lambda_UU[0][2][i] = -dldx[0][i] * el[2] - el[0] * dldx[2][i];
            d1_lambda_UU[1][1][i] = 2.0 * dldx[0][i] * el[0] + 2.0 * dldx[2][i] * el[2];
            d1_lambda_UU[1][2][i] = -dldx[1][i] * el[2] - el[1] * dldx[2][i];
            d1_lambda_UU[2][2][i] = 2.0 * dldx[0][i] * el[0] + 2.0 * dldx[1][i] * el[1]; 
            //Symmetries
            d1_lambda_UU[1][0][i] = d1_lambda_UU[0][1][i];
            d1_lambda_UU[2][0][i] = d1_lambda_UU[0][2][i];
            d1_lambda_UU[2][1][i] = d1_lambda_UU[1][2][i];
        }
        FOR3(i,j,k)
        {
            
            vars.d1_gamma_UU[i][j][k] = -2.0 * dHdx[k] * (TensorAlgebra::delta(i, j) + 2.0 * H * lambda_UU[i][j]) / ((1 + 2.0 * H) * (1 + 2.0 * H))
            + ((2.0 * dHdx[k] * lambda_UU[i][j]) + (2.0 * H * d1_lambda_UU[i][j][k]))/(1 + 2.0 * H);
            
            //vars.d1_gamma_UU[i][j][k] = 0;
            /*
            vars.d1_gamma_UU[i][j][k] = 2.0 * (dHdx[k] * el[i] * el[j] + H * dldx[i][k] * el[j] + H * el[i] * dldx[j][k])
            - 2.0 * pow(vars.lapse,-3) * vars.d1_lapse[k] * vars.shift[i] * vars.shift[j]
            + (vars.d1_shift[i][k] * vars.shift[j] + vars.shift[i] * vars.d1_shift[j][k])/ (vars.lapse * vars.lapse);
            */
        }
        //Second derivative of the spatial metric (can simplify further)
        FOR1(i)
        {
            FOR3(j,k,m)
            {
                vars.d2_gamma[i][j][k][m] = 
                    2.0 * (d2Hdx2[k][m] * el[i] * el[j] + dHdx[k] * dldx[i][m] * el[j] + dHdx[k] * el[i] * dldx[j][m]
                    + dHdx[m] * dldx[i][k] * el[j] + H * d2ldx2[i][k][m] * el[j] + H * dldx[i][k] * dldx[j][m]
                    + dHdx[m] * el[i] * dldx[j][k] + H * dldx[i][m] * dldx[j][k] + H * el[i] * d2ldx2[j][k][m]);
                //vars.d2_gamma[i][j][k][m] =0;
            }

        }
                                                
        //Calculate derivative of the Christoffel symbol (phys)
        FOR2(i, j)
        {
            FOR3(k, m, n)
            {
                vars.d1_chris_phys[i][j][k][m] = vars.d1_gamma_UU[i][n][m] * (vars.d1_gamma[k][n][j] + vars.d1_gamma[n][j][k] - vars.d1_gamma[j][k][n])
                                                + gamma_UU[i][n] * (vars.d2_gamma[k][n][j][m] + vars.d2_gamma[n][j][k][m] - vars.d2_gamma[j][k][n][m])/2.0;
                vars.d1_chris_phys[i][j][k][m] = 0;
            }
        }

        // calculate derivs of lapse and shift
        FOR1(i)
        {
            vars.d1_lapse[i] = -pow(vars.lapse, 3.0) * el_t *
                               (el_t * dHdx[i] + 2.0 * H * dltdx[i]);
        }

        // use the fact that shift^i = lapse^2 * shift_i
        FOR2(i, j)
        {
            vars.d1_shift[i][j] =
                2.0 * el_t * dHdx[j] * pow(vars.lapse, 2.0) * el[i] +
                4.0 * el_t * H * vars.lapse * vars.d1_lapse[j] * el[i] +
                2.0 * el_t * H * pow(vars.lapse, 2.0) * dldx[i][j] +
                2.0 * dltdx[j] * H * pow(vars.lapse, 2.0) * el[i];
        }
        data_t alpha2 = vars.lapse * vars.lapse;
        FOR3(i,j,k)
        {
            /*
            vars.d2_shift[i][j][k] = 4* el[i] * el_t * vars.d2_lapse[j][k] * vars.lapse *H + 
            2* vars.lapse *( vars.lapse *( dHdx[j] * dltdx[k] * el[i]  + 
            dHdx[j] * dldx[i][k] * el_t  +  d2Hdx2[j][k] * el[i] * el_t  + 
            dHdx[k] *( dltdx[j] * el[i]  +  dldx[i][j] * el_t ) + 
            dldx[i][k] * dltdx[j] *H +  dldx[i][j] * dltdx[k] *H + 
            d2ltdx2[j][k] * el[i] *H +  d2ldx2[i][j][k] * el_t *H) + 
            2* vars.d1_lapse[k] *( dHdx[j] * el[i] * el_t  +  dltdx[j] * el[i] *H + 
            dldx[i][j] * el_t *H)) + 
            4* vars.d1_lapse[j] *( el[i] * el_t * vars.d1_lapse[k] *H + 
            vars.lapse *( dHdx[k] * el[i] * el_t  +  dltdx[k] * el[i] *H + 
            dldx[i][k] * el_t *H));
            */
            vars.d2_shift[i][j][k] =
            4.0 * vars.lapse * vars.d1_lapse[k] * dHdx[j] * el[i] * el_t 
            + 2.0 * alpha2 * d2Hdx2[j][k] * el[i] * el_t
            + 2.0 * alpha2 * dHdx[j] * dldx[i][k] * el_t 
            + 2.0 * alpha2 * dHdx[j] * el[i] * dltdx[k]

            + 4.0 * vars.d1_lapse[k] * vars.d1_lapse[j] * H * el[i] * el_t 
            + 4.0 * vars.lapse * vars.d2_lapse[j][k] * H * el[i] * el_t
            + 4.0 * vars.lapse * vars.d1_lapse[j] * dHdx[k] * el[i] * el_t 
            + 4.0 * vars.lapse * vars.d1_lapse[j] * H * dldx[i][k] * el_t
            + 4.0 * vars.lapse * vars.d1_lapse[j] * H * el[i] * dltdx[k]

            + 4.0 * vars.d1_lapse[k] * vars.lapse * H * dldx[i][j] * el_t 
            + 2.0 * alpha2 * dHdx[k] * dldx[i][j] * el_t
            + 2.0 * alpha2 * H * d2ldx2[i][j][k] * el_t 
            + 2.0 * alpha2 * H * dldx[i][j] * dltdx[k]

            + 4.0 * vars.d1_lapse[k] * vars.lapse * H * el[i] * dltdx[j] 
            + 2.0 * alpha2 * dHdx[k] * el[i] * dltdx[j] 
            + 2.0 * alpha2 * H * dldx[i][k] * dltdx[j] 
            + 2.0 * alpha2 * H * el[i] * d2ltdx2[j][k];
        }
        // calculate the extrinsic curvature, using the fact that
        // 2 * lapse * K_ij = D_i \beta_j + D_j \beta_i - dgamma_ij dt
        // and dgamma_ij dt = 0 in chosen fixed gauge
        const auto chris_phys = compute_christoffel(vars.d1_gamma, gamma_UU);
        FOR2(i, j)
        {
            vars.K_tensor[i][j] = 0.0;
            FOR1(k)
            {
                vars.K_tensor[i][j] +=
                    vars.gamma[k][j] * vars.d1_shift[k][i] +
                    vars.gamma[k][i] * vars.d1_shift[k][j] +
                    (vars.d1_gamma[k][i][j] + vars.d1_gamma[k][j][i]) *
                        vars.shift[k];
                FOR1(m)
                {
                    vars.K_tensor[i][j] += -2.0 * chris_phys.ULL[k][i][j] *
                                           vars.gamma[k][m] * vars.shift[m];
                }
            }
            vars.K_tensor[i][j] *= 0.5 / vars.lapse;
        }
        vars.K = compute_trace(gamma_UU, vars.K_tensor);
        FOR3(i, j, k)
        {
            vars.d1_K_tensor[i][j][k] = 0.0;

            FOR1(m)
            {
                vars.d1_K_tensor[i][j][k] += vars.d1_gamma[m][j][k] * vars.d1_shift[m][i] + vars.gamma[m][j] * vars.d2_shift[m][i][k]
                                            + vars.d1_gamma[m][i][k] * vars.d1_shift[m][j] + vars.gamma[m][i] * vars.d2_shift[m][j][k]
                + (vars.d2_gamma[m][i][j][k] + vars.d2_gamma[m][j][i][k]) * vars.shift[m]
                + (vars.d1_gamma[m][i][j] + vars.d1_gamma[m][j][i]) * vars.d1_shift[m][k];
                FOR1(n)
                {
                    vars.d1_K_tensor[i][j][k] += -2.0 * (vars.d1_chris_phys[m][i][j][k] * vars.gamma[m][n] * vars.shift[n]
                    + chris_phys.ULL[m][i][j] * vars.d1_gamma[m][n][k] * vars.shift[n]
                    + chris_phys.ULL[m][i][j] * vars.gamma[m][n] * vars.d1_shift[n][k]);
                }

            }
            vars.d1_K_tensor[i][j][k] *= 0.5 / vars.lapse;
            vars.d1_K_tensor[i][j][k] += - vars.d1_lapse[k] / vars.lapse * vars.K_tensor[i][j];

            //vars.d1_K_tensor[i][j][k] = 0.0;
        }
                //Derivative of the trace, \partial_i K = \partial_i(gamma^jk K_jk)
        FOR1(i)
        {
            vars.d1_K[i] = 0;
            FOR2(j,k)
            {
                vars.d1_K[i] += vars.d1_gamma_UU[j][k][i] * vars.K_tensor[j][k] + gamma_UU[j][k] * vars.d1_K_tensor[j][k][i];
            }
        }

        //spatial riemann curvature tensor 
        
        FOR1(i)
        {
            FOR3(j,k,l)
            {
                vars.riemann_phys_ULLL[i][j][k][l] = vars.d1_chris_phys[i][l][j][k] - vars.d1_chris_phys[i][k][j][l];

                FOR1(m)
                {
                    vars.riemann_phys_ULLL[i][j][k][l] += chris_phys.ULL[m][l][j] * chris_phys.ULL[i][m][k] - chris_phys.ULL[m][k][j] * chris_phys.ULL[i][m][l]; 
                    //vars.riemann_phys_ULLL[i][j][k][l] = 0;
                }
            }
        }

        //spatial ricci tensor
        FOR2(i,j)
        {
            vars.ricci_phys[i][j] = 0;
            FOR1(k)
            {
                vars.ricci_phys[i][j] += vars.riemann_phys_ULLL[k][i][k][j];
            }
        }
    }

  protected:
    /// Work out the gradients of the quantities H and el appearing in the Kerr
    /// Schild solution
    template <class data_t>
    void get_KS_derivs(Tensor<1, data_t> &dHdx, Tensor<2, data_t> &d2Hdx2, Tensor<2, data_t> &dldx, Tensor<3,data_t> &d2ldx2,
                       Tensor<1, data_t> &dltdx, Tensor<2, data_t> &d2ltdx2, const data_t &H,
                       const Coordinates<data_t> &coords) const
    {
        // black hole params - mass M and boost v
        const double M = m_params.mass;
        const double a = m_params.spin;
        const double a2 = a * a;

        // work out where we are on the grid, and useful quantities
        Tensor<1, data_t> x;
        x[0] = coords.x;
        x[1] = coords.y;
        x[2] = coords.z;
        const double z = coords.z;
        const data_t rho = coords.get_radius();
        const data_t rho2 = rho * rho;

        // the Kerr Schild radius r
        const data_t r2 = 0.5 * (rho2 - a2) +
                          sqrt(0.25 * (rho2 - a2) * (rho2 - a2) + a2 * z * z);
        const data_t r = sqrt(r2);
        const data_t cos_theta = z / r;
        const data_t cos_theta2 = cos_theta * cos_theta;

        using namespace TensorAlgebra;
        // derivatives of r wrt actual grid coords
        Tensor<1, data_t> drhodx;
        FOR1(i) { drhodx[i] = x[i] / rho; }
   
        Tensor<2,data_t> d2rhodx2;
        FOR2(i,j) { d2rhodx2[i][j] = delta(i,j) / rho - x[i] * x[j] / rho2; }

        Tensor<1, data_t> drdx;
        FOR1(i)
        {
            drdx[i] =
                0.5 / r *
                (rho * drhodx[i] +
                 0.5 / sqrt(0.25 * (rho2 - a2) * (rho2 - a2) + a2 * z * z) *
                     (drhodx[i] * rho * (rho2 - a2) +
                      delta(i, 2) * 2.0 * a2 * z));
        }

        Tensor<2, data_t> d2rdx2;
        //Denomimator in drdx term, for simplification
        const data_t denom = pow(0.25 * (rho2 - a2) * (rho2 - a2) + a2 * z * z, 0.5);
        Tensor<1, data_t> d1_denom;
        FOR1(i)
        {
            d1_denom[i] = 0.5 * (drhodx[i] * rho * (rho2 - a2) + 2.0 * delta(i,2) * a2 * z) / denom;
        }

        FOR2(i,j)
        {
 
            d2rdx2[i][j] = -0.5/r * drdx[j] * drdx[i] + 0.5/r * (drhodx[i] * drhodx[j] + rho * d2rhodx2[i][j])
            + 0.25/r * (d2rhodx2[i][j] * rho * (rho2 - a2) + drhodx[i] * drhodx[j] * (rho2 - a2) + 2.0 * drhodx[i] * drhodx[j] * rho2
            + 2.0 * delta(i,2) * delta(j,2) * a2) * pow(0.25 * (rho2 - a2) * (rho2 - a2) + a2 * z * z, -0.5)
            -0.25/r * (drhodx[i] * rho * (rho2 - a2) + 2.0 * delta(i,2) * a2 * z) * 0.5 * (drhodx[j] * rho * (rho2 - a2) + 2.0 * delta(j,2) * a2 * z) * pow(0.25 * (rho2 - a2) * (rho2 - a2) + a2 * z * z, -1.5);   
        }

        Tensor<1, data_t> dcosthetadx;
        FOR1(i) { dcosthetadx[i] = -z / r2 * drdx[i] + delta(i, 2) / r; }

        Tensor<2, data_t> d2costhetadx2;
        FOR2(i, j) 
        {
            d2costhetadx2[i][j] = (1.0 / r2) * (- delta(i,2) * drdx[j] - delta(j,2) * drdx[i] + (2.0 / r) * z * drdx[i] * drdx[j] - z * d2rdx2[i][j]); 
        }

        FOR1(i)
        {
            dHdx[i] = H * (drdx[i] / r -
                           2.0 / (r2 + a2 * cos_theta2) *
                               (r * drdx[i] + a2 * cos_theta * dcosthetadx[i]));
        }

        FOR2(i, j)
        {
            d2Hdx2[i][j] = H * (drdx[j] / r - 2.0 / (r2 + a2 * cos_theta2) * (r * drdx[j] + a2 * cos_theta * dcosthetadx[j])) * (drdx[i] / r - 2.0 / (r2 + a2 * cos_theta2) * (r * drdx[i] + a2 * cos_theta * dcosthetadx[i]))
            + H * (d2rdx2[i][j] / r - drdx[i] * drdx[j] / r2 
            + 4.0 / (r2 + a2 * cos_theta2) *  (r * drdx[j] + a2 * dcosthetadx[j] * cos_theta ) * (r * drdx[i] + a2 * dcosthetadx[i] * cos_theta) / (r2 + a2 * cos_theta2)
            - 2.0 * (drdx[j] * drdx[i] + r * d2rdx2[i][j] + a2 * dcosthetadx[i] * dcosthetadx[j] + a2 * cos_theta * d2costhetadx2[i][j]) / (r2 + a2 * cos_theta2));
        }
        // note to use convention as in rest of tensors the last index is the
        // derivative index so these are d_i l_j
        FOR1(i)
        {
            // first the el_x comp
            dldx[0][i] =
                (x[0] * drdx[i] + r * delta(i, 0) + a * delta(i, 1) -
                 2.0 * r * drdx[i] * (r * x[0] + a * x[1]) / (r2 + a2)) /
                (r2 + a2);
            // now the el_y comp
            dldx[1][i] =
                (x[1] * drdx[i] + r * delta(i, 1) - a * delta(i, 0) -
                 2.0 * r * drdx[i] * (r * x[1] - a * x[0]) / (r2 + a2)) /
                (r2 + a2);
            // now the el_z comp
            dldx[2][i] = -x[2] * drdx[i] / r2 + delta(i, 2) / r;
        }

        FOR2(i,j)
        {
            //el_x component
            /*
            d2ldx2[0][i][j] = delta(j,0) * drdx[i] + x[0] * d2rdx2[i][j] + drdx[j] * delta(i,0) +
             ((-2.0 * drdx[i] * drdx[j] - 2.0 * r * d2rdx2[i][j]) * (r * x[0] + a * x[1]) - 2.0 * r * drdx[i] * (drdx[j] * x[0] + r * delta(j,0) + a * delta(j,1) )) / (r2 + a2) / (r2 + a2)
            + 8.0 * r2 * drdx[i] * drdx[j] * (r * x[0] + a * x[1]) * pow((r2 + a2),-3.0);
            */
            d2ldx2[0][i][j] = (delta(j,0) * drdx[i] + drdx[j] * delta(i,0)) / (r2 + a2)
                            - 2.0 * drdx[j] * r * (x[0] * drdx[i] + r * delta(i,0) + a * delta(i,1)) / (r2 + a2) / (r2 + a2) +
                            (-2.0 * drdx[j] * drdx[i] * (r * x[0] + a * x[1]) - 2.0 * r * d2rdx2[i][j] * (r * x[0] + a * x[1]) 
                            - 2.0 * r * drdx[i] * (drdx[j] * x[0] + r * delta(j,0) + a * delta(j,1))) / (r2 + a2) / (r2 + a2)
                            + 4.0 * r * drdx[i] * (r * x[0] + a * x[1]) * 2.0 * r * drdx[j] * pow((r2 + a2),-3.0);
            //el_y component
            /*
            d2ldx2[1][i][j] = delta(j,1) * drdx[i] + x[1] * d2rdx2[i][j] + drdx[j] * delta(i,0) +
            ((-2.0 * drdx[i] * drdx[j] - 2.0 * r * d2rdx2[i][j]) * (r * x[1] - a * x[0]) - 2.0 * r * drdx[i] * (drdx[j] * x[1] + r * delta(j,1) - a * delta(j,0) )) / (r2 + a2) / (r2 + a2)
            + 8.0 * r2 * drdx[i] * drdx[j] * (r * x[1] - a * x[0]) * pow((r2 + a2),-3.0);
            */
            d2ldx2[1][i][j] = (delta(j,1) * drdx[i] + drdx[j] * delta(i,1)) / (r2 + a2)
                            - 2.0 * drdx[j] * r * (x[1] * drdx[i] + r * delta(i,1) - a * delta(i,0)) / (r2 + a2) / (r2 + a2) + 
                            (-2.0 * drdx[j] * drdx[i] * (r * x[1] - a * x[0]) - 2.0 * r * d2rdx2[i][j] * (r * x[1] - a * x[0]) 
                            - 2.0 * r * drdx[i] * (drdx[j] * x[1] + r * delta(j,1) - a * delta(j,0))) / (r2 + a2) / (r2 + a2)
                            + 4.0 * r * drdx[i] * (r * x[1] - a * x[0]) * 2.0 * r * drdx[j] * pow((r2 + a2),-3.0);
            //el_z component
            d2ldx2[2][i][j] = -(delta(j,2) * drdx[i] + x[2] * d2rdx2[i][j] - 2.0 * x[2] * drdx[i] * drdx[j] /r + delta(i,2) * drdx[j]) / r2;
        }

        // then dltdi
        FOR1(i) { dltdx[i] = 0.0; }
        FOR2(i,j) { d2ltdx2[i][j] = 0.0; }        
    }

  public:
    // used to decide when to excise - ie when within the horizon of the BH
    // note that this is not templated over data_t
    double excise(const Cell<double> &current_cell) const
    {
        // black hole params - mass M and spin a
        const double M = m_params.mass;
        const double a = m_params.spin;
        const double a2 = a * a;

        // work out where we are on the grid
        const Coordinates<double> coords(current_cell, m_dx, m_params.center);
        const double x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const double r_plus = M + sqrt(M * M - a2);
        const double r_minus = M - sqrt(M * M - a2);

        const double outer_horizon =
            (x * x + y * y) / (2.0 * M * r_plus) + z * z / r_plus / r_plus;

        const double inner_horizon =
            (x * x + y * y) / (2.0 * M * r_minus) + z * z / r_minus / r_minus;

        // value less than 1 indicates we are within the horizon
        return sqrt(outer_horizon);
    }
};

#endif /* KERRSCHILDFIXEDBG_HPP_ */
