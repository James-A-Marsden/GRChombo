/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ADMISOKERRFIXEDBG_HPP_
#define ADMISOKERRFIXEDBG_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

//! Class which computes the initial conditions per arXiv 1401.1548
//! For a highly spinning BH in quasi isotropic coords
class ADMIsoKerrFixedBG
{
  public:
    //! Struct for the params of the  BH
    struct params_t
    {
        double mass = 1.0;                      //!<< The mass of the BH
        std::array<double, CH_SPACEDIM> center; //!< The center of the BH
        double spin = 0.5; //!< The spin 'a' in the z direction
    };

    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;
    const params_t m_params;

  protected:
    const double m_dx;

  public:
    ADMIsoKerrFixedBG(params_t a_params, double a_dx)
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

    /*
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
    }
    */
    /// Refer to Witek et al 1401.1548 for reference for
    /// Quasi Isotropic Kerr and variable conventions used here
    template <class data_t, template <typename> class vars_t>

    void compute_metric_background(vars_t<data_t> &vars,
                                   const Coordinates<double> &coords) const
    {
        // where am i?
        // const Coordinates<data_t> coords(current_cell, m_dx,
        // m_params.center);

        // black hole params - mass M and spin a
        const double M = m_params.mass;
        const double a = m_params.spin;
        const double a2 = a * a;

        // work out where we are on the grid
        // R is the quasi isotropic radial coord
        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const data_t x2 = x * x;
        const double y2 = y * y;
        const double z2 = z * z;
        const data_t R = coords.get_radius();
        const data_t R2 = R * R;
        // the radius in xy plane, subject to a floor
        const data_t rho2 = simd_max(x * x + y * y, 1e-8);
        const data_t rho = sqrt(rho2);

        // useful position quantities
        const data_t cos_theta = z / R;
        const data_t sin_theta = rho / R;

        const data_t cot_theta = cos_theta / sin_theta;
        const data_t cot_theta2 = cot_theta * cot_theta;

        const data_t cos_theta2 = cos_theta * cos_theta;
        const data_t sin_theta2 = sin_theta * sin_theta;
        const data_t cos_phi = x / rho;
        const data_t sin_phi = y / rho;

        const data_t cos_phi2 = cos_phi * cos_phi;
        const data_t sin_phi2 = sin_phi * sin_phi;

        const double r_plus = M + sqrt(M * M - a * a);
        const double r_minus = M - sqrt(M * M - a * a);

        // the Boyer Lindquist coord (eqn (46) 1401.1548)
        const data_t r_BL = R + 0.5 * r_plus + 0.0625 * r_plus * r_plus / R;
        const data_t r_BL2 = r_BL * r_BL;

        // find the quantities in eqns (45)
        const data_t Sigma = r_BL2 + a2 * cos_theta2;
        const data_t Delta = r_BL2 - 2.0 * M * r_BL + a2;
        // In the paper this is script 'A', not to be confused with A_ij
        const data_t AA = (r_BL2 + a2) * (r_BL2 + a2) - Delta * a2 * sin_theta2;

        const data_t Delta2 = Delta * Delta;
        const data_t Sigma2 = Sigma * Sigma;
        const data_t AA2 = AA * AA;

        // now the components per eqn (47) using psi04 = Sigma / R2;
        const data_t gamma_RR = Sigma * (1.0 + 0.25 * r_plus / R) *
                                (1.0 + 0.25 * r_plus / R) /
                                (R * (r_BL - r_minus));
        // gamma_pp = psi04 * AA * R2 / Sigma^2
        const data_t gamma_pp = AA / Sigma * sin_theta2;

        // Need to convert spherical metric to cartesian
        Tensor<2, data_t> spherical_g;
        FOR2(i, j) { spherical_g[i][j] = 0.0; }
        spherical_g[0][0] = gamma_RR;
        spherical_g[1][1] = Sigma; // gamma_tt = psi04 * R2 = Sigma
        spherical_g[2][2] = gamma_pp;

        // jacobian derivatives drdx, dthetadx etc
        using namespace TensorAlgebra;
        const Tensor<1, data_t> x_i = {x, y, z};
        const Tensor<1, data_t> drhodx = {x / rho, y / rho, 0.0};
        Tensor<1, data_t> dRdx = {x / R, y / R, z / R};

        // New derivatives of rho, R
        Tensor<2, data_t> d2rhodx2;
        FOR1(i)
        {
            d2rhodx2[0][i] = delta(i, 0) / rho - x * drhodx[i] / rho2;
            d2rhodx2[1][i] = delta(i, 1) / rho - y * drhodx[i] / rho2;
            d2rhodx2[2][i] = 0.0;
        }

        Tensor<2, data_t> d2Rdx2;

        FOR2(i, j) { d2Rdx2[i][j] = delta(i, j) / R - x_i[i] * dRdx[j] / R2; }

        Tensor<2, data_t> jac;
        Tensor<3, data_t> djacdx;
        jac[0][0] = x / R;
        jac[1][0] = cos_phi * z / R2;
        jac[2][0] = -y / rho2;
        jac[0][1] = y / R;
        jac[1][1] = sin_phi * z / R2;
        jac[2][1] = x / rho2;
        jac[0][2] = z / R;
        jac[1][2] = -rho / R2;
        jac[2][2] = 0.0;
        FOR2(i, j)
        {
            djacdx[0][i][j] = (delta(i, j) - x_i[i] * x_i[j] / R2) / R;
        }
        FOR1(i)
        {
            djacdx[1][0][i] =
                1.0 / R2 / rho * (delta(i, 0) * z + delta(i, 2) * x) +
                x * z / R2 / rho * (-drhodx[i] / rho - 2.0 * dRdx[i] / R);
            djacdx[1][1][i] =
                1.0 / R2 / rho * (delta(i, 1) * z + delta(i, 2) * y) +
                y * z / R2 / rho * (-drhodx[i] / rho - 2.0 * dRdx[i] / R);
            djacdx[1][2][i] = -drhodx[i] / R2 + 2.0 * rho * dRdx[i] / R2 / R;
            djacdx[2][0][i] = (-delta(i, 1) + 2.0 * y * drhodx[i] / rho) / rho2;
            djacdx[2][1][i] = (+delta(i, 0) - 2.0 * x * drhodx[i] / rho) / rho2;
            djacdx[2][2][i] = 0.0;
        }

        // Second order derivative of the Jacobian
        Tensor<4, data_t> d2jacdx2;

        FOR3(i, j, k)
        {
            d2jacdx2[0][i][j][k] = -delta(i, j) * dRdx[k] / R2 -
                                   delta(i, k) * x_i[j] / R2 / R -
                                   delta(j, k) * x_i[i] / R2 / R +
                                   3.0 * x_i[i] * x_i[j] * dRdx[k] / (R2 * R2);
        }

        FOR2(i, j)
        {
            d2jacdx2[1][0][i][j] =
                (-2.0 * dRdx[j] / R - drhodx[j] / rho) / R2 / rho *
                    (delta(i, 0) * z + delta(i, 2) * x -
                     x * z * drhodx[i] / rho - 2.0 * x * z * dRdx[i] / R) +
                1.0 / R2 / rho *
                    (delta(i, 0) * delta(j, 2) + delta(i, 2) * delta(j, 0) -
                     (delta(j, 0) * z + x * delta(j, 2)) * drhodx[i] / rho -
                     x * z * d2rhodx2[i][j] / rho +
                     x * z * drhodx[i] * drhodx[j] / rho2 -
                     2.0 * (delta(j, 0) * z + x * delta(j, 2)) * dRdx[i] / R -
                     2.0 * x * z * d2Rdx2[i][j] / R +
                     2.0 * x * z * dRdx[i] * dRdx[j] / R2);

            /*1.0 / R2 / rho * (delta(i,0) * delta(j,2) + delta(i,2) *
               delta(j,0))
                                    - (2.0 * dRdx[j]/R + drhodx[j]/rho) *
               (delta(i,0) * z + delta(i,2) * x -x * z * drhodx[i] / rho - 2.0 *
               dRdx[i] * x * z / R) / R2 / rho
                                    + (delta(j,0) * z + delta(j,2) * x) *
               (-drhodx[i] / rho - 2.0 * dRdx[i] / R)
                                    + x * z / R2 / rho * (-d2rhodx2[i][j] / rho
               - 2.0 * d2Rdx2[i][j] / R)
                                    + x * z /R2 / rho * (drhodx[i] * drhodx[j] /
               rho2 + 2.0 * dRdx[i] * dRdx[j] / R2)
                                    + (- 2.0 * x * z * dRdx[j] / R / R2 / rho -
               x * z * drhodx[j] / R2 / rho2) * (-drhodx[i] / rho - 2.0 *
               dRdx[i] / R);
                                    */

            d2jacdx2[1][1][i][j] =
                (-2.0 * dRdx[j] / R - drhodx[j] / rho) / R2 / rho *
                    (delta(i, 1) * z + delta(i, 2) * y -
                     y * z * drhodx[i] / rho - 2.0 * y * z * dRdx[i] / R) +
                1.0 / R2 / rho *
                    (delta(i, 1) * delta(j, 2) + delta(i, 2) * delta(j, 1) -
                     (delta(j, 1) * z + y * delta(j, 2)) * drhodx[i] / rho -
                     y * z * d2rhodx2[i][j] / rho +
                     y * z * drhodx[i] * drhodx[j] / rho2 -
                     2.0 * (delta(j, 1) * z + y * delta(j, 2)) * dRdx[i] / R -
                     2.0 * y * z * d2Rdx2[i][j] / R +
                     2.0 * y * z * dRdx[i] * dRdx[j] / R2);

            d2jacdx2[1][2][i][j] =
                (-d2rhodx2[i][j] + 2.0 * drhodx[i] * dRdx[j] / R +
                 2.0 * drhodx[j] * dRdx[i] / R + 2.0 * rho * d2Rdx2[i][j] / R -
                 6.0 * rho * dRdx[i] * dRdx[j] / R2) /
                R2;

            // djacdx[2][0][i] = (-delta(i, 1) + 2.0 * y * drhodx[i] / rho) /
            // rho2; djacdx[2][1][i] = (+delta(i, 0) - 2.0 * x * drhodx[i] /
            // rho) / rho2; djacdx[2][2][i] = 0.0;
            d2jacdx2[2][0][i][j] =
                2.0 *
                (delta(i, 1) * drhodx[j] + delta(j, 1) * drhodx[i] +
                 y * d2rhodx2[i][j] - 3.0 * y * drhodx[i] * drhodx[j] / rho) /
                rho2 / rho;

            d2jacdx2[2][1][i][j] =
                -2.0 *
                (delta(i, 0) * drhodx[j] + delta(j, 0) * drhodx[i] +
                 x * d2rhodx2[i][j] - 3.0 * x * drhodx[i] * drhodx[j] / rho) /
                rho2 / rho;

            d2jacdx2[2][2][i][j] = 0.0;
        }

        // derivs wrt R
        const data_t drBLdR = 1.0 - 0.0625 * r_plus * r_plus / R2;
        const data_t dDeltadR = drBLdR * (2.0 * r_BL - 2.0 * M);
        const data_t dSigmadR = 2.0 * r_BL * drBLdR;
        const data_t dAAdR =
            4.0 * drBLdR * r_BL * (r_BL2 + a2) - dDeltadR * a2 * sin_theta2;
        const data_t dgammaRRdR =
            gamma_RR * (dSigmadR / Sigma - 1.0 / R - drBLdR / (r_BL - r_minus) -
                        0.5 * r_plus / R2 / (1.0 + 0.25 * r_plus / R));
        const data_t dgammappdR = gamma_pp * (dAAdR / AA - dSigmadR / Sigma);

        // derivs wrt theta
        const data_t dSigmadtheta = -2.0 * a2 * cos_theta * sin_theta;
        const data_t dAAdtheta = -2.0 * Delta * a2 * sin_theta * cos_theta;
        const data_t dgammaRRdtheta = gamma_RR * dSigmadtheta / Sigma;
        const data_t dgammappdtheta =
            gamma_pp * (dAAdtheta / AA - dSigmadtheta / Sigma) +
            2.0 * cos_theta * sin_theta * AA / Sigma;

        // Calculate the gradients needed (wrt x, y, z)
        Tensor<1, data_t> drBLdx;
        Tensor<1, data_t> dSigmadx;
        Tensor<1, data_t> dDeltadx;
        Tensor<1, data_t> dAAdx;
        Tensor<1, data_t> dgammaRRdx;
        Tensor<1, data_t> dgammappdx;
        FOR1(i)
        {
            drBLdx[i] = jac[0][i] * drBLdR;
            dDeltadx[i] = jac[0][i] * dDeltadR;
            dSigmadx[i] = jac[0][i] * dSigmadR + jac[1][i] * dSigmadtheta;
            dAAdx[i] = jac[0][i] * dAAdR + jac[1][i] * dAAdtheta;
            dgammaRRdx[i] = jac[0][i] * dgammaRRdR + jac[1][i] * dgammaRRdtheta;
            dgammappdx[i] = jac[0][i] * dgammappdR + jac[1][i] * dgammappdtheta;
        }
        // const data_t gamma_pp = AA / Sigma * sin_theta2;
        // derivs wrt x (a bit janky, but should work) Note this is with respect
        // to x,y,z - not R, theta, phi

        Tensor<1, data_t> dsin_thetadx;
        FOR1(i) { dsin_thetadx[i] = drhodx[i] / R - rho * dRdx[i] / R2; }
        Tensor<1, data_t> dcos_thetadx;
        FOR1(i) { dcos_thetadx[i] = delta(i, 2) / R - z * dRdx[i] / R2; }
        Tensor<1, data_t> dcot_thetadx;
        FOR1(i) { dcot_thetadx[i] = delta(i, 2) / rho - drhodx[i] * z / rho2; }
        Tensor<1, data_t> dsin_phidx;
        FOR1(i) { dsin_phidx[i] = delta(i, 1) / rho - y * drhodx[i] / rho2; }
        Tensor<1, data_t> dcos_phidx;
        FOR1(i) { dcos_phidx[i] = delta(i, 0) / rho - x * drhodx[i] / rho2; }
        Tensor<1, data_t> drBLdR_ddx;
        FOR1(i) { drBLdR_ddx[i] = 0.125 * r_plus * r_plus * dRdx[i] / R2 / R; }

        Tensor<1, data_t> dDeltadR_ddx;
        FOR1(i)
        {
            dDeltadR_ddx[i] = drBLdR_ddx[i] * (2.0 * r_BL - 2.0 * M) +
                              drBLdR * (2.0 * drBLdx[i]);
        }
        Tensor<1, data_t> dSigmadR_ddx;
        FOR1(i)
        {
            dSigmadR_ddx[i] =
                2.0 * drBLdx[i] * drBLdR + 2.0 * r_BL * drBLdR_ddx[i];
        }
        Tensor<1, data_t> dSigmadtheta_ddx;
        FOR1(i)
        {
            dSigmadtheta_ddx[i] =
                -2.0 * a2 *
                (dsin_thetadx[i] * cos_theta + dcos_thetadx[i] * sin_theta);
        }
        Tensor<1, data_t> dAAdtheta_ddx;
        FOR1(i)
        {
            dAAdtheta_ddx[i] = -2.0 * a2 * Delta *
                                   (dsin_thetadx[i] * cos_theta +
                                    dcos_thetadx[i] * sin_theta) -
                               2.0 * dDeltadx[i] * a2 * sin_theta * cos_theta;
        }
        Tensor<1, data_t> dAAdR_ddx;
        FOR1(i)
        {
            dAAdR_ddx[i] =
                4.0 * (drBLdR_ddx[i] * r_BL2 * r_BL +
                       3.0 * drBLdR * r_BL2 * drBLdx[i] +
                       drBLdR_ddx[i] * r_BL * a2 + drBLdR * drBLdx[i] * a2) -
                dDeltadR_ddx[i] * a2 * sin_theta2 -
                dDeltadR * a2 * 2.0 * sin_theta * dsin_thetadx[i];
        }
        Tensor<1, data_t> dgammaRRdR_ddx;
        FOR1(i)
        {
            dgammaRRdR_ddx[i] =
                dgammaRRdx[i] * dgammaRRdR / gamma_RR +
                gamma_RR *
                    (dSigmadR_ddx[i] / Sigma - dSigmadx[i] * dSigmadR / Sigma2 +
                     dRdx[i] / R2 - drBLdR_ddx[i] / (r_BL - r_minus) +
                     drBLdR * drBLdx[i] / (r_BL - r_minus) / (r_BL - r_minus) +
                     r_plus * dRdx[i] / R2 / R / (1.0 + 0.25 * r_plus / R) -
                     0.125 * r_plus * r_plus / R2 / R2 * dRdx[i] /
                         (1.0 + 0.25 * r_plus / R) / (1.0 + 0.25 * r_plus / R));
        }
        Tensor<1, data_t> dgammaRRdtheta_ddx;
        FOR1(i)
        {
            dgammaRRdtheta_ddx[i] =
                dgammaRRdx[i] * dSigmadtheta / Sigma +
                gamma_RR * dSigmadtheta_ddx[i] / Sigma -
                gamma_RR * dSigmadtheta * dSigmadx[i] / Sigma2;
        }
        Tensor<1, data_t> dgammappdR_ddx;
        FOR1(i)
        {
            dgammappdR_ddx[i] =
                dgammappdx[i] * (dAAdR / AA - dSigmadR / Sigma) +
                gamma_pp *
                    (dAAdR_ddx[i] / AA - dAAdR * dAAdx[i] / AA2 -
                     dSigmadR_ddx[i] / Sigma + dSigmadR * dSigmadx[i] / Sigma2);
        }
        Tensor<1, data_t> dgammappdtheta_ddx;
        FOR1(i)
        {
            dgammappdtheta_ddx[i] =
                dgammappdx[i] * (dAAdtheta / AA - dSigmadtheta / Sigma) +
                gamma_pp * (dAAdtheta_ddx[i] / AA - dAAdtheta * dAAdx[i] / AA2 -
                            dSigmadtheta_ddx[i] / Sigma +
                            dSigmadtheta * dSigmadx[i] / Sigma2) +
                2.0 *
                    (dcos_thetadx[i] * sin_theta +
                     cos_theta * dsin_thetadx[i]) *
                    AA / Sigma +
                2.0 * cos_theta * sin_theta *
                    (dAAdx[i] / Sigma - AA * dSigmadx[i] / Sigma2);
        }
        /*
                const data_t dgammaRRdR =
            gamma_RR * (dSigmadR / Sigma - 1.0 / R - drBLdR / (r_BL - r_minus) -
                        0.5 * r_plus / R2 / (1.0 + 0.25 * r_plus / R));
        const data_t dgammappdR = gamma_pp * (dAAdR / AA - dSigmadR / Sigma);

        // derivs wrt theta
        const data_t dSigmadtheta = -2.0 * a2 * cos_theta * sin_theta;
        const data_t dAAdtheta = -2.0 * Delta * a2 * sin_theta * cos_theta;
        const data_t dgammaRRdtheta = gamma_RR * dSigmadtheta / Sigma;
        const data_t dgammappdtheta =
            gamma_pp * (dAAdtheta / AA - dSigmadtheta / Sigma) +
            2.0 * cos_theta * sin_theta * AA / Sigma;
        */

        // const data_t dgammappdtheta =
        // gamma_pp * (dAAdtheta / AA - dSigmadtheta / Sigma) +
        // 2.0 * cos_theta * sin_theta * AA / Sigma;
        Tensor<2, data_t> d2rBLdx2;
        Tensor<2, data_t> d2Sigmadx2;
        Tensor<2, data_t> d2Deltadx2;
        Tensor<2, data_t> d2AAdx2;
        Tensor<2, data_t> d2gammaRRdx2;
        Tensor<2, data_t> d2gammappdx2;
        FOR2(i, j)
        {

            // dDeltadx[i] = jac[0][i] * dDeltadR;
            // dSigmadx[i] = jac[0][i] * dSigmadR + jac[1][i] * dSigmadtheta;
            // dAAdx[i] = jac[0][i] * dAAdR + jac[1][i] * dAAdtheta;
            // dgammaRRdx[i] = jac[0][i] * dgammaRRdR + jac[1][i] *
            // dgammaRRdtheta; dgammappdx[i] = jac[0][i] * dgammappdR +
            // jac[1][i] * dgammappdtheta;
            d2rBLdx2[i][j] =
                djacdx[0][i][j] * drBLdR + jac[0][i] * drBLdR_ddx[j];
            d2Sigmadx2[i][j] = djacdx[0][i][j] * dSigmadR +
                               jac[0][i] * dSigmadR_ddx[j] +
                               djacdx[1][i][j] * dSigmadtheta +
                               jac[1][i] * dSigmadtheta_ddx[j];
            d2Deltadx2[i][j] =
                djacdx[0][i][j] * dDeltadR + jac[0][i] * dDeltadR_ddx[j];
            d2AAdx2[i][j] = djacdx[0][i][j] * dAAdR + jac[0][i] * dAAdR_ddx[j] +
                            djacdx[1][i][j] * dAAdtheta +
                            jac[1][i] * dAAdtheta_ddx[j];
            d2gammaRRdx2[i][j] = djacdx[0][i][j] * dgammaRRdR +
                                 jac[0][i] * dgammaRRdR_ddx[j] +
                                 djacdx[1][i][j] * dgammaRRdtheta +
                                 jac[1][i] * dgammaRRdtheta_ddx[j];
            d2gammappdx2[i][j] = djacdx[0][i][j] * dgammappdR +
                                 jac[0][i] * dgammappdR_ddx[j] +
                                 djacdx[1][i][j] * dgammappdtheta +
                                 jac[1][i] * dgammappdtheta_ddx[j];
        }
        // dgammappdx[i] = jac[0][i] * dgammappdR + jac[1][i] * dgammappdtheta;

        // populate ADM vars - lapse and shift
        // use analytic continuation for lapse within horizon
        data_t sign_lapse = (R - 0.25 * r_plus) / abs(R - 0.25 * r_plus);
        vars.lapse = sign_lapse * sqrt(Delta * Sigma / AA);

        // now the shift
        const data_t beta_phi = -2.0 * M * a * r_BL / AA;
        FOR1(i) { vars.shift[i] = 0.0; }
        vars.shift[0] = -y * beta_phi;
        vars.shift[1] = x * beta_phi;

        // spatial metric - convert spherical to cartesian
        FOR2(i, j)
        {
            vars.gamma[i][j] = 0.0;

            FOR2(k, m)
            {
                vars.gamma[i][j] += spherical_g[k][m] * jac[k][i] * jac[m][j];
            }
        }
        const auto gamma_UU = TensorAlgebra::compute_inverse_sym(vars.gamma);
        vars.gamma_UU = TensorAlgebra::compute_inverse_sym(vars.gamma);

        // Calculate partial derivative of spatial metric
        FOR3(i, j, k)
        {
            vars.d1_gamma[i][j][k] = djacdx[0][j][k] * jac[0][i] * gamma_RR +
                                     djacdx[0][i][k] * jac[0][j] * gamma_RR +
                                     djacdx[1][j][k] * jac[1][i] * Sigma +
                                     djacdx[1][i][k] * jac[1][j] * Sigma +
                                     djacdx[2][j][k] * jac[2][i] * gamma_pp +
                                     djacdx[2][i][k] * jac[2][j] * gamma_pp +
                                     jac[0][i] * jac[0][j] * dgammaRRdx[k] +
                                     jac[1][i] * jac[1][j] * dSigmadx[k] +
                                     jac[2][i] * jac[2][j] * dgammappdx[k];
        }

        // NEW - second derivative of the spatial metric

        FOR3(i, j, k)
        {
            FOR1(l)
            {
                vars.d2_gamma[i][j][k][l] =
                    d2jacdx2[0][j][k][l] * jac[0][i] * gamma_RR +
                    djacdx[0][j][k] * djacdx[0][i][l] * gamma_RR +
                    djacdx[0][j][k] * jac[0][i] * dgammaRRdx[l] +
                    d2jacdx2[0][i][k][l] * jac[0][j] * gamma_RR +
                    djacdx[0][i][k] * djacdx[0][j][l] * gamma_RR +
                    djacdx[0][i][k] * jac[0][j] * dgammaRRdx[l] +
                    d2jacdx2[1][j][k][l] * jac[1][i] * Sigma +
                    djacdx[1][j][k] * djacdx[1][i][l] * Sigma +
                    djacdx[1][j][k] * jac[1][i] * dSigmadx[l] +
                    d2jacdx2[1][i][k][l] * jac[1][j] * Sigma +
                    djacdx[1][i][k] * djacdx[1][j][l] * Sigma +
                    djacdx[1][i][k] * jac[1][j] * dSigmadx[l] +
                    d2jacdx2[2][j][k][l] * jac[2][i] * gamma_pp +
                    djacdx[2][j][k] * djacdx[2][i][l] * gamma_pp +
                    djacdx[2][j][k] * jac[2][i] * dgammappdx[l] +
                    d2jacdx2[2][i][k][l] * jac[2][j] * gamma_pp +
                    djacdx[2][i][k] * djacdx[2][j][l] * gamma_pp +
                    djacdx[2][i][k] * jac[2][j] * dgammappdx[l] +
                    djacdx[0][i][l] * jac[0][j] * dgammaRRdx[k] +
                    jac[0][i] * djacdx[0][j][l] * dgammaRRdx[k] +
                    jac[0][i] * jac[0][j] * d2gammaRRdx2[k][l] +
                    djacdx[1][i][l] * jac[1][j] * dSigmadx[k] +
                    jac[1][i] * djacdx[1][j][l] * dSigmadx[k] +
                    jac[1][i] * jac[1][j] * d2Sigmadx2[k][l] +
                    djacdx[2][i][l] * jac[2][j] * dgammappdx[k] +
                    jac[2][i] * djacdx[2][j][l] * dgammappdx[k] +
                    jac[2][i] * jac[2][j] * d2gammappdx2[k][l];
            }
        }

        FOR1(i)
        {
            vars.d1_gamma_UU[0][0][i] =
                -2.0 * drhodx[i] / rho * gamma_UU[0][0] +
                1.0 / rho2 *
                    (2.0 * (y * delta(i, 1) * rho2 + rho * drhodx[i] * y2) /
                         gamma_pp -
                     y2 * rho2 * dgammappdx[i] / gamma_pp / gamma_pp +
                     2.0 * x * delta(i, 0) / R2 / R2 / Sigma / gamma_RR *
                         (R2 * rho2 * Sigma + R2 * R2 * z2 * gamma_RR) +
                     x2 / R2 / R2 / Sigma / gamma_RR *
                         (2.0 * R * dRdx[i] * rho2 * Sigma +
                          2.0 * R2 * rho * drhodx[i] * Sigma +
                          R2 * rho2 * dSigmadx[i] +
                          4.0 * R2 * R * dRdx[i] * z2 * gamma_RR +
                          2.0 * R2 * R2 * z * delta(i, 2) * gamma_RR +
                          R2 * R2 * z2 * dgammaRRdx[i]) +
                     x2 * (R2 * rho2 * Sigma + R2 * R2 * z2 * gamma_RR) *
                         (-4.0 * dRdx[i] / R2 / R2 / R / Sigma / gamma_RR -
                          dSigmadx[i] / R2 / R2 / Sigma2 / gamma_RR -
                          dgammaRRdx[i] / R2 / R2 / Sigma / gamma_RR /
                              gamma_RR));

            vars.d1_gamma_UU[0][1][i] =
                (delta(i, 0) / x + delta(i, 1) / y - 2.0 * drhodx[i] / rho) *
                    gamma_UU[0][1] +
                x * y / rho2 *
                    (-2.0 * drhodx[i] * rho / gamma_pp +
                     rho2 * dgammappdx[i] / gamma_pp / gamma_pp +
                     (2.0 * dRdx[i] * R * rho2 * Sigma +
                      2.0 * R2 * drhodx[i] * rho * Sigma +
                      R2 * rho2 * dSigmadx[i] +
                      4.0 * R2 * R * dRdx[i] * z2 * gamma_RR +
                      2.0 * R2 * R2 * z * delta(i, 2) * gamma_RR +
                      R2 * R2 * z2 * dgammaRRdx[i]) /
                         R2 / R2 / Sigma / gamma_RR +
                     (-4.0 * dRdx[i] / R2 / R2 / R / Sigma / gamma_RR -
                      dSigmadx[i] / R2 / R2 / Sigma2 / gamma_RR -
                      dgammaRRdx[i] / R2 / R2 / Sigma / gamma_RR / gamma_RR) *
                         (R2 * rho2 * Sigma + R2 * R2 * z2 * gamma_RR));

            vars.d1_gamma_UU[0][2][i] =
                (delta(i, 0) / x + delta(i, 2) / z - 2.0 * dRdx[i] / R -
                 dSigmadx[i] / Sigma - dgammaRRdx[i] / gamma_RR) *
                    gamma_UU[0][2] -
                x * z / R2 / Sigma / gamma_RR *
                    (-dSigmadx[i] + 2.0 * R * dRdx[i] * gamma_RR +
                     R2 * dgammaRRdx[i]);

            vars.d1_gamma_UU[1][1][i] =
                -2.0 * drhodx[i] / rho * gamma_UU[1][1] +
                1.0 / rho2 *
                    (2.0 * (x * delta(i, 0) * rho2 + rho * drhodx[i] * x2) /
                         gamma_pp -
                     x2 * rho2 * dgammappdx[i] / gamma_pp / gamma_pp +
                     2.0 * y * delta(i, 1) / R2 / R2 / Sigma / gamma_RR *
                         (R2 * rho2 * Sigma + R2 * R2 * z2 * gamma_RR) +
                     y2 / R2 / R2 / Sigma / gamma_RR *
                         (2.0 * R * dRdx[i] * rho2 * Sigma +
                          2.0 * R2 * rho * drhodx[i] * Sigma +
                          R2 * rho2 * dSigmadx[i] +
                          4.0 * R2 * R * dRdx[i] * z2 * gamma_RR +
                          2.0 * R2 * R2 * z * delta(i, 2) * gamma_RR +
                          R2 * R2 * z2 * dgammaRRdx[i]) +
                     y2 * (R2 * rho2 * Sigma + R2 * R2 * z2 * gamma_RR) *
                         (-4.0 * dRdx[i] / R2 / R2 / R / Sigma / gamma_RR -
                          dSigmadx[i] / R2 / R2 / Sigma2 / gamma_RR -
                          dgammaRRdx[i] / R2 / R2 / Sigma / gamma_RR /
                              gamma_RR));

            vars.d1_gamma_UU[1][2][i] =
                gamma_UU[1][2] * (-2.0 * dRdx[i] / R - dSigmadx[i] / Sigma -
                                  dgammaRRdx[i] / gamma_RR + delta(i, 1) / y +
                                  delta(i, 2) / z - 2.0 * drhodx[i] / rho +
                                  2.0 * drhodx[i] / rho) -
                y * z / R2 / Sigma / gamma_RR *
                    (-dSigmadx[i] + 2.0 * dRdx[i] * R * gamma_RR +
                     R2 * dgammaRRdx[i]);

            vars.d1_gamma_UU[2][2][i] =
                gamma_UU[2][2] * (-2.0 * dRdx[i] / R - dSigmadx[i] / Sigma -
                                  dgammaRRdx[i] / gamma_RR) +
                1.0 / R2 / Sigma / gamma_RR *
                    (2.0 * z * delta(i, 2) * Sigma + z2 * dSigmadx[i] +
                     2.0 * R * dRdx[i] * rho2 * gamma_RR +
                     2.0 * R2 * rho * drhodx[i] * gamma_RR +
                     R2 * rho2 * dgammaRRdx[i]);

            vars.d1_gamma_UU[1][0][i] = vars.d1_gamma_UU[0][1][i];
            vars.d1_gamma_UU[2][0][i] = vars.d1_gamma_UU[0][2][i];
            vars.d1_gamma_UU[2][1][i] = vars.d1_gamma_UU[1][2][i];
        }

        // calculate derivs of lapse and shift
        // use analytic continuation of lapse within the horizon
        // (taken care of automatically by use of vars.lapse)
        FOR1(i)
        {
            vars.d1_lapse[i] =
                0.5 * vars.lapse *
                (dDeltadx[i] / Delta + dSigmadx[i] / Sigma - dAAdx[i] / AA);
        }

        // Extra derivative of the lapse

        FOR2(i, j)
        {
            vars.d2_lapse[i][j] =
                vars.d1_lapse[j] * vars.d1_lapse[i] / vars.lapse +
                0.5 * vars.lapse *
                    (d2Deltadx2[i][j] / Delta -
                     dDeltadx[i] * dDeltadx[j] / Delta2 +
                     d2Sigmadx2[i][j] / Sigma -
                     dSigmadx[i] * dSigmadx[j] / Sigma2 - d2AAdx2[i][j] / AA +
                     dAAdx[i] * dAAdx[j] / AA2);
        }

        FOR2(i, j) { vars.d1_shift[i][j] = 0.0; }

        FOR1(i)
        {
            vars.d1_shift[0][i] =
                vars.shift[0] *
                (drBLdx[i] / r_BL - dAAdx[i] / AA + delta(i, 1) / y);
            vars.d1_shift[1][i] =
                vars.shift[1] *
                (drBLdx[i] / r_BL - dAAdx[i] / AA + delta(i, 0) / x);
        }

        FOR3(i, j, k) { vars.d2_shift[i][j][k] = 0.0; }

        FOR2(i, j)
        {
            vars.d2_shift[0][i][j] =
                vars.d1_shift[0][j] * vars.d1_shift[0][i] / vars.shift[0] +
                vars.shift[0] *
                    (d2rBLdx2[i][j] / r_BL - drBLdx[i] * drBLdx[j] / r_BL2 -
                     d2AAdx2[i][j] / AA + dAAdx[i] * dAAdx[j] / AA2 -
                     delta(i, 1) * delta(j, 1) / (y * y));

            vars.d2_shift[1][i][j] =
                vars.d1_shift[1][j] * vars.d1_shift[1][i] / vars.shift[1] +
                vars.shift[1] *
                    (d2rBLdx2[i][j] / r_BL - drBLdx[i] * drBLdx[j] / r_BL2 -
                     d2AAdx2[i][j] / AA + dAAdx[i] * dAAdx[j] / AA2 -
                     delta(i, 0) * delta(j, 0) / (x * x));
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
        vars.K = compute_trace(vars.K_tensor, gamma_UU);

        FOR2(i, j)
        {
            FOR2(k, m)
            {
                vars.d1_chris_phys[i][j][k][m] = 0.0;

                FOR1(n)
                {
                    vars.d1_chris_phys[i][j][k][m] +=
                        0.5 * vars.d1_gamma_UU[i][n][m] *
                            (vars.d1_gamma[k][n][j] + vars.d1_gamma[n][j][k] -
                             vars.d1_gamma[j][k][n])

                        + 0.5 * gamma_UU[i][n] *
                              (vars.d2_gamma[k][n][j][m] +
                               vars.d2_gamma[n][j][k][m] -
                               vars.d2_gamma[j][k][n][m]);
                }
            }
        }
        FOR3(i, j, k)
        {
            vars.d1_K_tensor[i][j][k] = 0.0;

            FOR1(m)
            {
                vars.d1_K_tensor[i][j][k] +=
                    vars.d1_gamma[m][j][k] * vars.d1_shift[m][i] +
                    vars.gamma[m][j] * vars.d2_shift[m][i][k] +
                    vars.d1_gamma[m][i][k] * vars.d1_shift[m][j] +
                    vars.gamma[m][i] * vars.d2_shift[m][j][k] +
                    (vars.d2_gamma[m][i][j][k] + vars.d2_gamma[m][j][i][k]) *
                        vars.shift[m] +
                    (vars.d1_gamma[m][i][j] + vars.d1_gamma[m][j][i]) *
                        vars.d1_shift[m][k];
                FOR1(n)
                {
                    vars.d1_K_tensor[i][j][k] +=
                        -2.0 * (vars.d1_chris_phys[m][i][j][k] *
                                    vars.gamma[m][n] * vars.shift[n] +
                                chris_phys.ULL[m][i][j] *
                                    vars.d1_gamma[m][n][k] * vars.shift[n] +
                                chris_phys.ULL[m][i][j] * vars.gamma[m][n] *
                                    vars.d1_shift[n][k]);
                }
            }
            vars.d1_K_tensor[i][j][k] *= 0.5 / vars.lapse;
            vars.d1_K_tensor[i][j][k] +=
                -vars.d1_lapse[k] / vars.lapse * vars.K_tensor[i][j];

            // vars.d1_K_tensor[i][j][k] = 0.0;
        }
        // Derivative of the trace, \partial_i K = \partial_i(gamma^jk K_jk)
        FOR1(i)
        {
            vars.d1_K[i] = 0.0;
            FOR2(j, k)
            {
                vars.d1_K[i] +=
                    vars.d1_gamma_UU[j][k][i] * vars.K_tensor[j][k] +
                    gamma_UU[j][k] * vars.d1_K_tensor[j][k][i];
            }
        }

        // spatial riemann curvature tensor

        FOR1(i)
        {
            FOR3(j, k, l)
            {
                vars.riemann_phys_ULLL[i][j][k][l] =
                    vars.d1_chris_phys[i][l][j][k] -
                    vars.d1_chris_phys[i][k][j][l];

                FOR1(m)
                {
                    vars.riemann_phys_ULLL[i][j][k][l] +=
                        chris_phys.ULL[m][l][j] * chris_phys.ULL[i][m][k] -
                        chris_phys.ULL[m][k][j] * chris_phys.ULL[i][m][l];
                    // vars.riemann_phys_ULLL[i][j][k][l] = 0;
                }
            }
        }

        // spatial ricci tensor
        FOR2(i, j)
        {
            vars.ricci_phys[i][j] = 0;
            FOR1(k)
            {
                vars.ricci_phys[i][j] += vars.riemann_phys_ULLL[k][i][k][j];
            }
        }
    }

    // Other 3+1 ADM derivatives

  public:
    // used to decide when to excise - ie when within the horizon of the BH
    // note that this is not templated over data_t
    double excise(const Coordinates<double> &coords) const
    {
        // where am i?
        // const Coordinates<double> coords(current_cell, m_dx,
        // m_params.center);

        // black hole params - mass M and boost v
        // "boost" is the gamma factor for the boost
        const double M = m_params.mass;
        const double a = m_params.spin;

        // the quasi isotropic Kerr radius
        const double R = coords.get_radius();

        // compare this to horizon in quasi isotropic Kerr coords
        // which is r+/4
        const double r_horizon = 0.25 * (M + sqrt(M * M - a * a));

        return R / r_horizon;
    }
};

#endif /* ADMISOKERRFIXEDBG_HPP_ */