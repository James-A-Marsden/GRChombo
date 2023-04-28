/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(ISOTROPICBOOSTEDBH_HPP_)
#error "This file should only be included through IsotropicBoostedBH.hpp"
#endif

#ifndef ISOTROPICBOOSTEDBH_IMPL_HPP_
#define ISOTROPICBOOSTEDBH_IMPL_HPP_

// Compute the value of the initial vars on the grid
template <class data_t>
void IsotropicBoostedBH::compute(Cell<data_t> current_cell) const
{
    CCZ4Vars::VarsWithGauge<data_t> vars;
    VarsTools::assign(vars,
                      0.); // Set only the non-zero components explicitly below
    Coordinates<data_t> coords1(current_cell, m_dx, m_bh1_params.center);
    Coordinates<data_t> coords2(current_cell, m_dx, m_bh2_params.center);

    // work out where we are on the grid
    const double m1 = m_bh1_params.mass;
    const double v1 = m_bh1_params.momentum[0];
    const double gamma1 = 1. / sqrt(1. - (v1 * v1));
    const double m2 = m_bh2_params.mass;
    const double v2 = m_bh2_params.momentum[0];
    const double gamma2 = 1. / sqrt(1. - v2 * v2);

    ///////////////
    // 1st BH

    // boosts and coordinate objects
    double x = coords1.x * gamma1;
    double y = coords1.y; // set /tilde{t} to zero
    double z = coords1.z;         // coords.y+impact_parameter/2.;
    double r1 = sqrt(x * x + y * y + z * z);

    // first star physical variables
    double omega = 1. - m1 / (2. * r1);
    double psi = 1. + m1 / (2. * r1);
    double dx_psi = -m1 * x / (2. * pow(r1, 3));
    double dy_psi = -m1 * y / (2. * pow(r1, 3));

    double psi_omega = sqrt(pow(psi, 6) - v1 * v1 * omega * omega);
    double lapse_1 = omega * psi * psi / psi_omega / gamma1;
    double beta_x_1 =
        v1 * (pow(psi, 6) - pow(omega, 2)) / (psi_omega * psi_omega);
    double g_xx_1 = gamma1 * gamma1 * (psi_omega * psi_omega) / (psi * psi);
    double g_yy_1 = pow(psi, 4);
    double g_zz_1 = pow(psi, 4);
    double g_xx, g_yy, g_zz;

    double KLL_1[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double KLL_2[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double KLL[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double gammaLL[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double gammaUU[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double gammaUU_1[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double gammaUU_2[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double K1 = 0., K2 = 0.;
    gammaUU_1[0][0] = 1. / g_xx_1;
    gammaUU_1[1][1] = 1. / g_yy_1;
    gammaUU_1[2][2] = 1. / g_zz_1;

    KLL_1[0][0] =
        gamma1 * gamma1 * v1 * dx_psi *
        (2 * pow(psi, 6) * (psi + 2 * omega) - 2 * v1 * v1 * omega * omega) /
        pow(psi, 5) / psi_omega;
    KLL_1[1][1] = -2 * v1 * psi * omega * dx_psi / psi_omega;
    KLL_1[2][2] = KLL_1[1][1];
    KLL_1[0][1] = gamma1 * v1 * psi * dy_psi * (3 * omega + psi) / psi_omega;
    KLL_1[1][0] = KLL_1[0][1];
    KLL_1[2][0] = 0;
    KLL_1[0][2] = 0;
    KLL_1[2][1] = 0.;
    KLL_1[1][2] = 0.;
    //    FOR2(i,j) K1 += gammaUU_1[i][j]*KLL_1[i][j];

    /////////////////////
    // 2nd BH

    // boosts and coordinate objects
    x = coords2.x * gamma2;
    y = coords2.y; // set /tilde{t} to zero
    z = coords2.z;         // coords.y+impact_parameter/2.;
    double r2 = sqrt(x * x + y * y + z * z);

    // first star physical variables
    omega = 1. - m2 / (2. * r2);
    psi = 1. + m2 / (2. * r2);
    dx_psi = -m2 * x / (2. * pow(r2, 3));
    dy_psi = -m2 * y / (2. * pow(r2, 3));

    psi_omega = sqrt(pow(psi, 6) - v2 * v2 * omega * omega);
    double lapse_2 = omega * psi * psi / psi_omega / gamma2;
    double beta_x_2 =
        v2 * (pow(psi, 6) - pow(omega, 2)) / (psi_omega * psi_omega);
    double g_xx_2 = gamma2 * gamma2 * (psi_omega * psi_omega) / (psi * psi);
    double g_yy_2 = pow(psi, 4);
    double g_zz_2 = pow(psi, 4);

    gammaUU_2[0][0] = 1. / g_xx_2;
    gammaUU_2[1][1] = 1. / g_yy_2;
    gammaUU_2[2][2] = 1. / g_zz_2;

    KLL_2[0][0] =
        gamma2 * gamma2 * v2 * dx_psi *
        (2 * pow(psi, 6) * (psi + 2 * omega) - 2 * v2 * v2 * omega * omega) /
        pow(psi, 5) / psi_omega;
    KLL_2[1][1] = -2 * v2 * psi * omega * dx_psi / psi_omega;
    KLL_2[2][2] = KLL_2[1][1];
    KLL_2[0][1] = gamma2 * v2 * psi * dy_psi * (3 * omega + psi) / psi_omega;
    KLL_2[1][0] = KLL_2[0][1];
    KLL_2[2][0] = 0;
    KLL_2[0][2] = 0;
    KLL_2[2][1] = 0.;
    KLL_2[1][2] = 0.;

    ///////////////
    // Superpose solutions

    g_xx = g_xx_1 + g_xx_2 - 1.;
    g_yy = g_yy_1 + g_yy_2 - 1.;
    g_zz = g_zz_1 + g_zz_2 - 1.;
    gammaLL[0][0] = g_xx;
    gammaLL[1][1] = g_yy;
    gammaLL[2][2] = g_zz;
    gammaUU[0][0] = 1. / g_xx;
    gammaUU[1][1] = 1. / g_yy;
    gammaUU[2][2] = 1. / g_zz;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int m = 0; m < 3; m++)
            {
                for (int k = 0; k < 3; k++)
                {
                    KLL[i][j] += 0.5 * ((KLL_1[k][i] * gammaUU_1[k][m] +
                                         KLL_2[k][i] * gammaUU_2[k][m]) *
                                            gammaLL[j][m] +
                                        (KLL_1[k][j] * gammaUU_1[k][m] +
                                         KLL_2[k][j] * gammaUU_2[k][m]) *
                                            gammaLL[i][m]);
                }
            }
        }
    }

    double chi_ = pow(g_xx * g_yy * g_zz, -1. / 3.);

    vars.chi = chi_;
    vars.shift[0] = beta_x_1 + beta_x_2;
    vars.lapse = sqrt(vars.chi);
    //    vars.shift[0] = gammaUU[0][0]*(g_xx_1*beta_x_1+g_xx_2*beta_x_2);
    //    vars.lapse = 1/sqrt(lapse_1*lapse_1+lapse_2*lapse_2-1);

    double one_third = 1. / 3.;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            vars.K += KLL[i][j] * gammaUU[i][j];
        }
    }
    FOR2(i, j) vars.h[i][j] = vars.chi * gammaLL[i][j];
    FOR2(i, j)
    vars.A[i][j] = vars.chi * (KLL[i][j] - one_third * vars.K * gammaLL[i][j]);
    //vars.hww = vars.chi * gammaLL[2][2];
    //vars.Aww = vars.chi * (KLL[2][2] - one_third * vars.K * gammaLL[2][2]);

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* ISOTROPICBOOSTEDBH _IMPL_HPP_ */
