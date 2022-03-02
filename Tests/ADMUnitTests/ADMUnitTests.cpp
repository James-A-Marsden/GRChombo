/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Chombo includes
#include "BoxIterator.H"
#include "FArrayBox.H"

// Other includes
#include <iostream>

// Our includes
#include "BoxLoops.hpp"
#include "ADMTestsCompute.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SixthOrderDerivatives.hpp"
#include "UserVariables.hpp"
#include "Coordinates.hpp"
#include "ADMKerrSchildFixedBG.hpp"

// Chombo namespace
#include "UsingNamespace.H"

bool is_wrong(double value, double correct_value, std::string deriv_type)
{
    if (abs(value - correct_value) > 1e-10)
    {
        std::cout.precision(17);
        std::cout << "Test of " << deriv_type << " failed "
                  << " with value " << value << " instead of " << correct_value
                  << ".\n";
        return true;
    }
    else
    {
        std::cout << deriv_type << " " << abs(value - correct_value) << endl; 
        return false;
    }
}

template <class data_t> struct LocalVars
{
    Tensor<2, data_t> gamma;
    Tensor<2, data_t> K_tensor;
    data_t K;
    data_t lapse;
    Tensor<1, data_t> shift;
    Tensor<2, Tensor<1, data_t>> d1_gamma;
    Tensor<1, data_t> d1_lapse;
    Tensor<2, data_t> d1_shift;

    //Extra variables needed for tensor field case
    Tensor<2, data_t> gamma_UU;
    Tensor<1, data_t> d1_K;    
    Tensor<3, data_t> d1_K_tensor; 
    Tensor<2, data_t> d2_lapse;
    Tensor<3, data_t> d2_shift;
    Tensor<4, data_t> d2_gamma;
    Tensor<3, data_t> d1_gamma_UU;
    Tensor<4, data_t> d1_chris_phys; 
    Tensor<4, data_t> riemann_phys_ULLL;
    Tensor<2, data_t> ricci_phys;
};

int main()
{
    //const int num_cells = 512;
    const int num_cells = 100;
    // box is flat in y direction to make test cheaper
    IntVect domain_hi_vect(num_cells - 1, 0, num_cells - 1);
    Box box(IntVect(0, 0, 0), domain_hi_vect);
    Box ghosted_box(IntVect(-4, -4, -4),
                    IntVect(num_cells + 3, 4, num_cells + 3));

    FArrayBox in_fab(ghosted_box, NUM_VARS);
    FArrayBox out_fab(box, NUM_VARS);


    ADMKerrSchildFixedBG::params_t bg_params;
    bg_params.center[0] = 2.0;
    bg_params.center[1] = -5.0;
    bg_params.center[2] = 3.0;
    bg_params.mass = 1.0;
    bg_params.spin = 0.0;

    const double dx = 10.0 / num_cells;
    LocalVars<double> local_vars;

    ADMKerrSchildFixedBG kerr_bh(bg_params, dx);

    /*
    std::cout << "bh mass = " << bg_params.mass << "\n";
    std::cout << "BH centre \n";
    std::cout << bg_params.center[0] << " " << bg_params.center[1] << " " << bg_params.center[2] << "\n";
    */
    BoxIterator bit_ghost(ghosted_box);
    for (bit_ghost.begin(); bit_ghost.ok(); ++bit_ghost)
    {
        //const double x = (0.5 + bit_ghost()[0]) * dx;
        //const double z = (0.5 + bit_ghost()[2]) * dx;
        
        const Coordinates<double> ghost_coords(IntVect(bit_ghost()[0], bit_ghost()[1], bit_ghost()[2]), dx, bg_params.center);

        kerr_bh.compute_metric_background(local_vars, ghost_coords);

        in_fab(bit_ghost(), c_d1_gamma_UU11) = local_vars.gamma_UU[0][0];
        in_fab(bit_ghost(), c_d1_gamma_UU12) = local_vars.gamma_UU[0][1];
        in_fab(bit_ghost(), c_d1_gamma_UU13) = local_vars.gamma_UU[0][2];
        in_fab(bit_ghost(), c_d1_gamma_UU22) = local_vars.gamma_UU[1][1];
        in_fab(bit_ghost(), c_d1_gamma_UU23) = local_vars.gamma_UU[1][2];
        in_fab(bit_ghost(), c_d1_gamma_UU33) = local_vars.gamma_UU[2][2];

        in_fab(bit_ghost(), c_d2_gamma11) = local_vars.gamma[0][0];
        in_fab(bit_ghost(), c_d2_gamma12) = local_vars.gamma[0][1];
        in_fab(bit_ghost(), c_d2_gamma13) = local_vars.gamma[0][2];
        in_fab(bit_ghost(), c_d2_gamma22) = local_vars.gamma[1][1];
        in_fab(bit_ghost(), c_d2_gamma23) = local_vars.gamma[1][2];
        in_fab(bit_ghost(), c_d2_gamma33) = local_vars.gamma[2][2];

        in_fab(bit_ghost(), c_d1_K_tensor11) = local_vars.K_tensor[0][0];
        in_fab(bit_ghost(), c_d1_K_tensor12) = local_vars.K_tensor[0][1];
        in_fab(bit_ghost(), c_d1_K_tensor13) = local_vars.K_tensor[0][2];
        in_fab(bit_ghost(), c_d1_K_tensor22) = local_vars.K_tensor[1][1];
        in_fab(bit_ghost(), c_d1_K_tensor23) = local_vars.K_tensor[1][2];
        in_fab(bit_ghost(), c_d1_K_tensor33) = local_vars.K_tensor[2][2];

        in_fab(bit_ghost(), c_d2_lapse) = local_vars.lapse;

        in_fab(bit_ghost(), c_d1_shift1) = local_vars.shift[0];
        in_fab(bit_ghost(), c_d1_shift2) = local_vars.shift[1];
        in_fab(bit_ghost(), c_d1_shift3) = local_vars.shift[2];

        in_fab(bit_ghost(), c_d2_shift1) = local_vars.shift[0];
        in_fab(bit_ghost(), c_d2_shift2) = local_vars.shift[1];
        in_fab(bit_ghost(), c_d2_shift3) = local_vars.shift[2];
    }

    // Fourth order derivatives
    BoxLoops::loop(ADMTestsCompute<SixthOrderDerivatives>(dx), in_fab,
                   out_fab);

    BoxIterator bit(box);
    for (bit.begin(); bit.ok(); ++bit)
    {
        const double x = (0.5 + bit()[0]) * dx;
        const double z = (0.5 + bit()[2]) * dx;

        const Coordinates<double> coords(IntVect(bit()[0], bit()[1], bit()[2]), dx, bg_params.center);
        kerr_bh.compute_metric_background(local_vars, coords);
   
        bool error = false;
        
        error |= is_wrong(out_fab(bit(), c_d1_gamma_UU11), local_vars.d1_gamma_UU[0][0][0], "c_d1_gammaUU11");
        error |= is_wrong(out_fab(bit(), c_d1_gamma_UU12), local_vars.d1_gamma_UU[0][1][0], "c_d1_gammaUU12");
        error |= is_wrong(out_fab(bit(), c_d1_gamma_UU13), local_vars.d1_gamma_UU[0][2][0], "c_d1_gammaUU13");
        error |= is_wrong(out_fab(bit(), c_d1_gamma_UU22), local_vars.d1_gamma_UU[1][1][0], "c_d1_gammaUU22");
        error |= is_wrong(out_fab(bit(), c_d1_gamma_UU23), local_vars.d1_gamma_UU[1][2][0], "c_d1_gammaUU23");
        error |= is_wrong(out_fab(bit(), c_d1_gamma_UU33), local_vars.d1_gamma_UU[2][2][0], "c_d1_gammaUU33");

        error |= is_wrong(out_fab(bit(), c_d2_gamma11), local_vars.d2_gamma[0][0][0][1], "c_d2_gamma11");
        error |= is_wrong(out_fab(bit(), c_d2_gamma12), local_vars.d2_gamma[0][1][0][1], "c_d2_gamma11");
        error |= is_wrong(out_fab(bit(), c_d2_gamma13), local_vars.d2_gamma[0][2][0][1], "c_d2_gamma11");
        error |= is_wrong(out_fab(bit(), c_d2_gamma22), local_vars.d2_gamma[1][1][0][1], "c_d2_gamma11");
        error |= is_wrong(out_fab(bit(), c_d2_gamma23), local_vars.d2_gamma[1][2][0][1], "c_d2_gamma11");
        error |= is_wrong(out_fab(bit(), c_d2_gamma33), local_vars.d2_gamma[2][2][0][1], "c_d2_gamma11");

        error |= is_wrong(out_fab(bit(), c_d2_lapse), local_vars.d2_lapse[0][0], "c_d2_lapse");
        /*
        error |= is_wrong(out_fab(bit(), c_d1_shift1), local_vars.d1_shift[0][0], "c_d1_shift1");
        error |= is_wrong(out_fab(bit(), c_d1_shift2), local_vars.d1_shift[1][0], "c_d1_shift2");
        error |= is_wrong(out_fab(bit(), c_d1_shift3), local_vars.d1_shift[2][0], "c_d1_shift3");
        */

        error |= is_wrong(out_fab(bit(), c_d2_shift1), local_vars.d2_shift[0][0][0], "c_d2_shift1");
        error |= is_wrong(out_fab(bit(), c_d2_shift2), local_vars.d2_shift[1][0][0], "c_d2_shift2");
        error |= is_wrong(out_fab(bit(), c_d2_shift3), local_vars.d2_shift[2][0][0], "c_d2_shift3");

        error |= is_wrong(out_fab(bit(), c_d1_K_tensor11), local_vars.d1_K_tensor[0][0][2], "c_d1_K_tensor11");
        error |= is_wrong(out_fab(bit(), c_d1_K_tensor12), local_vars.d1_K_tensor[0][1][2], "c_d1_K_tensor12");
        error |= is_wrong(out_fab(bit(), c_d1_K_tensor13), local_vars.d1_K_tensor[0][2][2], "c_d1_K_tensor13");
        error |= is_wrong(out_fab(bit(), c_d1_K_tensor22), local_vars.d1_K_tensor[1][1][2], "c_d1_K_tensor22");
        error |= is_wrong(out_fab(bit(), c_d1_K_tensor23), local_vars.d1_K_tensor[1][2][2], "c_d1_K_tensor23");
        error |= is_wrong(out_fab(bit(), c_d1_K_tensor33), local_vars.d1_K_tensor[2][2][2], "c_d1_K_tensor33");


        //std::cout << out_fab(bit(), c_d1_gamma_UU111) << " gamma test" << endl;
        //std::cout << out_fab(bit(), c_d1_lapse1) << " numeric \n" << endl;
        //std::cout << local_vars.d1_lapse[0] << " analytic \n" << endl;
        if (error)
        {
            std::cout << "ADM unit tests NOT passed.\n";
            return error;
        }
    }

    // Sixth order derivatives
    BoxLoops::loop(ADMTestsCompute<SixthOrderDerivatives>(dx), in_fab,
                   out_fab);

    for (bit.begin(); bit.ok(); ++bit)
    {
        const double x = (0.5 + bit()[0]) * dx;
        const double z = (0.5 + bit()[2]) * dx;

        bool error = false;
        /*
        error |= is_wrong(out_fab(bit(), c_d1), 2 * x * (z - 0.5),
                          "diff1 (sixth order)");
        error |= is_wrong(out_fab(bit(), c_d2), 2 * x, "diff2 (sixth order)");
        error |= is_wrong(out_fab(bit(), c_d2_mixed), 2 * (z - 0.5),
                          "mixed diff2 (sixth order)");

        double correct_dissipation = (1. + z * (z - 1)) * pow(dx, 5) / 64;
        error |= is_wrong(out_fab(bit(), c_diss), correct_dissipation,
                          "dissipation (sixth order)");

        double correct_advec_down = -2 * z * (z - 1) - 3 * x * (2 * z - 1);
        error |= is_wrong(out_fab(bit(), c_advec_down), correct_advec_down,
                          "advection down (sixth order)");

        double correct_advec_up = 2 * z * (z - 1) + 3 * x * (2 * z - 1);
        error |= is_wrong(out_fab(bit(), c_advec_up), correct_advec_up,
                          "advection up (sixth order)");

        if (error)
        {
            std::cout << "ADM unit tests NOT passed.\n";
            return error;
        }
        */
    }

    std::cout << "ADM unit tests passed.\n";
    return 0;
}
