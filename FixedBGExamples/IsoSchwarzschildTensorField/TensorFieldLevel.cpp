/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "TensorFieldLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"

// For tag cells
#include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "CustomExtraction.hpp"
#include "ExcisionDiagnostics.hpp"
#include "ExcisionEvolution.hpp"
#include "FixedBGEvolution.hpp"
#include "FixedBGTensorField.hpp"
#include "InitialConditions.hpp"
#include "IsoSchwarzschildFixedBG.hpp"
#include "TensorDiagnostics.hpp"

// Initial data for field and metric variables
void TensorFieldLevel::initialData()
{
    CH_TIME("TensorFieldLevel::initialData");
    if (m_verbosity)
        pout() << "TensorFieldLevel::initialData " << m_level << endl;

    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then initial conditions for fields
    SetValue set_zero(0.0);
    IsoSchwarzschildFixedBG isoschwarzschild_bg(m_p.bg_params,
                                                m_dx); // just calculates chi

    auto compute_pack = make_compute_pack(set_zero, isoschwarzschild_bg);

    BoxLoops::loop(compute_pack, m_state_diagnostics, m_state_diagnostics,
                   INCLUDE_GHOST_CELLS);

    InitialConditions set_vars(m_p.center, m_p.bg_params, m_dx);

    BoxLoops::loop(set_zero, m_state_new, m_state_new, INCLUDE_GHOST_CELLS);
    BoxLoops::loop(set_vars, m_state_new, m_state_new, INCLUDE_GHOST_CELLS,
                   disable_simd());

    // excise within horizon, no simd
    BoxLoops::loop(
        ExcisionEvolution<FixedBGTensorField, IsoSchwarzschildFixedBG>(
            m_dx, m_p.center, isoschwarzschild_bg),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS, disable_simd());
}

// Things to do before outputting a plot file
void TensorFieldLevel::prePlotLevel()
{

    CH_TIME("TensorFieldLevel::prePlotLevel");
    // Now done in PostTimestep below
}

// Things to do in RHS update, at each RK4 step
void TensorFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    CH_TIME("TensorFieldLevel::specificEvalRHS");
    IsoSchwarzschildFixedBG isoschwarzschild_bg(m_p.bg_params, m_dx);
    FixedBGTensorField tensor_field(
        m_p.tensor_field_mass, m_p.damping_kappa, m_p.damping_is_active,
        m_p.dRGT_ij_is_active, m_p.dRGT_mass_is_active);

    // enforce continuous prescription for sigma as per arXiv:2104.06978
    // const int ratio = pow(2, 5 * (m_level - m_p.max_level));
    // const double sigma = m_p.sigma * double(ratio);

    FixedBGEvolution<FixedBGTensorField, IsoSchwarzschildFixedBG> my_matter(
        tensor_field, isoschwarzschild_bg, m_p.sigma, m_dx, m_p.center);

    // Calculate!
    BoxLoops::loop(my_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);

    // excise within horizon, no simd
    BoxLoops::loop(
        ExcisionEvolution<FixedBGTensorField, IsoSchwarzschildFixedBG>(
            m_dx, m_p.center, isoschwarzschild_bg),
        a_soln, a_rhs, INCLUDE_GHOST_CELLS, disable_simd());
}

void TensorFieldLevel::specificPostTimeStep()
{
    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new, SKIP_GHOST_CELLS,
                       disable_simd());

    int min_level = 0;
    bool calculate_diagnostics = at_level_timestep_multiple(min_level);
    if (calculate_diagnostics)
    {
        // Calculate the constraints
        fillAllEvolutionGhosts();

        FixedBGTensorField tensor_field(
            m_p.tensor_field_mass, m_p.damping_kappa, m_p.damping_is_active,
            m_p.dRGT_ij_is_active, m_p.dRGT_mass_is_active);
        IsoSchwarzschildFixedBG isoschwarzschild_bg(m_p.bg_params, m_dx);

        TensorDiagnostics<FixedBGTensorField, IsoSchwarzschildFixedBG>
            diagnostics(tensor_field, isoschwarzschild_bg, m_dx, m_p.center,
                        m_p.tensor_field_mass);

        // Calculate diagnostics
        BoxLoops::loop(diagnostics, m_state_new, m_state_diagnostics,
                       SKIP_GHOST_CELLS);

        // excise within horizon
        BoxLoops::loop(
            ExcisionDiagnostics<FixedBGTensorField, IsoSchwarzschildFixedBG>(
                m_dx, m_p.center, isoschwarzschild_bg, m_p.inner_r,
                m_p.outer_r),
            m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
            disable_simd());

        if (m_level == 0)
        {
            bool first_step = (m_time == 0.0);

            // Output the constraints integral
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
            double L2_Ham = amr_reductions.norm(c_primaryConstraintScalar);
            double L2_Mom = amr_reductions.norm(Interval(
                c_primaryConstraintVector1, c_primaryConstraintVector3));
            SmallDataIO constraints_file(m_p.output_path + "constraint_norms",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            constraints_file.remove_duplicate_time_data();
            if (first_step)
            {
                constraints_file.write_header_line({"L^2_Ham", "L^2_Mom"});
            }
            constraints_file.write_time_data_line({L2_Ham, L2_Mom});

            // Extract data along a line
            AMRInterpolator<Lagrange<4>> interpolator(
                m_gr_amr, m_p.origin, m_p.dx, m_p.boundary_params,
                m_p.verbosity);
            bool fill_ghosts = false;
            interpolator.refresh(fill_ghosts);
            m_gr_amr.fill_multilevel_ghosts(
                VariableType::evolution, Interval(c_fspatial11, c_fspatial11),
                min_level);
            int num_points = 100;
            CustomExtraction extraction(c_fspatial11, num_points, m_p.L,
                                        m_p.center, m_dt, m_time);
            extraction.execute_query(m_gr_amr.m_interpolator, "outputs");
        }
    }
}

void TensorFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level, m_p.regrid_length,
                                              m_p.center),
                   current_state, tagging_criterion);
}
