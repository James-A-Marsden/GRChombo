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
    fillAllGhosts();

    FixedBGTensorField tensor_field(
        m_p.tensor_field_mass, m_p.damping_kappa, m_p.damping_is_active,
        m_p.dRGT_ij_is_active, m_p.dRGT_mass_is_active);
    IsoSchwarzschildFixedBG isoschwarzschild_bg(m_p.bg_params, m_dx);

    TensorDiagnostics<FixedBGTensorField, IsoSchwarzschildFixedBG> diagnostics(
        tensor_field, isoschwarzschild_bg, m_dx, m_p.center,
        m_p.tensor_field_mass);

    // Calculate diagnostics
    BoxLoops::loop(diagnostics, m_state_new, m_state_diagnostics,
                   SKIP_GHOST_CELLS, disable_simd());

    // excise within horizon
    BoxLoops::loop(
        ExcisionDiagnostics<FixedBGTensorField, IsoSchwarzschildFixedBG>(
            m_dx, m_p.center, isoschwarzschild_bg, m_p.inner_r, m_p.outer_r),
        m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
        disable_simd());
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
}

void TensorFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level, m_p.regrid_length,
                                              m_p.center),
                   current_state, tagging_criterion, disable_simd());
}
