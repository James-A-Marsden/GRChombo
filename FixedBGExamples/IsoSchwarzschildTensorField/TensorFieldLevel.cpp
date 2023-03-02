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
#include "TensorFieldTaggingCriterion.hpp"

// Problem specific includes
#include "ExcisionDiagnostics.hpp"
#include "ExcisionEvolution.hpp"
#include "FixedBGDiagnostics.hpp"
#include "FixedBGEvolution.hpp"
#include "FixedBGFluxes.hpp"
#include "FixedBGTensorField.hpp"
#include "FluxExtraction.hpp"
#include "InitialConditions.hpp"
#include "IsoSchwarzschildFixedBG.hpp"
#include "SetRest.hpp"
#include "TensorPotential.hpp"
#include "TraceFieldRemoval.hpp"

// Things to do at each advance step, after the RK4 is calculated
void TensorFieldLevel::specificAdvance()
{

//    IsoSchwarzschildFixedBG isoschwarzschild_bg(m_p.bg_params, m_dx);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new, SKIP_GHOST_CELLS,
                       disable_simd());
}

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

    InitialConditions set_vars(m_p.potential_params.tensor_mass, m_p.center,
                               m_p.bg_params, m_dx, m_p.initial_constant);

    BoxLoops::loop(set_zero, m_state_new, m_state_new, INCLUDE_GHOST_CELLS);
    BoxLoops::loop(set_vars, m_state_new, m_state_new, INCLUDE_GHOST_CELLS,
                   disable_simd());

    // excise within horizon, no simd
    BoxLoops::loop(
        ExcisionEvolution<TensorFieldWithPotential, IsoSchwarzschildFixedBG>(
            m_dx, m_p.center, isoschwarzschild_bg),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS, disable_simd());
}

// Things to do before outputting a plot file
void TensorFieldLevel::prePlotLevel()
{

    fillAllGhosts();
    TensorPotential potential(m_p.potential_params);
    TensorFieldWithPotential tensor_field(potential);
    IsoSchwarzschildFixedBG isoschwarzschild_bg(m_p.bg_params, m_dx);

    FixedBGDiagnostics<TensorFieldWithPotential, IsoSchwarzschildFixedBG>
        diagnostics(tensor_field, isoschwarzschild_bg, m_dx, m_p.center);
    FixedBGFluxes<TensorFieldWithPotential, IsoSchwarzschildFixedBG>
        energy_fluxes(tensor_field, isoschwarzschild_bg, m_dx, m_p.center);

    // Enforce traceless Hmunu
    BoxLoops::loop(
        TraceFieldRemoval<TensorFieldWithPotential, IsoSchwarzschildFixedBG>(
            isoschwarzschild_bg, m_dx),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS, disable_simd());

    // Calculate diagnostics
    BoxLoops::loop(make_compute_pack(diagnostics, energy_fluxes), m_state_new,
                   m_state_diagnostics, SKIP_GHOST_CELLS, disable_simd());

    // excise within horizon
    BoxLoops::loop(
        ExcisionDiagnostics<TensorFieldWithPotential, IsoSchwarzschildFixedBG>(
            m_dx, m_p.center, isoschwarzschild_bg, m_p.inner_r, m_p.outer_r),
        m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
        disable_simd());
    fillAllGhosts();

}

// Things to do in RHS update, at each RK4 step
void TensorFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    SetValue set_zero(0.0);
    IsoSchwarzschildFixedBG isoschwarzschild_bg(m_p.bg_params, m_dx);
    TensorPotential potential(m_p.potential_params);
    TensorFieldWithPotential tensor_field(potential);
    FixedBGEvolution<TensorFieldWithPotential, IsoSchwarzschildFixedBG>
        my_matter(tensor_field, isoschwarzschild_bg, m_p.sigma, m_dx,
                  m_p.center);
    TraceFieldRemoval<TensorFieldWithPotential, IsoSchwarzschildFixedBG>
        make_traceless(isoschwarzschild_bg, m_dx);

    BoxLoops::loop(set_zero, a_soln, a_rhs, INCLUDE_GHOST_CELLS);
    BoxLoops::loop(my_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);

    // excise within horizon, no simd
    BoxLoops::loop(
        ExcisionEvolution<TensorFieldWithPotential, IsoSchwarzschildFixedBG>(
            m_dx, m_p.center, isoschwarzschild_bg),
        a_soln, a_rhs, INCLUDE_GHOST_CELLS, disable_simd());

    // Enforce traceless Hmunu
    // BoxLoops::loop(make_compute_pack(make_traceless), a_soln, a_soln,
    //               INCLUDE_GHOST_CELLS);
}

void TensorFieldLevel::specificPostTimeStep() {}

// enforce trace removal during RK4 substeps
void TensorFieldLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
/*
    IsoSchwarzschildFixedBG isoschwarzschild_bg(m_p.bg_params, m_dx);
    BoxLoops::loop(
        TraceFieldRemoval<TensorFieldWithPotential, IsoSchwarzschildFixedBG>(
            isoschwarzschild_bg, m_dx),
        a_soln, a_soln, INCLUDE_GHOST_CELLS, disable_simd());
*/
}

void TensorFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(TensorFieldTaggingCriterion(m_dx, m_level, m_p.regrid_length,
                                               m_p.center),
                   current_state, tagging_criterion, disable_simd());
}
