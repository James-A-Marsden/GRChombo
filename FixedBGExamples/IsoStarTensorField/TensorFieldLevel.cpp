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
#include "TensorPotential.hpp"
#include "ExcisionDiagnostics.hpp"
#include "ExcisionEvolution.hpp"
#include "FixedBGTensorField.hpp"
#include "FixedBGDiagnostics.hpp"
#include "FixedBGFluxes.hpp"
#include "FixedBGEvolution.hpp"
#include "FluxExtraction.hpp"
#include "InitialConditions.hpp"
#include "Sethbar.hpp"
#include "SetRest.hpp"
#include "IsoStarFixedBG.hpp"
#include "TraceFieldRemoval.hpp"

// Things to do at each advance step, after the RK4 is calculated
void TensorFieldLevel::specificAdvance()
{

    //Enforce the spatial rank 2 field component to be traceless
    IsoStarFixedBG isostar_bg(m_p.bg_params, m_dx); 
    BoxLoops::loop(TraceFieldRemoval<TensorFieldWithPotential, IsoStarFixedBG>(
        isostar_bg),m_state_new, m_state_new, SKIP_GHOST_CELLS, disable_simd());
                   
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
    IsoStarFixedBG isostar_bg(m_p.bg_params, m_dx); // just calculates chi
    //m_p.field_amplitude_re, m_p.field_amplitude_im,

    InitialConditions set_vars(m_p.potential_params.tensor_mass, m_p.center,
                              m_p.bg_params, m_dx, m_p.initial_constant);

    Sethbar set_hbar(m_p.potential_params.tensor_mass, m_p.center,
                            m_p.bg_params, m_dx, m_p.initial_constant);

    SetRest set_rest(m_p.potential_params.tensor_mass, m_p.center,
    m_p.bg_params, m_dx, m_p.initial_constant);

    auto compute_pack = make_compute_pack(set_zero, isostar_bg);

    BoxLoops::loop(compute_pack, m_state_diagnostics, m_state_diagnostics,
                   SKIP_GHOST_CELLS);
   
    //NO EXCISION 
    
    BoxLoops::loop(set_vars, m_state_new, m_state_new, INCLUDE_GHOST_CELLS, disable_simd());
    //fillAllGhosts();
    //BoxLoops::loop(set_hbar, m_state_new, m_state_new, INCLUDE_GHOST_CELLS, disable_simd());
    //fillAllGhosts();
    //BoxLoops::loop(set_rest, m_state_new, m_state_new, INCLUDE_GHOST_CELLS, disable_simd());

    // excise within horizon, no simd
    BoxLoops::loop(
        ExcisionEvolution<TensorFieldWithPotential, IsoStarFixedBG>(
            m_dx, m_p.center, isostar_bg),
        m_state_new, m_state_new, SKIP_GHOST_CELLS, disable_simd());


    ///DIAGNOSTICS 
}

// Things to do before outputting a plot file
void TensorFieldLevel::prePlotLevel() {

    fillAllGhosts();
    TensorPotential potential(m_p.potential_params);
    TensorFieldWithPotential tensor_field(potential);
    IsoStarFixedBG isostar_bg(m_p.bg_params, m_dx);

    FixedBGDiagnostics<TensorFieldWithPotential, IsoStarFixedBG>
        diagnostics(tensor_field, isostar_bg, m_dx, m_p.center);
    FixedBGFluxes<TensorFieldWithPotential,
                                    IsoStarFixedBG>
        energy_fluxes(tensor_field, isostar_bg, m_dx, m_p.center);
    BoxLoops::loop(make_compute_pack(diagnostics, energy_fluxes),
                    m_state_new, m_state_diagnostics, SKIP_GHOST_CELLS, disable_simd());
    // excise within horizon
    BoxLoops::loop(
        ExcisionDiagnostics<TensorFieldWithPotential, IsoStarFixedBG>(
            m_dx, m_p.center, isostar_bg, m_p.inner_r, m_p.outer_r),
        m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
        disable_simd());     
}

// Things to do in RHS update, at each RK4 step
void TensorFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{

    IsoStarFixedBG isostar_bg(m_p.bg_params, m_dx);
    BoxLoops::loop(TraceFieldRemoval<TensorFieldWithPotential, IsoStarFixedBG>(
        isostar_bg),a_soln, a_soln, SKIP_GHOST_CELLS, disable_simd());

    // Calculate MatterCCZ4 right hand side with matter_t = TensorField
    // We don't want undefined values floating around in the constraints so
    // zero these

    // enforce continuous prescription for sigma as per arXiv:2104.06978
    //const int ratio = pow(2, 5 * (m_level - m_p.max_level));
    //const double sigma = m_p.sigma * double(ratio);

    TensorPotential potential(m_p.potential_params);
    TensorFieldWithPotential tensor_field(potential);
    //TensorFieldWithPotential;
    //IsoStarFixedBG isostar_bg(m_p.bg_params, m_dx);
    FixedBGEvolution<TensorFieldWithPotential, IsoStarFixedBG> my_matter(
        tensor_field, isostar_bg, m_p.sigma, m_dx, m_p.center);
        //tensor_field, isostar_bg, sigma, m_dx, m_p.center);

    //TraceFieldRemoval<TensorFieldWithPotential, IsoStarFixedBG> make_traceless(isostar_bg); 

    BoxLoops::loop(my_matter, a_soln, a_rhs, SKIP_GHOST_CELLS);

    //Enforce traceless fspatial
    //BoxLoops::loop(make_compute_pack(make_traceless), a_soln, a_soln, INCLUDE_GHOST_CELLS);    


    // excise within horizon, no simd
    BoxLoops::loop(
        ExcisionEvolution<TensorFieldWithPotential, IsoStarFixedBG>(
            m_dx, m_p.center, isostar_bg),
        a_soln, a_rhs, SKIP_GHOST_CELLS, disable_simd());

    
}

void TensorFieldLevel::specificPostTimeStep()
{
    // At any level, but after the coarsest timestep
    double coarsest_dt = m_p.coarsest_dx * m_p.dt_multiplier;
    const double remainder = fmod(m_time, coarsest_dt);
    if (min(abs(remainder), abs(remainder - coarsest_dt)) < 1.0e-8)
    {
        // calculate the density of the PF, but excise the BH region completely
        fillAllGhosts();
        TensorPotential potential(m_p.potential_params);
        TensorFieldWithPotential tensor_field(potential);
        IsoStarFixedBG isostar_bg(m_p.bg_params, m_dx);

        FixedBGDiagnostics<TensorFieldWithPotential, IsoStarFixedBG>
            diagnostics(tensor_field, isostar_bg, m_dx, m_p.center);

        FixedBGFluxes<TensorFieldWithPotential,
                                       IsoStarFixedBG>
            energy_fluxes(tensor_field, isostar_bg, m_dx, m_p.center);

        BoxLoops::loop(make_compute_pack(diagnostics, energy_fluxes),
                       m_state_new, m_state_diagnostics, SKIP_GHOST_CELLS, disable_simd());

        // excise within horizon
        BoxLoops::loop(
            ExcisionDiagnostics<TensorFieldWithPotential, IsoStarFixedBG>(
                m_dx, m_p.center, isostar_bg, m_p.inner_r, m_p.outer_r),
            m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
            disable_simd());
    }

    // write out the integral after each coarse timestep
    /*
    if (m_level == 0)
    {
        bool first_step = (m_time == m_dt);

        // integrate the densities and write to a file
        AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
        double rho_sum = amr_reductions.sum(c_rho);
        double rhoJ_sum = amr_reductions.sum(c_rhoJ);

        SmallDataIO integral_file(m_p.integral_filename, m_dt, m_time,
                                  m_restart_time, SmallDataIO::APPEND,
                                  first_step);
        // remove any duplicate data if this is post restart
        integral_file.remove_duplicate_time_data();
        std::vector<double> data_for_writing = {rho_sum, rhoJ_sum};
        // write data
        if (first_step)
        {
            integral_file.write_header_line({"rho", "rhoJ"});
        }
        integral_file.write_time_data_line(data_for_writing);

        // Now refresh the interpolator and do the interpolation
        m_gr_amr.m_interpolator->refresh();
        FluxExtraction my_extraction(m_p.extraction_params, m_dt, m_time,
                                      m_restart_time);
        my_extraction.execute_query(m_gr_amr.m_interpolator);
    }
    */
}

// enforce trace removal during RK4 substeps
void TensorFieldLevel::specificUpdateODE(GRLevelData &a_soln,
                                      const GRLevelData &a_rhs, Real a_dt)
{
    IsoStarFixedBG isostar_bg(m_p.bg_params, m_dx);
    BoxLoops::loop(TraceFieldRemoval<TensorFieldWithPotential, IsoStarFixedBG>(
        isostar_bg),a_soln, a_soln, SKIP_GHOST_CELLS, disable_simd());

}



void TensorFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level, m_p.regrid_length,
                                              m_p.center),
                   current_state, tagging_criterion, disable_simd());
}
