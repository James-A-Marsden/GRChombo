/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "BoostedHeadOnLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "CCZ4RHS.hpp"

#include "ChiExtractionTaggingCriterion.hpp"
#include "ChiPunctureExtractionTaggingCriterion.hpp"
//Refinement threshold for the boosted BH case
//#include "BoostedPunctureTrackerTaggingCriterion.hpp"

#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "NewConstraints.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "PunctureTracker.hpp"
#include "SetValue.hpp"
#include "SixthOrderDerivatives.hpp"
#include "SmallDataIO.hpp"
#include "TraceARemoval.hpp"
#include "TwoPuncturesInitialData.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"

//Initial conditions
#include "IsotropicBoostedBH_bk.hpp"

// Things to do during the advance step after RK4 steps
void BoostedHeadOnLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck("NaNCheck in specific Advance: "), m_state_new,
                       m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

// This initial data uses an approximation for the metric which
// is valid for small boosts
void BoostedHeadOnLevel::initialData()
{

    CH_TIME("HeadOnLevel::initialData");
    if (m_verbosity)
        pout() << "HeadOnLevel::initialData " << m_level << endl;

        IsotropicBoostedBH isotropicBoostedBH(m_p.bh1_params, m_p.bh2_params,
                                            m_dx);
    BoxLoops::loop(make_compute_pack(SetValue(0.0), isotropicBoostedBH),
                    m_state_new, m_state_new, INCLUDE_GHOST_CELLS,disable_simd());

}

// Calculate RHS during RK4 substeps
void BoostedHeadOnLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                    const double a_time)
{
    // Enforce positive chi and alpha and trace free A
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate CCZ4 right hand side
    if (m_p.max_spatial_derivative_order == 4)
    {
        BoxLoops::loop(CCZ4RHS<MovingPunctureGauge, FourthOrderDerivatives>(
                           m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation),
                       a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        BoxLoops::loop(CCZ4RHS<MovingPunctureGauge, SixthOrderDerivatives>(
                           m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation),
                       a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// enforce trace removal during RK4 substeps
void BoostedHeadOnLevel::specificUpdateODE(GRLevelData &a_soln,
                                      const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void BoostedHeadOnLevel::preTagCells()
{
    // We only use chi in the tagging criterion so only fill the ghosts for chi
    fillAllGhosts(VariableType::evolution, Interval(c_chi, c_chi));
}

// specify the cells to tag
void BoostedHeadOnLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)
{
    
    if (m_p.track_punctures)
    {
        std::vector<double> puncture_masses;
        #ifdef USE_TWOPUNCTURES
        // use calculated bare masses from TwoPunctures
        puncture_masses = {m_tp_amr.m_two_punctures.mm,
                           m_tp_amr.m_two_punctures.mp};
        //puncture_masses = {m_p.bh1_params.mass, m_p.bh2_params.mass};
        //puncture_momentum = {m_p.bh1_params.momentum[0],
                            // m_p.bh2_params.momentum[0]};

        //punctures = m_bh_amr.m_puncture_tracker.get_puncture_coords();                           
        #else
        puncture_masses = {m_p.bh1_params.mass, m_p.bh2_params.mass};
        #endif /* USE_TWOPUNCTURES */

        auto puncture_coords =
            m_bh_amr.m_puncture_tracker.get_puncture_coords();
            
        BoxLoops::loop(ChiPunctureExtractionTaggingCriterion(
                           m_dx, m_level, m_p.max_level, m_p.extraction_params,
                           puncture_coords, m_p.activate_extraction,
                           m_p.track_punctures, puncture_masses),
                       current_state, tagging_criterion);
                       
        /*
        BoxLoops::loop(BoostedPunctureTrackerTaggingCriterion<FourthOrderDerivatives>(
                m_dx, m_level, m_p.puncture_tag_max_levels, m_p.extraction_params,
                punctures, m_p.activate_extraction, m_p.track_punctures,
                puncture_masses, puncture_momentum, m_p.bh_tagging_buffer,
                m_p.puncture_tag_min_separation),
                current_state, tagging_criterion);
        */
    }
    else
    {
        BoxLoops::loop(ChiExtractionTaggingCriterion(m_dx, m_level,
                                                     m_p.extraction_params,
                                                     m_p.activate_extraction),
                       current_state, tagging_criterion);
    }
}

void BoostedHeadOnLevel::specificPostTimeStep()
{
    CH_TIME("BoostedHeadOnLevel::specificPostTimeStep");

    bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main
    // bool first_step = (m_time == m_dt); // if not called in Main

    if (m_p.activate_extraction == 1)
    {
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);
        if (calculate_weyl)
        {
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();
            BoxLoops::loop(
                Weyl4(m_p.extraction_params.center, m_dx, m_p.formulation),
                m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

            // Do the extraction on the min extraction level
            if (m_level == min_level)
            {
                CH_TIME("WeylExtraction");
                // Now refresh the interpolator and do the interpolation
                // fill ghosts manually to minimise communication
                bool fill_ghosts = false;
                m_gr_amr.m_interpolator->refresh(fill_ghosts);
                m_gr_amr.fill_multilevel_ghosts(
                    VariableType::diagnostic, Interval(c_Weyl4_Re, c_Weyl4_Im),
                    min_level);
                WeylExtraction my_extraction(m_p.extraction_params, m_dt,
                                             m_time, first_step,
                                             m_restart_time);
                my_extraction.execute_query(m_gr_amr.m_interpolator);
            }
        }
    }

    if (m_p.calculate_constraint_norms)
    {
        fillAllGhosts();
        BoxLoops::loop(Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3)),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        if (m_level == 0)
        {
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
            double L2_Ham = amr_reductions.norm(c_Ham);
            double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3));
            SmallDataIO constraints_file(m_p.data_path + "constraint_norms",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            constraints_file.remove_duplicate_time_data();
            if (first_step)
            {
                constraints_file.write_header_line({"L^2_Ham", "L^2_Mom"});
            }
            constraints_file.write_time_data_line({L2_Ham, L2_Mom});
        }
    }

    // do puncture tracking on requested level
    if (m_p.track_punctures && m_level == m_p.puncture_tracking_level)
    {
        CH_TIME("PunctureTracking");
        // only do the write out for every coarsest level timestep
        int coarsest_level = 0;
        bool write_punctures = at_level_timestep_multiple(coarsest_level);
        m_bh_amr.m_puncture_tracker.execute_tracking(m_time, m_restart_time,
                                                     m_dt, write_punctures);
    }
}

#ifdef CH_USE_HDF5
// Things to do before a plot level - need to calculate the Weyl scalars
void BoostedHeadOnLevel::prePlotLevel()
{
    fillAllGhosts();
    if (m_p.activate_extraction == 1)
    {
        BoxLoops::loop(
            make_compute_pack(
                Weyl4(m_p.extraction_params.center, m_dx, m_p.formulation),
                Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3))),
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    }
}
#endif /* CH_USE_HDF5 */