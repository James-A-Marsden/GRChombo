/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TENSORFIELDLEVEL_HPP_
#define TENSORFIELDLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// Problem specific includes
#include "FixedBGTensorField.hpp"
#include "TensorPotential.hpp"

//!  A class for the evolution of a tensor field, minimally coupled to gravity
/*!
     The class takes some initial data for a tensor field (variables phi and Pi)
     and evolves it using the CCZ4 equations. It is possible to specify an
   initial period of relaxation for the conformal factor chi, for non analytic
   initial conditions (for example, a general field configuration at a moment of
   time symmetry assuming conformal flatness). \sa MatterCCZ4(),
   ConstraintsMatter(), TensorField(), RelaxationChi()
*/
class TensorFieldLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<TensorFieldLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // Typedef for tensor field
    typedef FixedBGTensorField<TensorPotential> TensorFieldWithPotential;
    // typedef FixedBGTensorField
    //     TensorFieldWithPotential;
    //! Things to do at the end of the advance step, after RK4 calculation
    virtual void specificAdvance();

    //! Initialize data for the field and metric variables
    virtual void initialData();

    //! routines to do before outputing plot file
    virtual void prePlotLevel();

    //! RHS routines used at each RK4 step
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time);

    /// Things to do after dt*rhs has been added to the solution
    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs,
                                   Real a_dt) override;

    //! To do after each timestep
    virtual void specificPostTimeStep();

    //! Tell Chombo how to tag cells for regridding
    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state);
};

#endif /* TENSORFIELDLEVEL_HPP_ */
