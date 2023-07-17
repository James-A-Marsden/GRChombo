/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "FixedBGSimulationParametersBase.hpp"
#include "GRParmParse.hpp"
// Problem specific includes:
#include "KerrSchildFixedBG.hpp"
#include "TensorPotential.hpp"

class SimulationParameters : public FixedBGSimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : FixedBGSimulationParametersBase(pp)
    {
        // read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // for regridding
        pp.load("regrid_length", regrid_length, L);
        pp.load("integral_filename", integral_filename);

        // Initial and Kerr data
        pp.load("bh_mass", bg_params.mass, 1.0);
        //        pp.load("bh_velocity", bg_params.velocity, 0.0);
        pp.load("bh_spin", bg_params.spin, 0.0);
        pp.load("bh_center", bg_params.center, center);
        // pp.load("field_amplitude_re", field_amplitude_re);
        // pp.load("field_amplitude_im", field_amplitude_im);
        pp.load("tensor_mass", potential_params.tensor_mass);
        pp.load("inner_r", inner_r, 1.0);
        pp.load("outer_r", outer_r, 0.75 * L);
        pp.load("initial_constant", initial_constant, 1.0);
        // pp.load("tensor_field_mass", tensor_mass, 1);
    }

    // Problem specific parameters
    // double field_amplitude_re, field_amplitude_im, regrid_length;
    double regrid_length;
    double inner_r, outer_r;
    double initial_constant;
    // double tensor_mass;
    std::string integral_filename;
    // Collection of parameters necessary for the sims
    KerrSchildFixedBG::params_t bg_params;
    TensorPotential::params_t potential_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
