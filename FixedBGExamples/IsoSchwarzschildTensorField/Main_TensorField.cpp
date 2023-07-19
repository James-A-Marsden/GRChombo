/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "parstream.H" //Gives us pout()
#include <iostream>

#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "MultiLevelTask.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "TensorFieldLevel.hpp"

int runGRChombo(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    char *in_file = argv[1];
    GRParmParse pp(argc - 2, argv + 2, NULL, in_file);
    SimulationParameters sim_params(pp);

    // The line below selects the problem that is simulated
    // (To simulate a different problem, define a new child of AMRLevel
    // and an associated LevelFactory)
    GRAMR gr_amr;
    DefaultLevelFactory<TensorFieldLevel> tensor_field_level_fact(gr_amr,
                                                                  sim_params);
    setupAMRObject(gr_amr, tensor_field_level_fact);

    // call this after amr object setup so grids known
    // and need it to stay in scope throughout run
    AMRInterpolator<Lagrange<4>> interpolator(
        gr_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);
    gr_amr.set_interpolator(&interpolator);

    // Run specificPostTimeStep at t=0
    auto task = [](GRAMRLevel *level)
    {
        if (level->time() == 0.)
            level->specificPostTimeStep();
    };

    // call 'now' if not restarting
    if (!pp.contains("restart_file"))
    {
        MultiLevelTaskPtr<> call_task(task);
        call_task.execute(gr_amr);
    }

    // Engage! Run the evolution
    gr_amr.run(sim_params.stop_time, sim_params.max_steps);
    gr_amr.conclude();

    CH_TIMER_REPORT();

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runGRChombo(argc, argv);

    if (status == 0)
        pout() << "GRChombo finished." << std::endl;
    else
        pout() << "GRChombo failed with return code " << status << std::endl;

    mainFinalize();
    return status;
}
