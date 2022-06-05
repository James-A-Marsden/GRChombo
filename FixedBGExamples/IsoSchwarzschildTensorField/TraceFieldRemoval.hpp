/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This class enforces A to be trace-free
#ifndef TRACEFIELDREMOVAL_HPP_
#define TRACEFIELDREMOVAL_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"
#include "FixedBGTensorField.hpp"
#include "ADMFixedBGVars.hpp"

template <class matter_t, class background_t> class TraceFieldRemoval
{

    template <class data_t> using Vars = typename matter_t::template Vars<data_t>;
    // Use the variable definition in the matter class
    //template <class data_t>
    //using MatterVars = typename matter_t::template Vars<data_t>;
    //template <class data_t> using Vars = FixedBGTensorField<TensorPotential>::template Vars<data_t>;
    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  protected:
    const background_t m_background; //the background
    //const vars_t vars; // The field variables


  public:

    TraceFieldRemoval(background_t a_background)
        : m_background(a_background)
    {
    }


    //void compute(const Cell<double> current_cell) const
    template <class data_t> void compute(const Cell<data_t> current_cell) const
    {
       
     // Struct for the non grid ADM vars

        //template <class data_t> using MetricVars = typename ADMFixedBGVars::template Vars<data_t>; 

        Vars<data_t> vars;
        //Tensor<2,data_t> fspatial; //Spatial component of the tensor field
        //Tensor<2,data_t> v; //Spatial rank 2 v field
        //const auto local_vars = current_cell.template load_vars<Vars>();
        //Load metric vars

        MetricVars<double> metric_vars;

        m_background.compute_metric_background(metric_vars, current_cell);
        const auto gamma_UU = TensorAlgebra::compute_inverse_sym(metric_vars.gamma);


        //Make fspatial trace free
        auto local_vars = current_cell.template load_vars<Vars>();
        TensorAlgebra::make_trace_free(local_vars.fspatial, metric_vars.gamma, gamma_UU);
        TensorAlgebra::make_trace_free(local_vars.v, metric_vars.gamma, gamma_UU);
        
        auto ftrace = TensorAlgebra::compute_trace(local_vars.fspatial, gamma_UU);
        auto vtrace = TensorAlgebra::compute_trace(local_vars.v, gamma_UU);
        //double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
 
        //local_vars.fspatial[0][0] -= (1.0 / 2.0) * metric_vars.gamma[0][0] * ftrace;
        //local_vars.fspatial[1][1] -= (1.0 / 2.0) * metric_vars.gamma[1][1] * ftrace;

        //local_vars.v[0][0] -= (1.0 / 2.0) * metric_vars.gamma[0][0] * vtrace;
        //local_vars.v[1][1] -= (1.0 / 2.0) * metric_vars.gamma[1][1] * vtrace;

        vars.fhat = local_vars.fhat;
        vars.w = local_vars.w; 
        FOR1(i)
        {
          vars.fbar[i] = local_vars.fbar[i];
          vars.q[i] = local_vars.q[i];

          FOR1(j)
          {
            vars.fspatial[i][j] = local_vars.fspatial[i][j];
            vars.v[i][j] = local_vars.v[i][j];
          }
        }

        current_cell.store_vars(vars);
    }
};

/*
template <class data_t>
template <typename mapping_function_t>
void TraceFieldRemoval::Vars<data_t>::enum_mapping(
    mapping_function_t mapping_function)
{
    VarsTools::define_symmetric_enum_mapping(mapping_function, GRInterval<c_fspatial11, c_fspatial33>(), h);
    VarsTools::define_symmetric_enum_mapping(mapping_function, c_fhat , fhat);
}
*/
#endif /* TRACEFIELDREMOVAL_HPP_ */