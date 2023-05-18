/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FIXEDBGTENSORFIELD_HPP_
#define FIXEDBGTENSORFIELD_HPP_

#include "ADMFixedBGVars.hpp"
#include "CCZ4Geometry.hpp"
#include "DefaultPotential.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "TensorPotential.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS, total num of components
#include "VarsTools.hpp"

template <class potential_t = TensorPotential> class FixedBGTensorField
{
  protected:
    //! The local copy of the potential
    potential_t my_potential;

    double m_tensor_field_mass;
    double m_damping_kappa;
    bool m_damping_is_active;
    bool m_dRGT_ij_is_active;
    bool m_dRGT_mass_is_active;

  public:
    //!  Constructor of class FixedBGTensorField, inputs are the matter
    // //!  parameters.
    FixedBGTensorField(const potential_t a_potential,
                       double a_tensor_field_mass, double a_damping_kappa,
                       bool a_damping_is_active, bool a_dRGT_ij_is_active,
                       bool a_dRGT_mass_is_active)
        : my_potential(a_potential), m_tensor_field_mass(a_tensor_field_mass),
          m_damping_kappa(a_damping_kappa),
          m_damping_is_active(a_damping_is_active),
          m_dRGT_ij_is_active(a_dRGT_ij_is_active),
          m_dRGT_mass_is_active(a_dRGT_mass_is_active)
    {
    }
    //! Structure containing the rhs variables for the matter fields
    template <class data_t> struct Vars
    {

        // Tensor field components
        Tensor<2, data_t> fspatial; // Spatial component of the tensor field
        Tensor<1, data_t> fbar;
        data_t fhat;

        // Conjugate components
        Tensor<2, data_t> v; // Spatial rank 2 v field

        // Constraint damping variables
        Tensor<1, data_t> Xspatial; // Rank 1
        data_t Xhat;                // Scalar

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            // Map the components defined in user parameters to their full
            // objects. Field variables
            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_fspatial11, c_fspatial33>(),
                fspatial);

            VarsTools::define_enum_mapping(
                mapping_function, GRInterval<c_fbar1, c_fbar3>(), fbar);

            VarsTools::define_enum_mapping(mapping_function, c_fhat, fhat);
            // conjugate variables

            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_v11, c_v33>(), v);

            // Constraint Damping variables
            VarsTools::define_enum_mapping(
                mapping_function, GRInterval<c_Xspatial1, c_Xspatial3>(),
                Xspatial);
            VarsTools::define_enum_mapping(mapping_function, c_Xhat, Xhat);
        }
    };

    //! Structure containing the rhs variables for the matter fields requiring
    //!  2nd derivs
    template <class data_t> struct Diff2Vars
    {
        // data_t phi;
        Tensor<2, data_t> fspatial; // Spatial component of the tensor field
        Tensor<1, data_t> fbar; // Half-projected component of the tensor field
        data_t fhat;            // Scalar part of the tensor field

        /// Defines the mapping between members of Vars and Chombo grid
        ///  variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_symmetric_enum_mapping(
                mapping_function, GRInterval<c_fspatial11, c_fspatial33>(),
                fspatial);
        }
    };

    // Struct for the non grid ADM vars
    template <class data_t>
    using MetricVars = typename ADMFixedBGVars::template Vars<data_t>;

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, including the potential
    template <class data_t, template <typename> class vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &vars, //!< the value of the variables
        const MetricVars<data_t>
            &metric_vars, //!< the value of the metric variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const Tensor<2, data_t>
            &gamma_UU, //!< the inverse metric (raised indices)
        const Tensor<3, data_t> &chris_phys_ULL)
        const; //!< the physical christoffel symbol

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, excluding the potential
    template <class data_t, template <typename> class vars_t>
    static void emtensor_excl_potential(
        emtensor_t<data_t> &out,    //!< the em tensor output
        const vars_t<data_t> &vars, //!< the value of the variables
        const MetricVars<data_t>
            &metric_vars, //!< the value of the metric variables
        const vars_t<Tensor<1, data_t>>
            &d1, //!< the value of the first deriv of phi
        const Tensor<2, data_t>
            &gamma_UU, //!< the inverse metric (raised indices).
        const Tensor<3, data_t>
            &chris_phys_ULL); //!< the physical christoffel symbol

    //! The function which adds in the RHS for the matter field vars,
    //! including the potential
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void matter_rhs(
        rhs_vars_t<data_t> &total_rhs, //!< value of the RHS for all vars
        const vars_t<data_t> &vars,    //!< value of the variables
        const MetricVars<data_t>
            &metric_vars, //!< the value of the metric variables
        const vars_t<Tensor<1, data_t>> &d1,       //!< value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, //!< value of the 2nd derivs
        const vars_t<data_t> &advec)
        const; //!< the value of the advection terms

    //! The function which calculates the RHS for the matter field vars
    //! excluding the potential
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void matter_rhs_excl_potential(
        rhs_vars_t<data_t> &rhs, //!< the value of the RHS terms for the sf vars
        const vars_t<data_t> &vars, //!< the values of all the variables
        const MetricVars<data_t>
            &metric_vars, //!< the value of the metric variables
        const vars_t<Tensor<1, data_t>> &d1,       //!< value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, //!< value of the 2nd derivs
        const vars_t<data_t> &advec) const;
};

#include "FixedBGTensorField.impl.hpp"

#endif /* FIXEDBGTENSORFIELD_HPP_ */
