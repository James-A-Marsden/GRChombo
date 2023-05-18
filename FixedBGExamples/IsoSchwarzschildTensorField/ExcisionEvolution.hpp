/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXCISIONEVOLUTION_HPP_
#define EXCISIONEVOLUTION_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Does excision for fixed BG BH solutions
//! Note that it is does not using simd so one must set disable_simd()
template <class matter_t, class background_t> class ExcisionEvolution
{
    // Use matter_t class
    using Vars = typename matter_t::template Vars<double>;
    // template <class data_t> using Vars = typename matter_t::template
    // Vars<data_t>;

  protected:
    const double m_dx;                              //!< The grid spacing
    const std::array<double, CH_SPACEDIM> m_center; //!< The BH center
    const FourthOrderDerivatives m_deriv;
    const background_t m_background;

  public:
    ExcisionEvolution(const double a_dx,
                      const std::array<double, CH_SPACEDIM> a_center,
                      background_t a_background)
        : m_dx(a_dx), m_deriv(m_dx), m_center(a_center),
          m_background(a_background)
    {
    }

    void compute(const Cell<double> current_cell) const
    // template <class data_t> void compute(const Cell<data_t> current_cell)
    // const
    {
        const Coordinates<double> coords(current_cell, m_dx, m_center);
        double horizon_distance = m_background.excise(current_cell);

        // if (horizon_distance < 0.5)
        if (horizon_distance < 0.95)
        {
            // the matter rhs vars within the excision zone
            // recalculate them - for now set to decay to zero
            Vars vars;

            // const auto d1 = m_deriv.template diff1<Vars>(current_cell);
            VarsTools::assign(vars, 0.0);
            // VarsTools::assign(d1, 0.0);

            /// EXCISE d1

            // assign values of rhs or vars in output box
            current_cell.store_vars(vars);
            // current_cell.store_vars(d1);
        } // else do nothing
        if (horizon_distance < 1.0)
        {
            // the matter rhs vars within the excision zone
            // recalculate them - for now set to decay to zero
            Vars vars;

            // const auto d1 = m_deriv.template diff1<Vars>(current_cell);

            /// EXCISE the damping variables inside horizon
            current_cell.store_vars(0.0, c_Xhat);
            current_cell.store_vars(0.0, c_Xspatial1);
            current_cell.store_vars(0.0, c_Xspatial2);
            current_cell.store_vars(0.0, c_Xspatial3);

            // assign values of rhs or vars in output box
            // current_cell.store_vars(vars);
            // current_cell.store_vars(d1);
        } // else do nothing
    }
};

#endif /* EXCISIONEVOLUTION_HPP_ */
