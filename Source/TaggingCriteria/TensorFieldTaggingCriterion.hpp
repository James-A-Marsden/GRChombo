/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TENSORFIELDTAGGINGCRITERION_HPP_
#define TENSORFIELDTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"

class TensorFieldTaggingCriterion
{
  protected:
    const double m_dx;
    const double m_L;
    const int m_level;
    const std::array<double, CH_SPACEDIM> m_center;

    const FourthOrderDerivatives m_deriv;

  public:
    TensorFieldTaggingCriterion(const double dx, const int a_level,
                                const double a_L,
                                const std::array<double, CH_SPACEDIM> a_center)
        : m_dx(dx), m_L(a_L), m_level(a_level), m_center(a_center),
          m_deriv(dx){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        data_t criterion = 0.0;
        data_t grid_criterion = 0.0;
        data_t field_criterion = 0.0;
        // make sure the inner part is regridded around the horizon
        // take L as the length of full grid, so tag inner 1/2
        // of it, which means inner \pm L/4
        // Regridding based on the fixed background
        double ratio = pow(2.0, -(m_level + 2.0));
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        const data_t max_abs_xy = simd_max(abs(coords.x), abs(coords.y));
        const data_t max_abs_xyz = simd_max(max_abs_xy, abs(coords.z));
        auto regrid = simd_compare_lt(max_abs_xyz, m_L * ratio);
        grid_criterion = simd_conditional(regrid, 100.0, criterion);

        // Regridding based on the tensor field and momentum
        // looks awful but works using GRChombo's current functions

        Tensor<1, data_t> d1_fspatial11;
        Tensor<1, data_t> d1_fspatial22;
        Tensor<1, data_t> d1_fspatial33;
        Tensor<1, data_t> d1_fspatial12;
        Tensor<1, data_t> d1_fspatial13;
        Tensor<1, data_t> d1_fspatial23;

        Tensor<1, data_t> d1_v11;
        Tensor<1, data_t> d1_v22;
        Tensor<1, data_t> d1_v33;
        Tensor<1, data_t> d1_v12;
        Tensor<1, data_t> d1_v13;
        Tensor<1, data_t> d1_v23;

        FOR(idir)
        m_deriv.diff1(d1_fspatial11, current_cell, idir, c_fspatial11);
        FOR(idir)
        m_deriv.diff1(d1_fspatial22, current_cell, idir, c_fspatial22);
        /*
        FOR(idir) m_deriv.diff1(d1_fspatial33, current_cell, idir,
        c_fspatial33); FOR(idir) m_deriv.diff1(d1_fspatial12, current_cell,
        idir, c_fspatial12); FOR(idir) m_deriv.diff1(d1_fspatial13,
        current_cell, idir, c_fspatial13); FOR(idir)
        m_deriv.diff1(d1_fspatial23, current_cell, idir, c_fspatial23);
        */

        FOR(idir) m_deriv.diff1(d1_v11, current_cell, idir, c_v11);
        FOR(idir) m_deriv.diff1(d1_v22, current_cell, idir, c_v22);
        /*
        FOR(idir) m_deriv.diff1(d1_v33, current_cell, idir, c_v33);
        FOR(idir) m_deriv.diff1(d1_v12, current_cell, idir, c_v12);
        FOR(idir) m_deriv.diff1(d1_v13, current_cell, idir, c_v13);
        FOR(idir) m_deriv.diff1(d1_v23, current_cell, idir, c_v23);
        */

        data_t mod_d1_fspatial = 0.0;
        data_t mod_d1_v = 0.0;

        FOR(idir)
        {

            mod_d1_fspatial +=
                d1_fspatial11[idir] * d1_fspatial11[idir] +
                d1_fspatial22[idir] *
                    d1_fspatial22
                        [idir]; // + d1_fspatial33[idir] * d1_fspatial33[idir]
                                //+d1_fspatial12[idir] * d1_fspatial12[idir] +
                                //d1_fspatial13[idir] * d1_fspatial13[idir] +
                                //d1_fspatial23[idir] * d1_fspatial23[idir];

            mod_d1_v +=
                d1_v11[idir] * d1_v11[idir] +
                d1_v22[idir] * d1_v22[idir]; // + d1_v33[idir] * d1_v33[idir]
                                             //+d1_v12[idir] * d1_v12[idir] +
                                             //d1_v13[idir] * d1_v13[idir] +
                                             //d1_v23[idir] * d1_v23[idir];
        }
        // temp
        // const double threshold_field = 0.0000001;
        const double threshold_field = 0.05;

        const data_t field_regrid =
            m_dx * (sqrt(mod_d1_fspatial / 2.0) + sqrt(mod_d1_v / 2.0)) /
            threshold_field;

        // Take the larger of the two regridding criteria

        // Write back into the flattened Chombo box
        // Just consider the tensor field values?

        // Only regrid near the middle of the box to save computing time

        auto regrid_centre = simd_compare_gt(max_abs_xyz, 32.0);
        field_criterion = simd_conditional(regrid_centre, 0.0, field_regrid);
        criterion = simd_max(grid_criterion, field_criterion);
        current_cell.store_vars(grid_criterion, 0.0);
    }
};

#endif /* TENSORFIELDTAGGINGCRITERION_HPP_ */
