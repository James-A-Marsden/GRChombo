/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOOSTEDPUNCTURETRACKERTAGGINGCRITERION_HPP_
#define BOOSTEDPUNCTURETRACKERTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SphericalExtraction.hpp"
#include "Tensor.hpp"
#include "VarsTools.hpp"

//! This class tags cells based on three criteria - the
//! value of the second derivs, the extraction regions
//! and the puncture horizons (which must be covered to
//! a given level
template <class derivative_t = FourthOrderDerivatives>
class BoostedPunctureTrackerTaggingCriterion
{
  protected:
    const double m_dx;
    const int m_level;
    const std::array<int, 2> m_puncture_max_levels;
    const bool m_activate_extraction;
    const bool m_track_punctures;
    const std::array<double, 2> m_puncture_masses;
    const std::array<double, 2> m_puncture_momentum;
    const SphericalExtraction::params_t m_params;
    const derivative_t m_deriv;
    const std::vector<std::array<double, CH_SPACEDIM>> m_puncture_coords;
    const double m_buffer;
    const double m_puncture_min_separation; // the minimum coordinate separation
                                            // to tag around the punctures
                                            // separately

  public:
    template <class data_t> struct Vars
    {
        data_t chi; //!< Conformal factor

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_chi, chi);
        }
    };

    // The constructor
    BoostedPunctureTrackerTaggingCriterion(
        const double dx, const int a_level,
        const std::array<int, 2> a_puncture_max_levels,
        const SphericalExtraction::params_t a_params,
        const std::vector<std::array<double, CH_SPACEDIM>> a_puncture_coords,
        const bool activate_extraction = false,
        const bool track_punctures = false,
        const std::array<double, 2> a_puncture_masses = {1.0, 1.0},
        const std::array<double, 2> a_puncture_momentum = {1.0, 1.0},
        const double a_buffer = 0.5,
        const double a_puncture_min_separation = 1e-3)
        : m_dx(dx), m_deriv(dx), m_params(a_params), m_level(a_level),
          m_puncture_max_levels(a_puncture_max_levels),
          m_track_punctures(track_punctures),
          m_activate_extraction(activate_extraction),
          m_puncture_masses(a_puncture_masses),
          m_puncture_momentum(a_puncture_momentum),
          m_puncture_coords(a_puncture_coords), m_buffer(a_buffer),
          m_puncture_min_separation(a_puncture_min_separation)
    {
        // check that the number of punctures is consistent
        CH_assert(m_puncture_coords.size() == 2);
    };

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // first test the gradients for regions of high curvature
        const auto d2 = m_deriv.template diff2<Vars>(current_cell);
        data_t mod_d2_chi = 0;
        FOR(idir, jdir)
        {
            mod_d2_chi += d2.chi[idir][jdir] * d2.chi[idir][jdir];
        }
        data_t criterion = m_dx * sqrt(mod_d2_chi);

        // if extracting weyl data at a given radius, enforce a given resolution
        // there
        if (m_activate_extraction)
        {
            for (int iradius = 0; iradius < m_params.num_extraction_radii;
                 ++iradius)
            {
                // regrid if within extraction level and not at required
                // refinement
                if (m_level < m_params.extraction_levels[iradius])
                {
                    const Coordinates<data_t> coords(current_cell, m_dx,
                                                     m_params.center);
                    const data_t r = coords.get_radius();
                    // add a 20% buffer to extraction zone so not too near to
                    // boundary
                    auto regrid = simd_compare_lt(
                        r, 1.2 * m_params.extraction_radii[iradius]);
                    criterion = simd_conditional(regrid, 100.0, criterion);
                }
            }
        }

        // next use the punctures to decide where to tag
        if (m_track_punctures)
        {
            double puncture_separation2 = 0.0;
            double factor = 0.0;
            FOR(idir)
            {
                double displacement =
                    m_puncture_coords[0][idir] - m_puncture_coords[1][idir];
                puncture_separation2 += displacement * displacement;
            }
            double puncture_separation = sqrt(puncture_separation2);

            const int merger_max_level =
                min(m_puncture_max_levels[0], m_puncture_max_levels[1]);

            if (puncture_separation > m_puncture_min_separation)
            {
                // punctures still far enough apart so tag around each one
                // separately
                for (int ipuncture = 0; ipuncture < 2; ++ipuncture)
                {

                    if ((m_level <= m_puncture_max_levels[ipuncture]))
                    {
                        // we want the 2nd and 3rd levels above
                        // puncture_max_level to be twice the size of the next
                        // finest level
                        const Coordinates<data_t> coords(
                            current_cell, m_dx, m_puncture_coords[ipuncture]);
                        const double v1 = m_puncture_momentum[ipuncture];
                        const double gamma1 = 1. / sqrt(1. - (v1 * v1));
                        const double xfactor =
                            m_puncture_masses[ipuncture] * gamma1 *
                            pow(1.3,
                                m_puncture_max_levels[ipuncture] - m_level);
                        const double yfactor =
                            m_puncture_masses[ipuncture] * gamma1 *
                            pow(1.3,
                                m_puncture_max_levels[ipuncture] - m_level);
                        const data_t max_abs_xy = simd_max(
                            abs(coords.x / xfactor), abs(coords.y / yfactor));
                        auto regrid = simd_compare_lt(max_abs_xy, 1.);
                        criterion = simd_conditional(regrid, 100.0, criterion);
                    }
                }
            }
            else
            {
                if (m_level > merger_max_level)
                {
                    // drop any finer levels after merger
                    // tagging for merger BH handled below
                    criterion = 0.0;
                }
            }

            double sum_masses = m_puncture_masses[0] + m_puncture_masses[1];
            // if punctures are close enough together tag cells at the
            // center of mass for the merger BH
            if ((puncture_separation < sum_masses + m_buffer) &&
                (m_level <= merger_max_level))
            {
                std::array<double, CH_SPACEDIM> center_of_mass;
                FOR(idir)
                {
                    center_of_mass[idir] =
                        (m_puncture_masses[0] * m_puncture_coords[0][idir] +
                         m_puncture_masses[1] * m_puncture_coords[1][idir]) /
                        sum_masses;
                }
                Coordinates<data_t> coords(current_cell, m_dx, center_of_mass);
                const data_t max_abs_xy =
                    simd_max(abs(coords.x), abs(coords.y));
#if CH_SPACEDIM == 3
                const data_t max_abs_xyz = simd_max(max_abs_xy, abs(coords.z));
#endif

                const double factor =
                    pow(2.0, min(merger_max_level - m_level - 1, 2));
#if CH_SPACEDIM == 2
                auto regrid2 = simd_compare_lt(
                    max_abs_xy, factor * (sum_masses + m_buffer));
#endif

#if CH_SPACEDIM == 3
                auto regrid2 = simd_compare_lt(
                    max_abs_xyz, factor * (sum_masses + m_buffer));
#endif
                criterion = simd_conditional(regrid2, 100.0, criterion);
            }
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* PUNCTURETRACKERTAGGINGCRITERION_HPP_ */
