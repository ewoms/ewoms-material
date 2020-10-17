// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the eWoms project.

  eWoms is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  eWoms is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with eWoms.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 * \copydoc Ewoms::ComputeFromReferencePhase
 */
#ifndef EWOMS_COMPUTE_FROM_REFERENCE_PHASE_HH
#define EWOMS_COMPUTE_FROM_REFERENCE_PHASE_HH

#include <ewoms/material/constraintsolvers/compositionfromfugacities.hh>

#include <ewoms/common/valgrind.hh>

#include <dune/common/fvector.hh>

namespace Ewoms {

/*!
 * \brief Computes all quantities of a generic fluid state if a
 *        reference phase has been specified.
 *
 * This makes it is possible to specify just one phase and let the
 * remaining ones be calculated by the constraint solver. This
 * constraint solver assumes thermodynamic equilibrium. It assumes the
 * following quantities to be set:
 *
 * - composition (mole+mass fractions) of the *reference* phase
 * - temperature of the *reference* phase
 * - saturations of *all* phases
 * - pressures of *all* phases
 *
 * after calling the solve() method the following quantities are
 * calculated in addition:
 *
 * - temperature of *all* phases
 * - density, molar density, molar volume of *all* phases
 * - composition in mole and mass fractions and molaries of *all* phases
 * - mean molar masses of *all* phases
 * - fugacity coefficients of *all* components in *all* phases
 * - if the setViscosity parameter is true, also dynamic viscosities of *all* phases
 * - if the setEnthalpy parameter is true, also specific enthalpies and internal energies of *all* phases
 */
template <class Scalar, class FluidSystem, class Evaluation = Scalar>
class ComputeFromReferencePhase
{
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    typedef Ewoms::CompositionFromFugacities<Scalar, FluidSystem, Evaluation> CompositionFromFugacities;
    typedef Dune::FieldVector<Evaluation, numComponents> ComponentVector;

public:
    /*!
     * \brief Computes all quantities of a generic fluid state if a
     *        reference phase has been specified.
     *
     * This makes it is possible to specify just one phase and let the
     * remaining ones be calculated by the constraint solver. This
     * constraint solver assumes thermodynamic equilibrium. It assumes the
     * following quantities to be set:
     *
     * - composition (mole+mass fractions) of the *reference* phase
     * - temperature of the *all* phases
     * - saturations of *all* phases
     * - pressures of *all* phases
     *
     * after calling the solve() method the following quantities are
     * calculated in addition:
     *
     * - temperature of *all* phases
     * - density, molar density, molar volume of *all* phases
     * - composition in mole and mass fractions and molaries of *all* phases
     * - mean molar masses of *all* phases
     * - fugacity coefficients of *all* components in *all* phases
     * - if the setViscosity parameter is true, also dynamic viscosities of *all* phases
     * - if the setEnthalpy parameter is true, also specific enthalpies and internal energies of *all* phases
     *
     * \param fluidState Thermodynamic state of the fluids
     * \param paramCache  Container for cache parameters
     * \param refPhaseIdx The phase index of the reference phase
     * \param setViscosity Specify whether the dynamic viscosity of
     *                     each phase should also be set.
     * \param setEnthalpy Specify whether the specific
     *                    enthalpy/internal energy of each phase
     *                    should also be set.
     */
    template <class FluidState>
    static void solve(FluidState& fluidState,
                      typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                      unsigned refPhaseIdx,
                      bool setViscosity,
                      bool setEnthalpy)
    {
        // compute the density and enthalpy of the
        // reference phase
        paramCache.updatePhase(fluidState, refPhaseIdx);
        fluidState.setDensity(refPhaseIdx,
                              FluidSystem::density(fluidState,
                                                   paramCache,
                                                   refPhaseIdx));

        if (setEnthalpy)
            fluidState.setEnthalpy(refPhaseIdx,
                                   FluidSystem::enthalpy(fluidState,
                                                         paramCache,
                                                         refPhaseIdx));

        if (setViscosity)
            fluidState.setViscosity(refPhaseIdx,
                                    FluidSystem::viscosity(fluidState,
                                                           paramCache,
                                                           refPhaseIdx));

        // compute the fugacities of all components in the reference phase
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            fluidState.setFugacityCoefficient(refPhaseIdx,
                                              compIdx,
                                              FluidSystem::fugacityCoefficient(fluidState,
                                                                               paramCache,
                                                                               refPhaseIdx,
                                                                               compIdx));
        }

        // compute all quantities for the non-reference phases
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (phaseIdx == refPhaseIdx)
                continue; // reference phase is already calculated

            ComponentVector fugVec;
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                const auto& fug = fluidState.fugacity(refPhaseIdx, compIdx);
                fugVec[compIdx] = Ewoms::decay<Evaluation>(fug);
            }

            CompositionFromFugacities::solve(fluidState, paramCache, phaseIdx, fugVec);

            if (setViscosity)
                fluidState.setViscosity(phaseIdx,
                                        FluidSystem::viscosity(fluidState,
                                                               paramCache,
                                                               phaseIdx));

            if (setEnthalpy)
                fluidState.setEnthalpy(phaseIdx,
                                       FluidSystem::enthalpy(fluidState,
                                                             paramCache,
                                                             phaseIdx));
        }
    }
};

} // namespace Ewoms

#endif
