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
 * \copydoc Ewoms::EclThconrLaw
 */
#ifndef EWOMS_ECL_THCONR_LAW_HH
#define EWOMS_ECL_THCONR_LAW_HH

#include "eclthconrlawparams.hh"

#include <ewoms/common/densead/math.hh>

namespace Ewoms
{
/*!
 * \ingroup material
 *
 * \brief Implements the total thermal conductivity relations specified by the ECL THCONR.
 */
template <class ScalarT,
          class FluidSystem,
          class ParamsT = EclThconrLawParams<ScalarT> >
class EclThconrLaw
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief Given a fluid state, return the total thermal conductivity [W/m^2 / (K/m)] of the porous
     *        medium.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation thermalConductivity(const Params& params,
                                          const FluidState& fluidState)
    {
        // THCONR + THCONSF approach.
        Scalar lambdaRef = params.referenceTotalThermalConductivity();
        static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;
        if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
            Scalar alpha = params.dTotalThermalConductivity_dSg();
            const Evaluation& Sg = Ewoms::decay<Evaluation>(fluidState.saturation(gasPhaseIdx));
            return lambdaRef*(1.0 - alpha*Sg);
        } else {
            return lambdaRef;
        }
    }
};
} // namespace Ewoms

#endif
