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
 * \copydoc Ewoms::FluidThermalConduction
 */
#ifndef EWOMS_FLUID_THERMAL_CONDUCTION_LAW_HH
#define EWOMS_FLUID_THERMAL_CONDUCTION_LAW_HH

#include "fluidthermalconductionlawparams.hh"

#include <ewoms/common/spline.hh>

#include <algorithm>

namespace Ewoms {
/*!
 * \ingroup material
 *
 * \brief Implements a thermal conduction law which just takes the conductivity of a given fluid phase.
 */
template <class FluidSystem,
          class ScalarT,
          int phaseIdx,
          class ParamsT = FluidThermalConductionLawParams<ScalarT> >
class FluidThermalConductionLaw
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief Given a fluid state, return the effective thermal conductivity [W/m^2 / (K/m)] of the porous
     *        medium.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation thermalConductivity(const Params& params EWOMS_UNUSED,
                                          const FluidState& fluidState)
    {
        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        paramCache.updatePhase(fluidState, phaseIdx);
        return FluidSystem::template thermalConductivity<FluidState, Evaluation>(fluidState,
                                                                                 paramCache,
                                                                                 phaseIdx);
    }
};
} // namespace Ewoms

#endif
