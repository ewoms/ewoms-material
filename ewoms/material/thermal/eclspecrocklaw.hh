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
 * \copydoc Ewoms::EclSpecrockLaw
 */
#ifndef EWOMS_ECL_SPECROCK_LAW_HH
#define EWOMS_ECL_SPECROCK_LAW_HH

#include "eclspecrocklawparams.hh"

#include <ewoms/common/densead/math.hh>

namespace Ewoms
{
/*!
 * \ingroup material
 *
 * \brief Implements the volumetric interior energy relations of rock used by ECL.
 *
 * This class uses the approach defined via SPECROCK keyword.
 */
template <class ScalarT,
          class ParamsT = EclSpecrockLawParams<ScalarT> >
class EclSpecrockLaw
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief Given a fluid state, compute the volumetric internal energy of the rock [W/m^3].
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation solidInternalEnergy(const Params& params, const FluidState& fluidState)
    {
        const auto& T = fluidState.temperature(/*phaseIdx=*/0);
        return params.internalEnergyFunction().eval(T, /*extrapolate=*/true);
    }
};
} // namespace Ewoms

#endif
