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
 * \copydoc Ewoms::NullSolidEnergyLaw
 */
#ifndef EWOMS_NULL_SOLID_ENERGY_LAW_HH
#define EWOMS_NULL_SOLID_ENERGY_LAW_HH

#include <ewoms/common/densead/math.hh>

namespace Ewoms
{
/*!
 * \ingroup material
 *
 * \brief Implements a solid energy storage law which just returns 0.
 */
template <class ScalarT>
class NullSolidEnergyLaw
{
public:
    typedef int Params;
    typedef ScalarT Scalar;

    /*!
     * \brief Given a fluid state, compute the volumetric internal energy of the solid
     *        matrix [W/m^3].
     *
     * This solid energy law simply returns 0.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation solidInternalEnergy(const Params& params EWOMS_UNUSED, const FluidState& fluidState EWOMS_UNUSED)
    { return 0.0; }
};
} // namespace Ewoms

#endif
