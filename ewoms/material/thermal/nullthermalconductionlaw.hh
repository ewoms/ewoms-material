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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Ewoms::NullThermalConductionLaw
 */
#ifndef EWOMS_NULL_THERMAL_CONDUCTION_LAW_HH
#define EWOMS_NULL_THERMAL_CONDUCTION_LAW_HH

#include <ewoms/common/unused.hh>
#include <ewoms/common/exceptions.hh>

namespace Ewoms {

/*!
 * \ingroup material
 *
 * \brief Implements a dummy law for thermal conduction to which isothermal models
 *        can fall back to.
 *
 * This law just returns 0 unconditionally.
 */
template <class ScalarT>
class NullThermalConductionLaw
{
public:
    typedef int Params;
    typedef ScalarT Scalar;

    /*!
     * \brief Given a fluid state, return the effective thermal conductivity [W/m^2 / (K/m)] of the porous
     *        medium.
     *
     * If this method is called an exception is thrown at run time.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation thermalConductivity(const Params& params EWOMS_UNUSED,
                                          const FluidState& fluidState EWOMS_UNUSED)
    { return 0.0; }
};
} // namespace Ewoms

#endif
