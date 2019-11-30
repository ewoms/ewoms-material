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
 * \copydoc Ewoms::FluidThermalConductionParams
 */
#ifndef EWOMS_FLUID_THERMAL_CONDUCTION_LAW_PARAMS_HH
#define EWOMS_FLUID_THERMAL_CONDUCTION_LAW_PARAMS_HH

namespace Ewoms {
/*!
 * \brief Parameters for the thermal conduction law which just takes the conductivity of a given fluid phase.
 */
template <class ScalarT>
class FluidThermalConductionLawParams
{
    // do not copy!
    FluidThermalConductionLawParams(const FluidThermalConductionLawParams&)
    {}

public:
    typedef ScalarT Scalar;

    FluidThermalConductionLawParams()
    { }

};

} // namespace Ewoms

#endif
