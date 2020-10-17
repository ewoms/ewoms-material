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
 *
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system assuming
 *        thermodynamic equilibrium.
 */
#ifndef EWOMS_IMMISCIBLE_FLUID_STATE_HH
#define EWOMS_IMMISCIBLE_FLUID_STATE_HH

#include "modularfluidstate.hh"

#include <ewoms/common/valgrind.hh>

#include <algorithm>

#include <string.h>

namespace Ewoms {

/*!
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system assuming
 *        thermodynamic equilibrium.
 */
template <class Scalar, class FluidSystem, bool storeEnthalpy=true>
class ImmiscibleFluidState;

// specialization for the enthalpy enabled case
template <class Scalar, class FluidSystem>
class ImmiscibleFluidState<Scalar, FluidSystem, true>
    : public ModularFluidState<Scalar,
                               FluidSystem::numPhases,
                               FluidSystem::numComponents,
                               FluidStateExplicitPressureModule<Scalar, FluidSystem::numPhases, ImmiscibleFluidState<Scalar, FluidSystem, true> >,
                               FluidStateEquilibriumTemperatureModule<Scalar, FluidSystem::numPhases, ImmiscibleFluidState<Scalar, FluidSystem, true> >,
                               FluidStateImmiscibleCompositionModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, true> >,
                               FluidStateExplicitFugacityModule<Scalar, FluidSystem::numPhases, FluidSystem::numComponents, ImmiscibleFluidState<Scalar, FluidSystem, true> >,
                               FluidStateExplicitSaturationModule<Scalar, FluidSystem::numPhases, ImmiscibleFluidState<Scalar, FluidSystem, true> >,
                               FluidStateExplicitDensityModule<Scalar, FluidSystem::numPhases, ImmiscibleFluidState<Scalar, FluidSystem, true> >,
                               FluidStateExplicitViscosityModule<Scalar, FluidSystem::numPhases, ImmiscibleFluidState<Scalar, FluidSystem, true> >,
                               FluidStateExplicitEnthalpyModule<Scalar, FluidSystem::numPhases, ImmiscibleFluidState<Scalar, FluidSystem, true> > >
{
public:
    ImmiscibleFluidState()
    {}
};

// specialization for the enthalpy disabled case
template <class Scalar, class FluidSystem>
class ImmiscibleFluidState<Scalar, FluidSystem, false>
    : public ModularFluidState<Scalar,
                               FluidSystem::numPhases,
                               FluidSystem::numComponents,
                               FluidStateExplicitPressureModule<Scalar, FluidSystem::numPhases, ImmiscibleFluidState<Scalar, FluidSystem, false> >,
                               FluidStateEquilibriumTemperatureModule<Scalar, FluidSystem::numPhases, ImmiscibleFluidState<Scalar, FluidSystem, false> >,
                               FluidStateImmiscibleCompositionModule<Scalar, FluidSystem, ImmiscibleFluidState<Scalar, FluidSystem, false> >,
                               FluidStateExplicitFugacityModule<Scalar, FluidSystem::numPhases, FluidSystem::numComponents, ImmiscibleFluidState<Scalar, FluidSystem, false> >,
                               FluidStateExplicitSaturationModule<Scalar, FluidSystem::numPhases, ImmiscibleFluidState<Scalar, FluidSystem, false> >,
                               FluidStateExplicitDensityModule<Scalar, FluidSystem::numPhases, ImmiscibleFluidState<Scalar, FluidSystem, false> >,
                               FluidStateExplicitViscosityModule<Scalar, FluidSystem::numPhases, ImmiscibleFluidState<Scalar, FluidSystem, false> >,
                               FluidStateNullEnthalpyModule<Scalar, FluidSystem::numPhases, ImmiscibleFluidState<Scalar, FluidSystem, false> > >
{
public:
    ImmiscibleFluidState()
    {}
};
} // namespace Ewoms

#endif
