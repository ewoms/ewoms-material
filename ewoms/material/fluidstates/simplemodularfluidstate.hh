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
 * \copydoc Ewoms::SimpleModularFluidState
 */
#ifndef EWOMS_SIMPLE_MODULAR_FLUID_STATE_HH
#define EWOMS_SIMPLE_MODULAR_FLUID_STATE_HH

#include "fluidstatepressuremodules.hh"
#include "fluidstatetemperaturemodules.hh"
#include "fluidstatecompositionmodules.hh"
#include "fluidstatefugacitymodules.hh"
#include "fluidstatesaturationmodules.hh"
#include "fluidstatedensitymodules.hh"
#include "fluidstateviscositymodules.hh"
#include "fluidstateenthalpymodules.hh"
#include "modularfluidstate.hh"

#include <type_traits>

namespace Ewoms {
// this macro is a small hack to prevent death-through verbosity
#define EWOMS_SMFS SimpleModularFluidState<ScalarT, \
                                         numPhasesV, \
                                         numComponentsV,    \
                                         FluidSystem,       \
                                         storePressure,     \
                                         storeTemperature,  \
                                         storeComposition,  \
                                         storeFugacity,     \
                                         storeSaturation,   \
                                         storeDensity,      \
                                         storeViscosity,    \
                                         storeEnthalpy>

/*!
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system assuming
 *        thermodynamic equilibrium.
 *
 * This class uses simpler and slightly less flexible template parameters as
 * ModularFluidState. Except for this, it is identical.
 */
template <class ScalarT,
          unsigned numPhasesV,
          unsigned numComponentsV,
          class FluidSystem, // only needed if the compositional stuff enabled
          bool storePressure,
          bool storeTemperature,
          bool storeComposition,
          bool storeFugacity,
          bool storeSaturation,
          bool storeDensity,
          bool storeViscosity,
          bool storeEnthalpy>
class SimpleModularFluidState
    : public ModularFluidState<ScalarT, numPhasesV, numComponentsV,
                               typename std::conditional<storePressure,
                                                         FluidStateExplicitPressureModule<ScalarT, numPhasesV, EWOMS_SMFS>,
                                                         FluidStateNullPressureModule<ScalarT> >::type,
                               typename std::conditional<storeTemperature,
                                                         FluidStateExplicitTemperatureModule<ScalarT, numPhasesV, EWOMS_SMFS>,
                                                         FluidStateNullTemperatureModule<ScalarT> >::type,
                               typename std::conditional<storeComposition,
                                                         FluidStateExplicitCompositionModule<ScalarT, FluidSystem, EWOMS_SMFS>,
                                                         FluidStateNullCompositionModule<ScalarT> >::type,
                               typename std::conditional<storeFugacity,
                                                         FluidStateExplicitFugacityModule<ScalarT, numPhasesV, numComponentsV, EWOMS_SMFS>,
                                                         FluidStateNullFugacityModule<ScalarT> >::type,
                               typename std::conditional<storeSaturation,
                                                         FluidStateExplicitSaturationModule<ScalarT, numPhasesV, EWOMS_SMFS>,
                                                         FluidStateNullSaturationModule<ScalarT> >::type,
                               typename std::conditional<storeDensity,
                                                         FluidStateExplicitDensityModule<ScalarT, numPhasesV, EWOMS_SMFS>,
                                                         FluidStateNullDensityModule<ScalarT, numPhasesV, EWOMS_SMFS> >::type,
                               typename std::conditional<storeViscosity,
                                                         FluidStateExplicitViscosityModule<ScalarT, numPhasesV, EWOMS_SMFS>,
                                                         FluidStateNullViscosityModule<ScalarT, numPhasesV, EWOMS_SMFS> >::type,
                               typename std::conditional<storeEnthalpy,
                                                         FluidStateExplicitEnthalpyModule<ScalarT, numPhasesV, EWOMS_SMFS>,
                                                         FluidStateNullEnthalpyModule<ScalarT, numPhasesV, EWOMS_SMFS> >::type
                               >
{};

// we don't need the death-prevention macro anymore
#undef EWOMS_SMFS

} // namespace Ewoms

#endif
