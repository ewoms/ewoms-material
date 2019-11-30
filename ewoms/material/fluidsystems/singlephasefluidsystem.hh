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
 * \copydoc Ewoms::SinglePhaseFluidSystem
 */
#ifndef EWOMS_SINGLE_PHASE_FLUIDSYSTEM_HH
#define EWOMS_SINGLE_PHASE_FLUIDSYSTEM_HH

#include "basefluidsystem.hh"
#include "nullparametercache.hh"

#include <ewoms/material/fluidsystems/liquidphase.hh>
#include <ewoms/material/fluidsystems/gasphase.hh>
#include <ewoms/material/components/simpleh2o.hh>
#include <ewoms/material/components/h2o.hh>
#include <ewoms/material/components/n2.hh>
#include <ewoms/material/components/tabulatedcomponent.hh>

#include <ewoms/common/unused.hh>

#include <limits>
#include <cassert>

namespace Ewoms {

/*!
 * \ingroup Fluidsystems
 *
 * \brief A fluid system for single phase models.
 *
 * The fluid is defined as a template parameter. For existing
 * components the Ewoms::LiquidPhase<Component> and
 * Ewoms::GasPhase<Component> may be used.
 */
template <class Scalar, class Fluid>
class SinglePhaseFluidSystem
    : public BaseFluidSystem<Scalar, SinglePhaseFluidSystem<Scalar, Fluid> >
{
    typedef SinglePhaseFluidSystem<Scalar, Fluid> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;

public:
    //! \copydoc BaseFluidSystem::ParameterCache
    template <class Evaluation>
    struct ParameterCache : public Ewoms::NullParameterCache<Evaluation>
    {};

    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 1;

    //! \copydoc BaseFluidSystem::phaseName
    static const char* phaseName(unsigned phaseIdx EWOMS_OPTIM_UNUSED)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return Fluid::name();
    }

    //! \copydoc BaseFluidSystem::isLiquid
    static bool isLiquid(unsigned /*phaseIdx*/)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        return Fluid::isLiquid();
    }

    //! \copydoc BaseFluidSystem::isCompressible
    static bool isCompressible(unsigned /*phaseIdx*/)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        // let the fluid decide
        return Fluid::isCompressible();
    }

    //! \copydoc BaseFluidSystem::isIdealGas
    static bool isIdealMixture(unsigned /*phaseIdx*/)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        // we assume immisibility
        return true;
    }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealGas(unsigned /*phaseIdx*/)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        // let the fluid decide
        return Fluid::isIdealGas();
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 1;

    //! \copydoc BaseFluidSystem::componentName
    static const char* componentName(unsigned compIdx EWOMS_OPTIM_UNUSED)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        return Fluid::name();
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(unsigned /*compIdx*/)
    {
        //assert(0 <= compIdx && compIdx < numComponents);

        return Fluid::molarMass();
    }

    /*!
     * \brief Critical temperature of a component [K].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalTemperature(unsigned /*compIdx*/)
    {
        //assert(0 <= compIdx && compIdx < numComponents);

        return Fluid::criticalTemperature();
    }

    /*!
     * \brief Critical pressure of a component [Pa].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalPressure(unsigned /*compIdx*/)
    {
        //assert(0 <= compIdx && compIdx < numComponents);

        return Fluid::criticalPressure();
    }

    /*!
     * \brief The acentric factor of a component [].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar acentricFactor(unsigned /*compIdx*/)
    {
        //assert(0 <= compIdx && compIdx < numComponents);

        return Fluid::acentricFactor();
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    //! \copydoc BaseFluidSystem::init
    static void init()
    { }

    //! \copydoc BaseFluidSystem::density
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval density(const FluidState& fluidState,
                           const ParameterCache<ParamCacheEval>& /*paramCache*/,
                           unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto& T = Ewoms::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Ewoms::decay<LhsEval>(fluidState.pressure(phaseIdx));
        return Fluid::density(T, p);
    }

    //! \copydoc BaseFluidSystem::viscosity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval viscosity(const FluidState& fluidState,
                             const ParameterCache<ParamCacheEval>& /*paramCache*/,
                             unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto& T = Ewoms::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Ewoms::decay<LhsEval>(fluidState.pressure(phaseIdx));
        return Fluid::viscosity(T, p);
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval fugacityCoefficient(const FluidState& /*fluidState*/,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx,
                                       unsigned compIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIdx && compIdx < numComponents);

        if (phaseIdx == compIdx)
            // TODO (?): calculate the real fugacity coefficient of
            // the component in the fluid. Probably that's not worth
            // the effort, since the fugacity coefficient of the other
            // component is infinite anyway...
            return 1.0;
        return std::numeric_limits<Scalar>::infinity();
    }

    //! \copydoc BaseFluidSystem::enthalpy
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval enthalpy(const FluidState& fluidState,
                            const ParameterCache<ParamCacheEval>& /*paramCache*/,
                            unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto& T = Ewoms::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Ewoms::decay<LhsEval>(fluidState.pressure(phaseIdx));
        return Fluid::enthalpy(T, p);
    }

    //! \copydoc BaseFluidSystem::thermalConductivity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval thermalConductivity(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto& T = Ewoms::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Ewoms::decay<LhsEval>(fluidState.pressure(phaseIdx));
        return Fluid::thermalConductivity(T, p);
    }

    //! \copydoc BaseFluidSystem::heatCapacity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval heatCapacity(const FluidState& fluidState,
                                const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto& T = Ewoms::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Ewoms::decay<LhsEval>(fluidState.pressure(phaseIdx));
        return Fluid::heatCapacity(T, p);
    }
};

} // namespace Ewoms

#endif
