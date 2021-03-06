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
 * \copydoc Ewoms::LiquidPhase
 */
#ifndef EWOMS_LIQUID_PHASE_HH
#define EWOMS_LIQUID_PHASE_HH

namespace Ewoms {

/*!
 * \ingroup Fluidsystems
 * \brief Represents the liquid phase of a single (pseudo-) component.
 */
template <class Scalar, class ComponentT>
class LiquidPhase
{
public:
    //! \copydoc GasPhase::Component
    typedef ComponentT Component;

    //! \copydoc GasPhase::name
    static const char* name()
    { return Component::name(); }

    //! \copydoc GasPhase::isLiquid
    static bool isLiquid()
    { return true; }

    //! \copydoc GasPhase::isCompressible
    static bool isCompressible()
    { return Component::liquidIsCompressible(); }

    //! \copydoc GasPhase::isIdealGas
    static bool isIdealGas()
    { return false; /* we're a liquid! */ }

    //! \copydoc GasPhase::molarMass
    static Scalar molarMass()
    {  return Component::molarMass(); }

    //! \copydoc GasPhase::criticalTemperature
    static Scalar criticalTemperature()
    {  return Component::criticalTemperature(); }

    //! \copydoc GasPhase::criticalPressure
    static Scalar criticalPressure()
    {  return Component::criticalPressure(); }

    //! \copydoc GasPhase::tripleTemperature
    static Scalar tripleTemperature()
    {  return Component::tripleTemperature(); }

    //! \copydoc GasPhase::triplePressure
    static Scalar triplePressure()
    { return Component::triplePressure(); }

    //! \copydoc GasPhase::vaporPressure
    template <class Evaluation>
    static Evaluation vaporPressure(const Evaluation& temperature)
    {  return Component::vaporPressure(temperature); }

    //! \copydoc GasPhase::density
    template <class Evaluation>
    static Evaluation density(const Evaluation& temperature, const Evaluation& pressure)
    {  return Component::liquidDensity(temperature, pressure); }

    //! \copydoc GasPhase::pressure
    template <class Evaluation>
    static Evaluation pressure(const Evaluation& temperature, const Evaluation& density)
    {  return Component::liquidPressure(temperature, density); }

    //! \copydoc GasPhase::enthalpy
    template <class Evaluation>
    static const Evaluation enthalpy(const Evaluation& temperature, const Evaluation& pressure)
    {  return Component::liquidEnthalpy(temperature, pressure); }

    //! \copydoc GasPhase::internalEnergy
    template <class Evaluation>
    static const Evaluation internalEnergy(const Evaluation& temperature, const Evaluation& pressure)
    { return Component::liquidInternalEnergy(temperature, pressure); }

    //! \copydoc GasPhase::viscosity
    template <class Evaluation>
    static Evaluation viscosity(const Evaluation& temperature, const Evaluation& pressure)
    {  return Component::liquidViscosity(temperature, pressure); }

    //! \copydoc GasPhase::thermalConductivity
    template <class Evaluation>
    static Evaluation thermalConductivity(const Evaluation& temperature, const Evaluation& pressure)
    { return Component::liquidThermalConductivity(temperature, pressure); }

    //! \copydoc GasPhase::heatCapacity
    template <class Evaluation>
    static Evaluation heatCapacity(const Evaluation& temperature, const Evaluation& pressure)
    { return Component::liquidHeatCapacity(temperature, pressure); }
};
} // namespace Ewoms

#endif
