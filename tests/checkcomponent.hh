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
 * \copydoc checkComponent
 */
#ifndef EWOMS_CHECK_COMPONENT_HH
#define EWOMS_CHECK_COMPONENT_HH

#include <ewoms/common/unused.hh>
#include <dune/common/classname.hh>

#include <iostream>
#include <string>

/*!
 * \brief Ensures that a class which represents a chemical components adheres to the
 *        components API.
 *
 * Note that this does *not* imply that the methods are implemented or even make sense...
 */
template <class Component, class Evaluation>
void checkComponent()
{
    std::cout << "Testing component '" << Dune::className<Component>() << "'\n";

    // make sure the necessary typedefs exist
    typedef typename Component::Scalar Scalar;

    // make sure the necessary constants are exported
    bool isTabulated EWOMS_UNUSED = Component::isTabulated;

    // test for the gas-phase functions
    Evaluation T=0, p=0;
    while (0) {
        { bool b EWOMS_UNUSED = Component::gasIsCompressible(); }
        { bool b EWOMS_UNUSED = Component::gasIsIdeal(); }
        { bool b EWOMS_UNUSED = Component::liquidIsCompressible(); }
        { std::string s EWOMS_UNUSED = Component::name(); }
        { Scalar M EWOMS_UNUSED = Component::molarMass(); }
        { Scalar Tc EWOMS_UNUSED = Component::criticalTemperature(); }
        { Scalar pc EWOMS_UNUSED = Component::criticalPressure(); }
        { Scalar Tt EWOMS_UNUSED = Component::tripleTemperature(); }
        { Scalar pt EWOMS_UNUSED = Component::triplePressure(); }
        { Evaluation pv EWOMS_UNUSED = Component::vaporPressure(T); }
        { Evaluation rho EWOMS_UNUSED = Component::gasDensity(T, p); }
        { Evaluation rho EWOMS_UNUSED = Component::liquidDensity(T, p); }
        { Evaluation h EWOMS_UNUSED = Component::gasEnthalpy(T, p); }
        { Evaluation h EWOMS_UNUSED = Component::liquidEnthalpy(T, p); }
        { Evaluation u EWOMS_UNUSED = Component::gasInternalEnergy(T, p); }
        { Evaluation u EWOMS_UNUSED = Component::liquidInternalEnergy(T, p); }
        { Evaluation mu EWOMS_UNUSED = Component::gasViscosity(T, p); }
        { Evaluation mu EWOMS_UNUSED = Component::liquidViscosity(T, p); }
        { Evaluation lambda EWOMS_UNUSED = Component::gasThermalConductivity(T, p); }
        { Evaluation lambda EWOMS_UNUSED = Component::liquidThermalConductivity(T, p); }
        { Evaluation cp EWOMS_UNUSED = Component::gasHeatCapacity(T, p); }
        { Evaluation cp EWOMS_UNUSED = Component::liquidHeatCapacity(T, p); }
    }
    std::cout << "----------------------------------\n";
}

#endif
