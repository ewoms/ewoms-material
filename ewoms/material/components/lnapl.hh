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
 * \copydoc Ewoms::LNAPL
 */
#ifndef EWOMS_LNAPL_HH
#define EWOMS_LNAPL_HH

#include "component.hh"

#include <ewoms/common/unused.hh>

namespace Ewoms {
/*!
 * \ingroup Components
 *
 * \brief A simple implementation of a LNAPL, e.g. a kind of oil
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class LNAPL : public Component<Scalar, LNAPL<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the iso-octane.
     */
    static const char* name()
    { return "LNAPL"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of iso-octane.
     */
    static Scalar molarMass()
    { return 0.11423; }

    /*!
     * \brief Returns true iff the liquid phase is assumed to be compressible
     */
    static bool liquidIsCompressible()
    { return false; }

    /*!
     * \brief Rough estimate of the density of oil \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidDensity(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
    { return 692.0; }

    /*!
     * \brief Rough estimate of the viscosity of oil in \f$\mathrm{[Pa*s]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidViscosity(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
    { return 0.005; }

    /*!
     * \brief The enthalpy of iso-octane at a given pressure and temperature \f$\mathrm{[J/kg]}\f$.
     *
     * We simply use the value of iso-octane here.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidEnthalpy(const Evaluation& temperature,
                                     const Evaluation& pressure EWOMS_UNUSED)
    {
        return 240.0/molarMass() * temperature; // [J/kg]
    }

    /*!
     * \brief Specific isobaric heat capacity \f$[J/(kg K)]\f$ of liquid iso-octane.
     *
     * We simply use the value of iso-octane here.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidHeatCapacity(const Evaluation& temperature EWOMS_UNUSED,
                                         const Evaluation& pressure EWOMS_UNUSED)
    {
        return 240.0/molarMass();
    }

    /*!
     * \brief Specific heat conductivity of liquid TCE \f$\mathrm{[W/(m K)]}\f$.
     *
     * \todo The value returned here is a guess which does not necessarily correspond to reality in any way!
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidThermalConductivity(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
    {
        return 0.3; // TODO: guess
    }
};

} // namespace Ewoms

#endif
